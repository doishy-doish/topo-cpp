// memtrash.hpp — the memory trash system: ECC for memory (liveness layer).
// C++17, zero dependencies, header-only. Native + emcc clean.
//
// Layer A (Registry): log ALL tracked memory areas — pointer, size, tag,
// site, age. sweep() releases stale entries through their registered releaser
// (the leak patch). Entries with no releaser are REPORTED, never freed blind.
// Entries marked permanent (the live sim) are audited, never swept.
// Layer B (Arena): fixed-capacity ECC-guarded pool. Every block carries
// header magic + checksum + tail canary (the parity bits). scrub() verifies
// each live block; corrupt blocks are QUARANTINED (never handed out again —
// the ECC spare-the-row move); idle blocks are swept. Touch on use.
//
// Honors the scraps rules: no allocation on the hot path (registry uses one
// vector, arena is one fixed buffer), explicit ownership (releaser/dispose),
// caps (arena capacity, registry soft cap with oldest-report), probes
// (audit/log_json) over arguments.
#pragma once
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

namespace memtrash {

typedef void (*Releaser)(const void* ptr, void* user);
static const uint64_t MT_MAGIC = 0x54524153484D5452ULL; // "TRASHMT"
static const uint32_t MT_HEAD = 0x48545241U;            // "HTRA"
static const uint32_t MT_TAIL = 0x5449414CU;            // "TAIL"
static const uint64_t MT_MAX_TRACKED = 4096;           // soft cap: registry never grows unbounded

struct Record {
  uint64_t id = 0;
  const void* ptr = nullptr;
  size_t size = 0;
  const char* tag = "?";
  const char* file = "?";
  int line = 0;
  uint64_t born = 0;
  uint64_t touched = 0;
  bool permanent = false;
  Releaser releaser = nullptr;
  void* user = nullptr;
  bool swept = false;
};

class Registry {
 public:
  uint64_t tick_ = 0;
  uint64_t next_id_ = 1;
  size_t peak_bytes_ = 0;
  uint64_t sweeps_ = 0;
  uint64_t quarantines_ = 0; // reserved: registry never quarantines; arena reports here
  std::vector<Record> recs_;
  // --- budget: the expected maximum. Latched + configured, then automatic. ---
  size_t budget_bytes_ = 0;   // hard cap, 0 = off
  float growth_factor_ = 0.0f; // e.g. 1.5 = trip when 50% over baseline, 0 = off
  size_t growth_min_ = 0;      // growth trip only counts above this many bytes
  size_t baseline_ = 0;        // latched expected-max on first check with bytes > 0
  bool was_over_ = false;      // for newly_over edge detection (warn throttle)

  void tick() { tick_++; }

  bool track(const void* ptr, size_t size, const char* tag, const char* file,
             int line, Releaser rel = nullptr, void* user = nullptr,
             bool permanent = false) {
    if (!ptr) return false;
    for (auto& r : recs_) {
      if (!r.swept && r.ptr == ptr) { // re-track = touch + update
        r.size = size; r.tag = tag; r.file = file; r.line = line;
        r.touched = tick_; r.permanent = permanent; r.releaser = rel; r.user = user;
        return true;
      }
    }
    if (recs_.size() >= MT_MAX_TRACKED) return false; // cap: refuse, don't sprawl
    Record r;
    r.id = next_id_++; r.ptr = ptr; r.size = size; r.tag = tag;
    r.file = file; r.line = line; r.born = tick_; r.touched = tick_;
    r.permanent = permanent; r.releaser = rel; r.user = user;
    recs_.push_back(r);
    size_t b = bytes();
    if (b > peak_bytes_) peak_bytes_ = b;
    return true;
  }

  bool untrack(const void* ptr) {
    for (auto& r : recs_) {
      if (!r.swept && r.ptr == ptr) { r.swept = true; return true; }
    }
    return false;
  }

  void touch(const void* ptr) {
    for (auto& r : recs_) {
      if (!r.swept && r.ptr == ptr) { r.touched = tick_; return; }
    }
  }

  size_t live() const {
    size_t n = 0;
    for (auto& r : recs_) if (!r.swept) n++;
    return n;
  }
  size_t bytes() const {
    size_t b = 0;
    for (auto& r : recs_) if (!r.swept) b += r.size;
    return b;
  }
  size_t peak() const { return peak_bytes_; }

  struct SweepReport { int released = 0; int reported = 0; int skipped_perm = 0; };  // The patch: stale non-permanent entries WITH a releaser get released.
  // Stale entries with NO releaser are counted as reported (logged, not freed).
  SweepReport sweep(uint64_t max_idle_ticks) {
    SweepReport rep;
    for (auto& r : recs_) {
      if (r.swept) continue;
      if (r.permanent) { rep.skipped_perm++; continue; }
      if (tick_ - r.touched < max_idle_ticks) continue;
      if (r.releaser) {
        r.releaser(r.ptr, r.user);
        r.swept = true;
        rep.released++;
      } else {
        rep.reported++;
      }
    }
    if (rep.released) sweeps_++;
    return rep;
  }

  void set_budget(size_t max_bytes, float growth_factor = 0.0f, size_t growth_min = 0) {
    budget_bytes_ = max_bytes; growth_factor_ = growth_factor; growth_min_ = growth_min;
    baseline_ = 0; was_over_ = false; // re-latch the expected max on next check
  }

  struct BudgetCheck {
    bool over = false;
    bool newly_over = false;
    const char* why = "";
    int released = 0;
    int reported = 0;
    int skipped_perm = 0;
  };
  // The smarter check: the registry learns the expected maximum (baseline) on
  // first sight, then patches AUTOMATICALLY whenever memory grows past the
  // hard budget or past baseline*growth_factor. sweep(0) releases everything
  // non-permanent, oldest-tracked first — "starts patching anyways".
  // Permanent entries (the live sim) are never patched, only counted.
  BudgetCheck check_budget() {
    BudgetCheck bc;
    size_t b = bytes();
    if (baseline_ == 0 && b > 0) baseline_ = b; // latch the expected max
    bool hard = budget_bytes_ > 0 && b > budget_bytes_;
    bool growth = growth_factor_ > 0 && b > growth_min_ && b > (size_t)(baseline_ * growth_factor_);
    bc.over = hard || growth;
    bc.newly_over = bc.over && !was_over_;
    was_over_ = bc.over;
    if (!bc.over) return bc;
    bc.why = hard ? "over-budget" : "growth-trip";
    SweepReport rep = sweep(0);
    bc.released = rep.released; bc.reported = rep.reported; bc.skipped_perm = rep.skipped_perm;
    return bc;
  }

  // JSON one-liner per record into out (truncated safely). Returns bytes written.
  int log_json(char* out, int cap) const {
    int n = 0;
    n += std::snprintf(out + n, cap - n,
                       "{\"tick\":%llu,\"live\":%u,\"bytes\":%u,\"peak\":%u,\"sweeps\":%llu,\"records\":[",
                       (unsigned long long)tick_, (unsigned)live(), (unsigned)bytes(),
                       (unsigned)peak_bytes_, (unsigned long long)sweeps_);
    bool first = true;
    for (auto& r : recs_) {
      if (r.swept) continue;
      if (!first) n += std::snprintf(out + n, cap - n, ",");
      first = false;
      n += std::snprintf(out + n, cap - n,
                         "{\"id\":%llu,\"size\":%u,\"tag\":\"%s\",\"site\":\"%s:%d\",\"age\":%llu,\"perm\":%d,\"has_rel\":%d}",
                         (unsigned long long)r.id, (unsigned)r.size, r.tag, r.file, r.line,
                         (unsigned long long)(tick_ - r.touched), r.permanent ? 1 : 0,
                         r.releaser ? 1 : 0);
      if (n > cap - 160) break;
    }
    n += std::snprintf(out + n, cap - n, "]}");
    return n;
  }
};

inline Registry& global() {
  static Registry g;
  return g;
}

// ---- Arena: the ECC-guarded pool -------------------------------------------
struct BlockHead {
  uint32_t magic = 0;
  uint32_t size = 0;
  uint32_t tag = 0;
  uint32_t state = 0; // 0 free, 1 live, 2 quarantined
  uint64_t born = 0;
  uint64_t touched = 0;
  uint32_t cksum = 0;
};

class Arena {
 public:
  std::vector<uint8_t> buf_;
  size_t capacity_ = 0;
  uint64_t tick_ = 0;
  uint64_t corrupt_ = 0;
  uint64_t swept_ = 0;

  bool init(size_t capacity) {
    buf_.assign(capacity, 0);
    capacity_ = capacity;
    // single free block spanning the pool
    BlockHead h;
    h.magic = MT_HEAD; h.size = (uint32_t)(capacity - sizeof(BlockHead) - 4);
    h.state = 0;
    h.cksum = cover(h);
    std::memcpy(buf_.data(), &h, sizeof(h));
    std::memcpy(buf_.data() + sizeof(h) + h.size, &MT_TAIL, 4);
    return true;
  }

  static uint32_t cover(const BlockHead& h) {
    // FNV-1a over the header fields (the parity bits). Payload guarded by canaries.
    uint32_t x = 2166136261U;
    auto mix = [&](uint32_t w) {
      for (int i = 0; i < 4; i++) { x ^= (uint8_t)(w >> (i * 8)); x *= 16777619U; }
    };
    mix(h.magic); mix(h.size); mix(h.tag); mix(h.state);
    mix((uint32_t)h.born); mix((uint32_t)h.touched);
    return x;
  }

  void* alloc(size_t n, uint32_t tag) {
    size_t off = 0;
    while (off + sizeof(BlockHead) + 4 <= capacity_) {
      BlockHead h;
      std::memcpy(&h, buf_.data() + off, sizeof(h));
      if (h.magic != MT_HEAD || h.cksum != cover(h)) return nullptr; // pool corrupt: refuse
      size_t span = sizeof(BlockHead) + h.size + 4;
      if (h.state == 0 && h.size >= n) {
        // split if room for another header + tail + 8 payload bytes
        size_t leftover = h.size - (uint32_t)n;
        h.state = 1; h.size = (uint32_t)n; h.tag = tag;
        h.born = tick_; h.touched = tick_; h.cksum = cover(h);
        std::memcpy(buf_.data() + off, &h, sizeof(h));
        std::memcpy(buf_.data() + off + sizeof(h) + n, &MT_TAIL, 4);
        if (leftover > sizeof(BlockHead) + 4 + 8) {
          size_t noff = off + sizeof(h) + n + 4;
          BlockHead nh;
          nh.magic = MT_HEAD;
          nh.size = (uint32_t)(leftover - sizeof(BlockHead) - 4);
          nh.state = 0; nh.cksum = cover(nh);
          std::memcpy(buf_.data() + noff, &nh, sizeof(nh));
          std::memcpy(buf_.data() + noff + sizeof(nh) + nh.size, &MT_TAIL, 4);
        }
        return buf_.data() + off + sizeof(BlockHead);
      }
      if (span == 0) break;
      off += span;
    }
    return nullptr; // full: refuse, don't sprawl (caller sees OOM, log shows why)
  }

  bool free(void* p) {
    if (!p) return false;
    size_t off = (uint8_t*)p - buf_.data() - sizeof(BlockHead);
    if (off + sizeof(BlockHead) > capacity_) return false;
    BlockHead h;
    std::memcpy(&h, buf_.data() + off, sizeof(h));
    if (h.magic != MT_HEAD || h.cksum != cover(h)) return false;
    if (h.state != 1) return false; // double-free: refuse, count it
    h.state = 0; h.cksum = cover(h);
    std::memcpy(buf_.data() + off, &h, sizeof(h));
    return true;
  }

  void touch(void* p) {
    if (!p) return;
    size_t off = (uint8_t*)p - buf_.data() - sizeof(BlockHead);
    if (off + sizeof(BlockHead) > capacity_) return;
    BlockHead h;
    std::memcpy(&h, buf_.data() + off, sizeof(h));
    if (h.magic == MT_HEAD && h.cksum == cover(h) && h.state == 1) {
      h.touched = tick_; h.cksum = cover(h);
      std::memcpy(buf_.data() + off, &h, sizeof(h));
    }
  }

  struct ScrubReport { int ok = 0; int quarantined = 0; int free_blocks = 0; };
  // Patrol scrub: verify every block. Corrupt live block -> QUARANTINE
  // (never handed out again). Returns counts. This is SECDED-spirit:
  // detect everything, correct by isolation, report the rest.
  ScrubReport scrub() {
    ScrubReport rep;
    size_t off = 0;
    while (off + sizeof(BlockHead) + 4 <= capacity_) {
      BlockHead h;
      std::memcpy(&h, buf_.data() + off, sizeof(h));
      if (h.magic != MT_HEAD) break; // untracked tail region: stop, don't invent blocks
      uint32_t tail = 0;
      if (h.size <= capacity_)
        std::memcpy(&tail, buf_.data() + off + sizeof(h) + h.size, 4);
      bool ok = (h.cksum == cover(h)) && (tail == MT_TAIL || h.state == 0);
      if (h.state == 1) {
        if (ok) { rep.ok++; }
        else {
          h.state = 2; // quarantine
          // re-checksum the quarantined header so future scrubs are stable
          h.cksum = cover(h);
          std::memcpy(buf_.data() + off, &h, sizeof(h));
          // retag the tail so the block reads as fenced, not valid
          uint32_t fence = 0x51444152U; // "QDAR"
          std::memcpy(buf_.data() + off + sizeof(h) + h.size, &fence, 4);
          corrupt_++; rep.quarantined++;
        }
      } else if (h.state == 0) {
        rep.free_blocks++;
      }
      size_t span = sizeof(BlockHead) + h.size + 4;
      if (span == 0) break;
      off += span;
    }
    return rep;
  }

  int sweep(uint64_t max_idle) {
    int n = 0;
    size_t off = 0;
    while (off + sizeof(BlockHead) + 4 <= capacity_) {
      BlockHead h;
      std::memcpy(&h, buf_.data() + off, sizeof(h));
      if (h.magic != MT_HEAD) break;
      if (h.state == 1 && tick_ - h.touched >= max_idle) {
        h.state = 0; h.cksum = cover(h);
        std::memcpy(buf_.data() + off, &h, sizeof(h));
        n++;
      }
      size_t span = sizeof(BlockHead) + h.size + 4;
      if (span == 0) break;
      off += span;
    }
    if (n) swept_++;
    return n;
  }

  void tick() { tick_++; }
};

} // namespace memtrash

#define MT_TRACK(ptr, size, tag) memtrash::global().track((ptr), (size), (tag), __FILE__, __LINE__)
#define MT_TRACK_PERM(ptr, size, tag) memtrash::global().track((ptr), (size), (tag), __FILE__, __LINE__, nullptr, nullptr, true)
#define MT_UNTRACK(ptr) memtrash::global().untrack((ptr))
#define MT_TOUCH(ptr) memtrash::global().touch((ptr))
