// memtrash_test.cpp — native proof for CI (g++ -O3 -std=c++17, no deps).
// Proves: registry track/sweep (patch via releaser, report-not-free without),
// permanent skip, arena alloc/scrub detects planted corruption + quarantines,
// arena sweep frees idle blocks, double-free refused, log_json names areas.
// Exit 0 = PASS, nonzero = FAIL (CI reads it).
#include <cstdio>
#include <cstring>
#include "memtrash.hpp"

static int failures = 0;
#define CHECK(cond, name) do { \
  if (!(cond)) { failures++; std::printf("FAIL %s\n", name); } \
  else { std::printf("ok %s\n", name); } } while (0)

static int released_count = 0;
static void test_releaser(const void* p, void* u) {
  (void)p; (void)u; released_count++;
}

int main() {
  // ---- Registry ----
  {
    memtrash::Registry r;
    int a = 1, b = 2, c = 3;
    CHECK(r.track(&a, 100, "chunk-full", "world.cpp", 10, test_releaser, nullptr, false), "reg: track leakable");
    CHECK(r.track(&b, 200, "gpu-mesh", "render.cpp", 20, test_releaser, nullptr, false), "reg: track second");
    CHECK(r.track(&c, 300, "sim-core", "sim.cpp", 30, nullptr, nullptr, true), "reg: track permanent");
    CHECK(r.live() == 3, "reg: 3 live");
    CHECK(r.bytes() == 600, "reg: bytes add up");
    for (int i = 0; i < 10; i++) r.tick();
    r.touch(&a); // a is used; b leaks
    for (int i = 0; i < 10; i++) r.tick(); // b idle 20 ticks, a idle 10
    memtrash::Registry::SweepReport rep = r.sweep(15);
    CHECK(rep.released == 1 && released_count == 1, "reg: sweep patches the leak via releaser");
    CHECK(rep.skipped_perm == 1, "reg: permanent skipped");
    CHECK(r.live() == 2, "reg: live counts hold");
    // ambiguous: no releaser -> reported, NOT freed
    int d = 4;
    r.track(&d, 50, "unknown", "???", 0);
    for (int i = 0; i < 100; i++) r.tick();
    r.touch(&a); // a is still in use: the sweep must spare it
    memtrash::Registry::SweepReport rep2 = r.sweep(5);
    CHECK(rep2.reported == 1 && rep2.released == 0, "reg: ambiguous reported, not freed");
    char buf[2048];
    int n = r.log_json(buf, sizeof(buf));
    CHECK(n > 0 && std::strstr(buf, "sim-core") != nullptr, "reg: log names the areas");
  }
  // ---- Budget: the expected maximum, then automatic patching ----
  {
    memtrash::Registry r;
    int a = 1, b = 2, c = 3;
    r.track(&a, 100, "temp-a", "t.cpp", 1, test_releaser);
    r.track(&b, 100, "temp-b", "t.cpp", 2, test_releaser);
    r.track(&c, 100, "sim-core", "t.cpp", 3, nullptr, nullptr, true);
    int before = released_count;
    r.set_budget(250); // expected max 250, but 300 live
    memtrash::Registry::BudgetCheck bc = r.check_budget();
    CHECK(bc.over && bc.newly_over, "budget: over-budget trip fires");
    CHECK(bc.released == 2 && released_count - before == 2, "budget: auto-patch releases both temp blocks");
    CHECK(bc.skipped_perm == 1 && r.live() == 1, "budget: permanent sim survives the patch");
    memtrash::Registry::BudgetCheck bc2 = r.check_budget();
    CHECK(!bc2.over && !bc2.newly_over, "budget: quiet once back under");
    // growth tripwire: latch baseline, then grow 2x
    int d = 4, e = 5;
    r.track(&d, 100, "temp-d", "t.cpp", 4, test_releaser);
    r.set_budget(0, 1.5f, 50); // no hard cap; trip at 1.5x baseline over 50B
    memtrash::Registry::BudgetCheck g1 = r.check_budget(); // latches 200 (sim 100 + d 100)
    CHECK(!g1.over, "budget: baseline latch is quiet");
    r.track(&e, 150, "temp-e", "t.cpp", 5, test_releaser); // 350 > 300 = 1.5x of 200
    memtrash::Registry::BudgetCheck g2 = r.check_budget();
    CHECK(g2.over && g2.why && std::strcmp(g2.why, "growth-trip") == 0, "budget: growth trip fires");
    CHECK(g2.released == 2, "budget: growth auto-patches the temps");
  }
  // ---- Arena ----
  {
    memtrash::Arena ar;
    CHECK(ar.init(4096), "arena: init");
    void* p1 = ar.alloc(64, 7);
    void* p2 = ar.alloc(128, 8);
    CHECK(p1 && p2, "arena: allocs succeed");
    ar.tick(); ar.touch(p1);
    memtrash::Arena::ScrubReport s0 = ar.scrub();
    CHECK(s0.ok == 2 && s0.quarantined == 0, "arena: clean scrub");
    // plant corruption: flip a payload byte inside p2 (canary span still valid,
    // so corrupt the header checksum path instead: poke p2's header tag)
    uint8_t* raw = (uint8_t*)p2 - sizeof(memtrash::BlockHead);
    raw[8] ^= 0xFF; // damage the tag field -> checksum must fail
    memtrash::Arena::ScrubReport s1 = ar.scrub();
    CHECK(s1.quarantined == 1, "arena: scrub detects + quarantines corruption");
    memtrash::Arena::ScrubReport s2 = ar.scrub();
    CHECK(s2.quarantined == 0 && s2.ok == 1, "arena: quarantine is stable, healthy block untouched");
    // idle sweep: p1 touched at tick 1, advance far, sweep frees it; quarantined block stays
    for (int i = 0; i < 50; i++) ar.tick();
    int freed = ar.sweep(10);
    CHECK(freed == 1, "arena: sweep patches idle leak");
    CHECK(!ar.free(p1), "arena: double-free refused");
    void* p3 = ar.alloc(64, 9);
    CHECK(p3 != nullptr, "arena: pool still serves after quarantine");
  }
  std::printf(failures ? "RESULT FAIL %d\n" : "RESULT PASS\n", failures);
  return failures ? 1 : 0;
}
