# MEMORY TRASH SYSTEM — ECC for memory (the big idea, written down)

Date: 2026-11-05. Status: designed + built (C++ `memtrash.hpp`, JS `memtrash.js`,
tests, wired into topo-cpp/topo-toy). Offered to andra-voxel as a drop-in.

## 1. The problem, plainly

Trusting code is not a memory strategy. Every non-trivial program leaks a little:
a forgotten `delete`, a rebuilt solver whose old buffers nobody untracked, a GPU
buffer orphaned when its chunk unloaded, a cache with no eviction. The failures
are silent — no crash, just slow growth until the tab dies or the run OOMs.
Reviews catch some. They never catch all, because the leak is in the path nobody
re-read.

## 2. The idea: ECC, but for liveness

ECC RAM does not trust bits. Every word carries extra parity bits; a patrol
scrub reads each word, checks parity, corrects a single-bit error on the spot
(SECDED), and flags the word uncorrectable when parity can't cover it. Nobody
asks the program for permission. The correction happens *below* the program.

The memory trash system does the same thing one level up, for *liveness*:

- every allocation carries a header — its parity bits (magic canaries +
  checksum over the header, tag, birth tick, owner site);
- a patrol scrub walks every live block, verifies the header, and
  - **corrects** what it can: idle blocks past their age get released through
    their registered releaser (the leak patch); corrupt blocks get quarantined
    so they can never be handed out again (the ECC "spare the row" move);
  - **flags** what it can't fix: ambiguous ownership is logged, never freed
    blindly. A blind free is corruption, not correction.
- every allocation is logged in a registry: pointer, size, tag, site
  (file:line), age. The whole memory map is inspectable at any time
  (`audit()` / `mt_log()` / `window.__memtrash.log()`).

So: **log all the memory areas, scrub them on a patrol, patch what is safe to
patch, quarantine the rest, and leave a paper trail.** That is the ECC analogy,
and it holds up: detect (canary+checksum), correct (releaser / quarantine),
scrub (patrol), report (log). What it does NOT do is also ECC-honest: parity
can't fix a double-bit error, and this can't fix a logic bug. A use-after-free
in program logic, a wrong-size computation, a design that genuinely needs 2 GB
— those stay human problems. The system converts *silent* leaks into *logged,
attributed, usually auto-fixed* events. That is the whole claim. Nothing more.

## 3. What andra-voxel already does (the scraps method — learned, not copied)

Pulled the 4 MB single-file game and studied the discipline, not the code:

- **Rule 5: no allocation on the hot path.** Module-level scratch arrays
  (21 KiB visited map, 85 KiB stack, per-pass scratch) allocated once at load,
  reused per chunk. Passes that allocate per chunk were measured, then killed.
- **Streaming pools.** Chunk block/height/floor-cache/flora arrays are recycled
  on unload (`_fullPool`, `_heightPool`, `_floorPool`, `_floraPool`) instead of
  re-allocated, "so travelling the world does not churn the JS heap."
  Wasm linear memory never grows; typed-array views are created once and reused.
- **Hard caps everywhere.** `VILLAGER_CAP = 12` alive, board caps
  (512 KiB blob, 400 votes, 300 notes), chest/sign/edit caps, skin map capped
  at 128 with oldest-evict, voice PCM cache 8 MB with FIFO evict, lazy
  DECtalk load + `unload()` that drops the module AND the cache.
- **Ownership callbacks.** `unloadFar(..., onUnload)` keeps the renderer in
  sync; `this.dirty.delete(c)` on unload "else it leaks in `dirty`" — the
  codebase names its own leak sites in comments.
- **Debug surface.** `memAudit()` (chunks, meshes, orphans, sheds, pool len),
  `memShed(reason)` with 60 s warn throttle, heap-growth tripwire
  (8 ticks, >256 MB, >1.35x baseline → shed). Probes over arguments:
  `__dbg.*` is how claims get checked.
- **Fixed arenas with documented geometry** (`ARENA_BYTES`, water sim slots,
  post-arena regions marked NOT FREE) instead of malloc sprawl.

That is a box-of-scraps memory system: fixed buffers, pools, caps, lazy
load/unload, and a debug readout — all verifiable in-page, no external
toolchain. The memory trash system below is that discipline made *portable and
automatic*: same rules, but the registry/scrub/sweep enforces them instead of
relying on every author remembering.

## 4. The two layers

**Layer A — Registry (log all memory areas).** `track(ptr, size, tag, site,
releaser?, permanent?)` / `untrack(ptr)` / `touch(ptr)` / `audit()` /
`sweep(maxIdleTicks)` / `log_json()`. Sweep calls the releaser for stale
non-permanent entries, then untracks. Permanent entries (the live sim) are
never swept — they are audited, not patched. Ambiguous entries (no releaser)
are reported, never freed.

**Layer B — Arena (ECC-guarded pool for new allocations).**
Fixed-capacity bump pool. Block layout: `[header | payload | tail magic]`.
Header = {magic, size, tag, state, born, checksum}. Scrub verifies magic +
checksum per live block; corrupt → state=QUARANTINED (never reused, counted,
logged). `sweep()` releases idle blocks. `touch()` on use. This is where
"actually patches" lives for new code: the pool owns its memory, so the patch
is safe by construction.

JS mirror (`memtrash.js`): same shape for GPU/pool/chunk resources —
`track(key, {size, tag, dispose, permanent})`, orphan predicate, age sweep,
`log()`. Explicit `dispose` (like andra-voxel's `onUnload`), never
FinalizationRegistry guessing.

## 5. Wiring (what "properly" means)

- topo-cpp: `memtrash.hpp` at root, `memtrash_test.cpp` proving
  track/sweep/corrupt-quarantine, `memguard` CI job (native compile + run).
  `wasm/sim.cpp` tracks the live sim buffers (permanent) and exports
  `topo_mt_live/bytes/peak/sweep/log` (emcc `EXPORTED_FUNCTIONS` updated in
  the same push — the lesson from the VOLFRAC-namespace miss: header and
  wiring ship atomically).
- topo-toy: `vendor/memtrash.js` + `window.__memtrash`; solve buffers tracked,
  MEM readout next to the engine tag, sweep on rebuild.
- andra-voxel: NOT touched (agent3's game, infra-on-request only). Offered as
  a drop-in: replace/augment `memAudit` with the registry, keep every cap.

## 6. On needing wandbox (the honest answer)

The question is fair. Two different things got mixed together:

1. **Genuinely needed a second iron.** The topo claim was *bit-identity*
   across JS, C++, and wasm ("same trajectory to 4 decimals on three machines").
   You cannot prove that on one machine. Wandbox/Godbolt/CI runners were the
   second and third machines: different compilers, different CPUs checking the
   same bytes. andra-voxel never needed that because it never made a
   cross-iron claim — it runs in one place (the browser) and its probes run
   there too.
2. **Wastefully done.** Too many blind pushes, too many round-trips, reading
   logs I couldn't reach (workflow logs are admin-gated; artifacts aren't —
   learned late). The scraps fix, applied from here on: verify in-workspace
   first (JS tests run here in ms), then ONE CI run as the second iron, with
   artifacts (`verify-log.txt`, `shots/`) as the readout. That is exactly how
   this memory system shipped: JS self-test here, C++ proof on our own public
   CI, no third-party pastebins.

Box of scraps, cave rules, adopted: fixed buffers, pools, caps, lazy
load/unload, probes over arguments, and no trip outside the workspace unless
the claim itself requires a second machine.

## 7. The smarter check (budget + growth tripwire, added same night)

Knowing the expected maximum makes the system proactive instead of merely
reactive. `set_budget(maxBytes, growthFactor, growthMin)` + `check_budget()`:

- the registry **latches the baseline** (expected max) on first sight;
- every check after that trips when live bytes pass the hard budget
  (`over-budget`) OR pass `baseline * growthFactor` above `growthMin`
  (`growth-trip`, the andra-voxel 1.35x-tripwire idea made explicit);
- on a trip it **patches anyways**: `sweep(0)` releases everything
  non-permanent, oldest-tracked first. The live sim is permanent — audited,
  never patched. If nothing is releasable the trip stays `over` and the
  report says so honestly (can't patch the live sim; report it);
- `newly_over` edge detection so patrols warn once per trip instead of
  spamming while over-budget with nothing left to release.

The page runs this as a 30 s patrol (`startShed` now checks age AND budget):
expected max learned live, 256 MB hard cap, 1.5x growth over 32 MB minimum.
Sheds log to console as `[memtrash] shed {...}`. C++ and JS share the exact
same semantics; both are covered by tests (hard-cap auto-patch, baseline
latch quiet, growth-trip auto-patch, permanent survives).
