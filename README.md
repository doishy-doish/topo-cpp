# topo-cpp

C++ strangler ports of the [topo-toy](../topo-toy) SIMP topology optimizer â C++17, zero dependencies.

## the rule between the two repos

Both repos implement the same bench contracts: `topo-bench-v1` (cantilever 120x20, 50% volume, 40 OC steps), `topo-bench-v2` (3D cantilever 40x12x6), plus the corpus (`topo-bench-bridge-v1`, `topo-bench-mbb-v1`, `topo-bench-lbracket-v1`). CI builds with `g++ -O3`, runs each bench, and differential-tests the compliance trajectory + tip deflection against the JS reference (tip 1%, trajectory 2%). **Any minimization that breaks the structure fails CI instead of silently drifting.**

```
g++ -O3 -std=c++17 -o topo main.cpp
./topo            # prints candidate.json
python3 bench/check.py bench/reference.json candidate.json

g++ -O3 -std=c++17 -o mbb mbb.cpp && ./mbb | tee candidate-mbb.json
python3 bench/check_corpus.py bench/reference-mbb.json candidate-mbb.json
# bridge, lbracket: same pattern
```

The JS repo is the reference implementation; this repo is the performance lineage. Timing from both is printed side by side on every run â the speedup number is the scoreboard.

Scoreboard (different iron, honest note): 2D cantilever C++ 1413ms vs JS 4020ms; 3D C++ 719ms/4 threads vs JS 3463ms.

## unified wasm core

`wasm/sim.cpp` is the single SIMP core (2D corpus + 3D cantilever) behind a C
API: solve/step accessors plus opt-in pro controls
(`topo2d_set_projection/get_projection`, `topo2d_set_filter/get_filter`,
`topo2d_set_penal_cont/get_penal`, `topo2d_compute_stress/topo2d_stress`)
and memory telemetry (`topo_mt_live/bytes/peak/tick/sweep/log`,
`set_budget/check_budget`). CI compiles it twice: natively (all 5 contracts
green, trajectories identical to 4 decimals) and with Emscripten to
SINGLE_FILE `build/sim.js` (~80-104KB), smoke-tested in node and shipped to
the toy repo as `vendor/sim-wasm.js`. `test_pro.cpp` proves the pro paths
(penal schedule, finite nonneg stress) in CI.

## memguard

`memtrash.hpp` (Registry: track/untrack/sweep/log_json + Arena: canary +
checksum scrub, quarantine, idle sweep) with `MEMTRASH.md` design doc and a
dedicated CI job (hard-cap + growth-trip + permanent-survives tests). The
in-page mirror is the toy repo's `vendor/memtrash.js`.
