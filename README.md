# topo-cpp

C++ strangler port of the [topo-toy](../topo-toy) SIMP topology optimizer — C++17, zero dependencies.

## the rule between the two repos

Both repos implement the same bench contract (`bench/reference.json`, `topo-bench-v1`): cantilever 120x20, 50% volume, 40 OC steps. CI builds with `g++ -O3`, runs the bench, and differential-tests the compliance trajectory + tip deflection against the JS reference (2% / 1% tolerances). **Any minimization that breaks the structure fails CI instead of silently drifting.**

```
g++ -O3 -std=c++17 -o topo main.cpp
./topo            # prints candidate.json
python3 bench/check.py bench/reference.json candidate.json
```

The JS repo is the reference implementation; this repo is the performance lineage. Timing from both is printed side by side on every run — the speedup number is the scoreboard.
