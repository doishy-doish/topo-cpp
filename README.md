# topo-cpp

C++ strangler ports of the [topo-toy](../topo-toy) SIMP topology optimizer — C++17, zero dependencies.

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

The JS repo is the reference implementation; this repo is the performance lineage. Timing from both is printed side by side on every run — the speedup number is the scoreboard.

Scoreboard (different iron, honest note): 2D cantilever C++ 1413ms vs JS 4020ms; 3D C++ 719ms/4 threads vs JS 3463ms.
