#!/usr/bin/env python3
"""Differential test: candidate.json (C++ bench) vs reference.json (JS bench).
Fails if structure drifted (compliance trajectory or tip out of tolerance).
Usage: python3 bench/check.py bench/reference.json candidate.json
"""
import json, sys

ref = json.load(open(sys.argv[1]))
got = json.load(open(sys.argv[2]))

assert got.get("contract") == ref["contract"] == "topo-bench-v1", "contract mismatch"
fails = []

def cmp(name, a, b, tol):
    err = abs(a - b) / max(1e-12, abs(b))
    ok = err <= tol
    print(f"{name:28s} ref={b:<14.4f} got={a:<14.4f} err={err*100:6.2f}% {'OK' if ok else 'FAIL'}")
    if not ok:
        fails.append(name)

print("== tip (full density) ==")
cmp("tip_full_density", got["tip_full_density"], ref["tip_full_density"], 0.01)

print("== compliance trajectory ==")
for i, (a, b) in enumerate(zip(got["compliance_trajectory"], ref["compliance_trajectory"])):
    cmp(f"iter {[1,5,10,15,20,25,30,35,40][i]}", a, b, 0.02)

print("== timing (informational only) ==")
print(f"JS reference : {ref.get('js_time_ms')} ms")
print(f"C++ candidate: {got.get('cpp_time_ms')} ms")
jt, ct = ref.get("js_time_ms") or 0, got.get("cpp_time_ms") or 0
if jt and ct:
    print(f"speedup: {jt/ct:.1f}x")

if fails:
    print(f"\nFAIL: {len(fails)} checks out of tolerance: {fails}")
    sys.exit(1)
print("\nPASS: structure matches reference within tolerance.")
