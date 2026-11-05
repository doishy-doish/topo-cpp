#!/usr/bin/env python3
"""Differential test for topo-bench-v2 (3D). Same shape as bench/check.py."""
import json, sys

ref = json.load(open(sys.argv[1]))
got = json.load(open(sys.argv[2]))

assert got.get("contract") == ref["contract"] == "topo-bench-v2", "contract mismatch"
fails = []

def cmp(name, a, b, tol):
    err = abs(a - b) / max(1e-12, abs(b))
    ok = err <= tol
    print(f"{name:28s} ref={b:<14.4f} got={a:<14.4f} err={err*100:6.2f}% {'OK' if ok else 'FAIL'}")
    if not ok:
        fails.append(name)

print("== tip (full density) ==")
cmp("tip_full_density", got["tip_full_density"], ref["tip_full_density"], 0.005)
print("== compliance trajectory ==")
for i, (a, b) in enumerate(zip(got["compliance_trajectory"], ref["compliance_trajectory"])):
    cmp(f"iter {i + 1}", a, b, 0.02)
print("== timing/threads (informational) ==")
print(f"JS reference : {ref.get('js_time_ms')} ms")
print(f"C++ candidate: {got.get('cpp_time_ms')} ms on {got.get('threads')} threads")
if fails:
    print(f"\nFAIL: {fails}")
    sys.exit(1)
print("\nPASS: 3D structure matches reference within tolerance.")
