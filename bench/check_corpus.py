#!/usr/bin/env python3
"""Generic differential test for corpus presets (bridge/mbb/lbracket).
Compares tip (1%) and every compliance-trajectory entry (2%).
Usage: python3 bench/check_corpus.py bench/reference-<preset>.json candidate-<preset>.json
"""
import json, sys
ref = json.load(open(sys.argv[1]))
got = json.load(open(sys.argv[2]))
assert got.get("contract") == ref["contract"], "contract mismatch: %r vs %r" % (got.get("contract"), ref["contract"])
assert got.get("preset") == ref["preset"], "preset mismatch"
assert got.get("nelx") == ref["nelx"] and got.get("nely") == ref["nely"], "dims mismatch"
fails = []
def cmp(name, a, b, tol):
    err = abs(a - b) / max(1e-12, abs(b))
    ok = err <= tol
    print("%-28s ref=%-14.4f got=%-14.4f err=%6.2f%% %s" % (name, b, a, err*100, "OK" if ok else "FAIL"))
    if not ok: fails.append(name)
print("== tip (full density) ==")
cmp("tip_full_density", got["tip_full_density"], ref["tip_full_density"], 0.01)
print("== compliance trajectory (%d steps) ==" % len(ref["compliance_trajectory"]))
assert len(got["compliance_trajectory"]) == len(ref["compliance_trajectory"]), "trajectory length mismatch"
for i, (a, b) in enumerate(zip(got["compliance_trajectory"], ref["compliance_trajectory"])):
    cmp("iter %d" % (i+1), a, b, 0.02)
print("== timing (informational) ==")
print("JS reference : %s ms" % ref.get("js_time_ms"))
print("C++ candidate: %s ms" % got.get("cpp_time_ms"))
if fails:
    print("\nFAIL: %s" % fails)
    sys.exit(1)
print("\nPASS: %s matches reference within tolerance." % ref["preset"])
