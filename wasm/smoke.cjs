// node wasm/smoke.cjs — loads build/sim.js (SINGLE_FILE, MODULARIZE), runs short
// trajectories for every preset + 3D, diffs against bench/*.json. Keep it quick.
const fs = require("fs");
const path = require("path");
const ROOT = path.join(__dirname, "..");

function loadRef(p) { return JSON.parse(fs.readFileSync(path.join(ROOT, p), "utf8")); }
function close(a, b, tol) {
  const err = Math.abs(a - b) / Math.max(1e-12, Math.abs(b));
  return { ok: err <= tol, err };
}

(async () => {
  const createModule = require("../build/sim.js");
  const M = await createModule();
  const F = (n, ret, args) => M.cwrap(n, ret, args);
  const api = {
    b2: F("topo2d_build", "number", ["number"]),
    nel2: F("topo2d_nel", "number", []),
    x2: F("topo2d_x", "number", []),
    ones2: F("topo2d_ones", null, []),
    reset2: F("topo2d_reset", null, []),
    verify2: F("topo2d_verify", null, ["number"]),
    step2: F("topo2d_step", "number", []),
    b3: F("topo3d_build", "number", []),
    nel3: F("topo3d_nel", "number", []),
    x3: F("topo3d_x", "number", []),
    ones3: F("topo3d_ones", null, []),
    reset3: F("topo3d_reset", null, []),
    verify3: F("topo3d_verify", null, ["number"]),
    step3: F("topo3d_step", "number", []),
  };
  const fails = [];
  const out = M._malloc(16);
  const cases = [
    { p: 0, ref: "bench/reference.json" },
    { p: 1, ref: "bench/reference-bridge.json" },
    { p: 2, ref: "bench/reference-mbb.json" },
    { p: 3, ref: "bench/reference-lbracket.json" },
  ];
  for (const c of cases) {
    const ref = loadRef(c.ref);
    const nel = api.b2(c.p);
    if (nel !== ref.nelx * ref.nely) { fails.push(ref.preset + ": nel " + nel); continue; }
    api.ones2();
    api.verify2(out);
    const tip = M.HEAPF64[out / 8];
    const r = close(tip, ref.tip_full_density, 0.01);
    console.log(ref.preset + " tip ref=" + ref.tip_full_density + " wasm=" + tip.toFixed(4) +
      " err=" + (r.err * 100).toFixed(2) + "% " + (r.ok ? "OK" : "FAIL"));
    if (!r.ok) fails.push(ref.preset + ": tip");
    api.reset2();
    // bench/reference.json stores SAMPLED iters [1,5,10,...]; corpus refs store every step.
    const want = (c.p === 0) ? [1, 5, 10] : [1, 2, 3];
    const maxStep = want[want.length - 1];
    for (let i = 1; i <= maxStep; i++) {
      const got = api.step2();
      const k = want.indexOf(i);
      if (k < 0) continue;
      const rr = close(got, ref.compliance_trajectory[k], 0.02);
      console.log("  iter " + i + " ref=" + ref.compliance_trajectory[k] + " wasm=" + got.toFixed(4) +
        " err=" + (rr.err * 100).toFixed(2) + "% " + (rr.ok ? "OK" : "FAIL"));
      if (!rr.ok) fails.push(ref.preset + ": iter " + i);
    }
  }
  // 3D cantilever
  {
    const ref = loadRef("bench/reference3d.json");
    const nel = api.b3();
    if (nel !== ref.nelx * ref.nely * ref.nelz) fails.push("3d: nel " + nel);
    else {
      api.ones3();
      api.verify3(out);
      const tip = M.HEAPF64[out / 8];
      const r = close(tip, ref.tip_full_density, 0.01);
      console.log("cantilever3d tip ref=" + ref.tip_full_density + " wasm=" + tip.toFixed(4) +
        " err=" + (r.err * 100).toFixed(2) + "% " + (r.ok ? "OK" : "FAIL"));
      if (!r.ok) fails.push("3d: tip");
      api.reset3();
      for (let i = 0; i < 2; i++) {
        const got = api.step3();
        const rr = close(got, ref.compliance_trajectory[i], 0.02);
        console.log("  3d iter " + (i + 1) + " ref=" + ref.compliance_trajectory[i] + " wasm=" + got.toFixed(4) +
          " err=" + (rr.err * 100).toFixed(2) + "% " + (rr.ok ? "OK" : "FAIL"));
        if (!rr.ok) fails.push("3d: iter " + (i + 1));
      }
    }
  }
  M._free(out);
  if (fails.length) { console.log("SMOKE FAIL: " + fails.join("; ")); process.exit(1); }
  console.log("SMOKE PASS: wasm core matches all references.");
})().catch((e) => { console.error("SMOKE ERROR: " + (e && e.stack || e)); process.exit(1); });
