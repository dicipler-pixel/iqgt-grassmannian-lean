// SCRIPT: GRASSMANN-MIXER-CORE-TEST
// Kill conditions written before the run (each must match the paper or the Python checks):
// C1 Jacobi eigen: residual ||A V - V diag(w)|| < 1e-10 on random 6x6 and 12x12
// C2 polar family geodesic from y=0.6 at 35 deg turns at y = 0.32981248 (paper S8), |err| < 1e-5
// C3 gap-narrowing geodesic from y=0.8 at 50 deg turns at y = -1.0969 (paper 8.6), |err| < 1e-3
// C4 five swaps vs 10-cycle on 20 points: sum sin^2 = 5 and 5; lengths 3.512407 and 2.896405 (Python S3), |err| < 1e-5
// C5 speed never exceeds the ceiling on 200 random feed points
// C6 envelope: alpha = pi/2, mix = 0 sits on the bound (ratio 1 within 1e-12); random mix falls inside (ratio <= 1)
// C7 tail stays under its Davis-Kahan bound in the corridor setup: subspace angle median < tip median at inner gap 0.01
const C = require('./core.js');
const out = {};
{ const r = C.rng(5); let worst = 0;
  for (const n of [6, 12]) for (let t = 0; t < 20; t++) { const A = C.randSym(n, r), e = C.eigh(A);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { let s = 0; for (let k = 0; k < n; k++) s += A[i][k] * e.V[k][j]; worst = Math.max(worst, Math.abs(s - e.V[i][j] * e.w[j])); } }
  out.C1 = [worst, worst < 1e-10]; }
{ const f = C.families.polar, p = C.geodesic(f.a, f.b, 0, 0.6, 35, 2.0, 20000); const ymin = Math.min(...p.map(q => q[1])); out.C2 = [ymin, Math.abs(ymin - 0.32981248) < 1e-5]; }
{ const f = C.families.narrow, p = C.geodesic(f.a, f.b, 0, 0.8, 50, 3.0, 30000); const ymin = Math.min(...p.map(q => q[1])); out.C3 = [ymin, Math.abs(ymin + 1.0969) < 1e-3]; }
{ const a = C.permStats([2, 2, 2, 2, 2], 20), b = C.permStats([10], 20);
  out.C4 = [[a.s2, b.s2, a.L, b.L], Math.abs(a.s2 - 5) < 1e-9 && Math.abs(b.s2 - 5) < 1e-9 && Math.abs(a.L - 3.512407) < 1e-5 && Math.abs(b.L - 2.896405) < 1e-5]; }
{ const F = C.makeFeed(11), r = C.rng(3); let worst = 0;
  for (let t = 0; t < 200; t++) { const sp = C.speed(F, 0.05 + r() * 0.5, 0.5 + r() * 2, -1 + 2 * r(), 1 + Math.floor(r() * 3)); worst = Math.max(worst, sp.v / sp.ceil); }
  out.C5 = [worst, worst <= 1]; }
{ const r = C.rng(9); const on = C.envelopePoint(r, 2, 4, Math.PI / 2, 0); let mx = 0; for (let t = 0; t < 500; t++) mx = Math.max(mx, C.envelopePoint(r, 2, 4, r() * 6.28, r()).ratio);
  out.C6 = [[on.ratio, mx], Math.abs(on.ratio - 1) < 1e-12 && mx <= 1 + 1e-12]; }
{ const F = C.makeFeed(4), r = C.rng(7); const tips = [], tails = [];
  for (let t = 0; t < 60; t++) { const rows = C.tipTail(C.feedOp(F, 0.01, 3, 0), 0.0015, r); tips.push(rows[1].vec); tails.push(rows[2].sub); }
  const med = a => a.slice().sort((x, y) => x - y)[a.length >> 1];
  out.C7 = [[med(tips), med(tails)], med(tails) < med(tips)]; }
for (const [k, [v, ok]] of Object.entries(out)) console.log(k, ok ? 'PASS' : 'FAIL', JSON.stringify(v));
