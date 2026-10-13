// SCRIPT: GRASSMANN-MIXER-CORE-V2
// Math core for the Grassmannian mixer. Pure functions, no DOM. Tested by core_test.js against the paper's numbers.
const Core = (() => {
  // seeded RNG (mulberry32) so every run of the page is reproducible
  function rng(seed) { let a = seed >>> 0; return () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
  function gauss(r) { let u = 0, v = 0; while (u === 0) u = r(); v = r(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
  const zeros = (n, m) => Array.from({ length: n }, () => new Float64Array(m));
  const eye = n => { const A = zeros(n, n); for (let i = 0; i < n; i++) A[i][i] = 1; return A; };
  function matmul(A, B) { const n = A.length, m = B[0].length, p = B.length, C = zeros(n, m); for (let i = 0; i < n; i++) for (let k = 0; k < p; k++) { const a = A[i][k]; if (a) for (let j = 0; j < m; j++) C[i][j] += a * B[k][j]; } return C; }
  const T = A => { const n = A.length, m = A[0].length, B = zeros(m, n); for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) B[j][i] = A[i][j]; return B; };
  // cyclic Jacobi for real symmetric matrices; returns ascending eigenvalues and column eigenvectors
  function eigh(A0) {
    const n = A0.length, A = A0.map(r => Float64Array.from(r)), V = eye(n);
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += A[i][j] * A[i][j];
      if (off < 1e-30) break;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-300) continue;
        const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) { const akp = A[k][p], akq = A[k][q]; A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq; }
        for (let k = 0; k < n; k++) { const apk = A[p][k], aqk = A[q][k]; A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk; }
        for (let k = 0; k < n; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
      }
    }
    const idx = [...Array(n).keys()].sort((a, b) => A[a][a] - A[b][b]);
    return { w: idx.map(i => A[i][i]), V: V.map(r => Float64Array.from(idx.map(i => r[i]))) };
  }
  const cols = (V, ids) => V.map(r => Float64Array.from(ids.map(i => r[i])));
  // singular values of a small matrix via eigenvalues of M^T M
  function svals(M) { const e = eigh(matmul(T(M), M)).w; return e.map(x => Math.sqrt(Math.max(0, x))); }
  // principal angles between column spaces of orthonormal U (n x k) and W (n x k), radians, descending
  function principal(U, W) { const s = svals(matmul(T(U), W)); return s.map(x => Math.acos(Math.min(1, x))).sort((a, b) => b - a); }
  const deg = r => r * 180 / Math.PI;
  function randSym(n, r, scale = 1) { const A = zeros(n, n); for (let i = 0; i < n; i++) for (let j = i; j < n; j++) { const g = gauss(r) * scale; A[i][j] = g; A[j][i] = g; } return A; }
  function randOrth(n, r) { return eigh(randSym(n, r)).V; }
  function opnorm(A) { const e = eigh(matmul(T(A), A)).w; return Math.sqrt(Math.max(...e)); }
  // the feed: a spectrum with a 3-level cluster (inner gap) under an outer gap, in a random frame, plus a sweep direction
  function makeFeed(seed) { const r = rng(seed); return { Q: randOrth(6, r), H1: randSym(6, r, 0.35), r }; }
  function feedOp(F, inner, outer, s) {
    const w = [0, inner, 2 * inner, 2 * inner + outer, 2 * inner + outer + 0.6, 2 * inner + outer + 1.3];
    const D = zeros(6, 6); for (let i = 0; i < 6; i++) D[i][i] = w[i];
    const A = matmul(matmul(F.Q, D), T(F.Q));
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) A[i][j] += s * F.H1[i][j];
    return A;
  }
  // clusters: consecutive levels whose relative gap is below thr fuse into one subspace
  function clusters(w, thr) { const span = (w[w.length - 1] - w[0]) || 1; const out = [[0]]; for (let i = 1; i < w.length; i++) { if ((w[i] - w[i - 1]) / span < thr) out[out.length - 1].push(i); else out.push([i]); } return out; }
  // tip and tail between two independent noisy halves of the same operator
  function tipTail(A, noise, r) {
    const n = A.length, Ea = randSym(n, r, noise), Eb = randSym(n, r, noise);
    const add = (X, Y) => X.map((row, i) => row.map((v, j) => v + Y[i][j]));
    const a = eigh(add(A, Ea)), b = eigh(add(A, Eb)), base = eigh(A);
    const rows = [];
    for (let k = 1; k < n; k++) {
      const vec = principal(cols(a.V, [k - 1]), cols(b.V, [k - 1]))[0];
      const sub = principal(cols(a.V, [...Array(k).keys()]), cols(b.V, [...Array(k).keys()]))[0];
      const gap = (base.w[k] - base.w[k - 1]) / ((base.w[n - 1] - base.w[0]) || 1);
      rows.push({ k, vec: deg(vec), sub: deg(sub), gap });
    }
    return rows;
  }
  // metric speed of the rank-k projector along s, and the ceiling sqrt(k)||H1||/gap of Theorem 5.2
  function speed(F, inner, outer, s, k) {
    const h = 1e-5, P = x => { const e = eigh(feedOp(F, inner, outer, x)); const U = cols(e.V, [...Array(k).keys()]); return { P: matmul(U, T(U)), w: e.w }; };
    const p1 = P(s + h), p0 = P(s - h), pc = P(s); let g = 0;
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) { const d = (p1.P[i][j] - p0.P[i][j]) / (2 * h); g += 0.5 * d * d; }
    const ceil = Math.sqrt(k) * opnorm(F.H1) / (pc.w[k] - pc.w[k - 1]);
    return { v: Math.sqrt(g), ceil };
  }
  // geodesics of a diagonal 2D metric g = diag(a(x,y), b(x,y)), RK4; returns path and the Clairaut quantity a*sin^2(theta)
  function geodesic(a, b, x0, y0, thetaDeg, T1, steps, sgn = -1) {
    const d = 1e-6;
    const ax = (x, y) => (a(x + d, y) - a(x - d, y)) / (2 * d), ay = (x, y) => (a(x, y + d) - a(x, y - d)) / (2 * d);
    const bx = (x, y) => (b(x + d, y) - b(x - d, y)) / (2 * d), by = (x, y) => (b(x, y + d) - b(x, y - d)) / (2 * d);
    const f = s => { const [x, y, u, v] = s, A = a(x, y), B = b(x, y);
      const xa = -(ax(x, y) * u * u + 2 * ay(x, y) * u * v - bx(x, y) * v * v) / (2 * A);
      const ya = -(-ay(x, y) * u * u + 2 * bx(x, y) * u * v + by(x, y) * v * v) / (2 * B);
      return [u, v, xa, ya]; };
    const th = thetaDeg * Math.PI / 180;
    let s = [x0, y0, Math.sin(th) / Math.sqrt(a(x0, y0)), sgn * Math.cos(th) / Math.sqrt(b(x0, y0))];
    const dt = T1 / steps, path = [[s[0], s[1]]];
    for (let i = 0; i < steps; i++) {
      const k1 = f(s), k2 = f(s.map((v, j) => v + dt / 2 * k1[j])), k3 = f(s.map((v, j) => v + dt / 2 * k2[j])), k4 = f(s.map((v, j) => v + dt * k3[j]));
      s = s.map((v, j) => v + dt / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
      if (!s.every(Number.isFinite) || Math.abs(s[1]) > 50) break;
      path.push([s[0], s[1]]);
    }
    return path;
  }
  // three two-band families from the paper: metric components in (x, y)
  const families = {
    polar: { label: 'Subspace stops rotating (polar stratum)', a: (x, y) => 9 * Math.sin(y) ** 2 / 4, b: () => 0.25, y0: 0.6, view: [-0.1, 2.6, -0.15, 0.85], line: 0,
      gap: () => 2, note: 'g∥ = (9/4) sin²y falls to zero on y = 0: the whole line is one projector.' },
    narrow: { label: 'Gap narrows (guiding line)', a: (x, y) => 9 / (4 * (1 + 6.25 * y * y)), b: (x, y) => 6.25 / (4 * (1 + 6.25 * y * y) ** 2), y0: 0.8, view: [-0.1, 4.2, -1.6, 1.6], line: 0,
      gap: (x, y) => 2 * Math.sqrt(1 + 6.25 * y * y), note: 'Gap 2√(1+6.25y²) is narrowest on y = 0, where the index is largest.' },
    fold: { label: 'Fold (collapse across the line)', a: (x, y) => x * x * Math.sin(y + 1.1) ** 2 / 4 + 1e-9, b: () => 0.25, y0: 0.6, view: [-1.2, 1.2, -0.3, 0.9], line: 'x0',
      gap: () => 2, note: 'det g → 0 on x = 0 while the projector keeps moving along it.' }
  };
  // permutations: eigenphases of a permutation from its cycle type; principal angles between graph(I) and graph(pi)
  function cycleAngles(cycles, n) { // cycles: array of lengths; remaining points fixed
    const th = []; let used = 0;
    for (const l of cycles) { for (let j = 0; j < l; j++) { const phi = 2 * Math.PI * j / l; const p = phi > Math.PI ? 2 * Math.PI - phi : phi; th.push(p / 2); } used += l; }
    for (let i = used; i < n; i++) th.push(0);
    return th;
  }
  function permStats(cycles, n) {
    const th = cycleAngles(cycles, n), moved = cycles.reduce((a, l) => a + (l > 1 ? l : 0), 0);
    const s2 = th.reduce((a, t) => a + Math.sin(t) ** 2, 0), t2 = th.reduce((a, t) => a + t * t, 0);
    const gram = []; for (const l of cycles) for (let j = 0; j < l; j++) gram.push({ l, j, val: 2 * Math.sin(Math.PI * j / l) ** 2 });
    return { th, dH: moved / n, s2, L: Math.sqrt(t2), perPoint: t2 / n, gram };
  }
  // curvature bound at a projector: blocks X, Y (k x m complex); g = Re Tr(X Y*), Omega = -2 Im Tr(X Y*)
  function envelopePoint(r, k, m, alpha, mix) {
    const cplx = () => ({ re: Array.from({ length: k * m }, () => gauss(r)), im: Array.from({ length: k * m }, () => gauss(r)) });
    const X = cplx(), Z = cplx();
    const ca = Math.cos(alpha), sa = Math.sin(alpha);
    const Y = { re: X.re.map((v, i) => ca * v - sa * X.im[i] + mix * Z.re[i]), im: X.im.map((v, i) => ca * v + sa * X.re[i] + mix * Z.im[i]) };
    const tr = (A, B) => { let re = 0, im = 0; for (let i = 0; i < A.re.length; i++) { re += A.re[i] * B.re[i] + A.im[i] * B.im[i]; im += A.im[i] * B.re[i] - A.re[i] * B.im[i]; } return { re, im }; };
    const xy = tr(X, Y), xx = tr(X, X).re, yy = tr(Y, Y).re;
    const g12 = xy.re, Om = -2 * xy.im, bound = 2 * Math.sqrt(Math.max(0, xx * yy - g12 * g12));
    const cosG = g12 / Math.sqrt(xx * yy);
    return { Om, bound, ratio: Math.abs(Om) / bound, thetaG: Math.acos(Math.max(-1, Math.min(1, cosG))) };
  }
  return { rng, gauss, eigh, principal, cols, deg, randSym, randOrth, makeFeed, feedOp, clusters, tipTail, speed, geodesic, families, cycleAngles, permStats, envelopePoint, matmul, T, opnorm };
})();

// ---------- v2 additions: complex bands (Appendices B and C) and Fisher speed from fidelity (Lemma 5.1) ----------
Object.assign(Core, (() => {
  const C = Core;
  // gauge-fix a complex vector: component idx real and positive, unit norm
  function gauge(re, im, idx) { const a = re[idx], b = im[idx], m = Math.hypot(a, b) || 1, c = a / m, s = -b / m;
    const r2 = re.map((x, i) => x * c - im[i] * s), i2 = im.map((y, i) => re[i] * s + y * c); const n = Math.hypot(...r2, ...i2);
    return { re: r2.map(x => x / n), im: i2.map(x => x / n) }; }
  // lowest eigenvector of a complex Hermitian H = {re, im} through the real embedding [[A,-B],[B,A]]
  function lowestVec(H, idx) { const n = H.re.length, M = Array.from({ length: 2 * n }, () => new Float64Array(2 * n));
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { M[i][j] = H.re[i][j]; M[i + n][j + n] = H.re[i][j]; M[i][j + n] = -H.im[i][j]; M[i + n][j] = H.im[i][j]; }
    const e = C.eigh(M), re = [], im = []; for (let i = 0; i < n; i++) { re.push(e.V[i][0]); im.push(e.V[i + n][0]); }
    if (idx === undefined) { let best = 0; re.forEach((x, i) => { if (Math.hypot(x, im[i]) > Math.hypot(re[best], im[best])) best = i; }); idx = best; }
    return Object.assign(gauge(re, im, idx), { idx, lam: e.w[0] }); }
  const cdot = (u, v) => { let re = 0, im = 0; for (let i = 0; i < u.re.length; i++) { re += u.re[i] * v.re[i] + u.im[i] * v.im[i]; im += u.re[i] * v.im[i] - u.im[i] * v.re[i]; } return { re, im }; };
  // tensor Q_ij = <d_i v_perp | d_j v_perp> of the lowest band of Hf(p), p in R^2
  function bandQ(Hf, p, h = 1e-5) { const v = lowestVec(Hf(p)), idx = v.idx, d = [];
    for (let i = 0; i < 2; i++) { const pp = p.slice(), pm = p.slice(); pp[i] += h; pm[i] -= h; const a = lowestVec(Hf(pp), idx), b = lowestVec(Hf(pm), idx);
      const dv = { re: a.re.map((x, j) => (x - b.re[j]) / (2 * h)), im: a.im.map((x, j) => (x - b.im[j]) / (2 * h)) }; const o = cdot(v, dv);
      d.push({ re: dv.re.map((x, j) => x - (v.re[j] * o.re - v.im[j] * o.im)), im: dv.im.map((x, j) => x - (v.re[j] * o.im + v.im[j] * o.re)) }); }
    return [[cdot(d[0], d[0]), cdot(d[0], d[1])], [cdot(d[1], d[0]), cdot(d[1], d[1])]]; }
  // one random pair of tangent directions V = a.d, W = b.d at the point; returns the same fields as envelopePoint
  function pairFromQ(Q, a, b) { const q = (x, y) => { let re = 0, im = 0; for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) { re += x[i] * y[j] * Q[i][j].re; im += x[i] * y[j] * Q[i][j].im; } return { re, im }; };
    const xy = q(a, b), xx = q(a, a).re, yy = q(b, b).re, g12 = xy.re, Om = -2 * xy.im, bound = 2 * Math.sqrt(Math.max(0, xx * yy - g12 * g12));
    return { Om, bound, ratio: Math.abs(Om) / bound, thetaG: Math.acos(Math.max(-1, Math.min(1, g12 / Math.sqrt(xx * yy)))) }; }
  const r2 = Math.SQRT2;
  const SIG = { x: { re: [[0, 1], [1, 0]], im: [[0, 0], [0, 0]] }, y: { re: [[0, 0], [0, 0]], im: [[0, -1], [1, 0]] }, z: { re: [[1, 0], [0, -1]], im: [[0, 0], [0, 0]] } };
  const SPN = { x: { re: [[0, 1 / r2, 0], [1 / r2, 0, 1 / r2], [0, 1 / r2, 0]], im: [[0, 0, 0], [0, 0, 0], [0, 0, 0]] },
                y: { re: [[0, 0, 0], [0, 0, 0], [0, 0, 0]], im: [[0, -1 / r2, 0], [1 / r2, 0, -1 / r2], [0, 1 / r2, 0]] },
                z: { re: [[1, 0, 0], [0, 0, 0], [0, 0, -1]], im: [[0, 0, 0], [0, 0, 0], [0, 0, 0]] } };
  const comb = (cs, Ms) => { const n = Ms[0].re.length, out = { re: [], im: [] }; for (let i = 0; i < n; i++) { out.re.push(new Float64Array(n)); out.im.push(new Float64Array(n));
    for (let j = 0; j < n; j++) Ms.forEach((M, k) => { out.re[i][j] += cs[k] * M.re[i][j]; out.im[i][j] += cs[k] * M.im[i][j]; }); } return out; };
  const dirH = ops => p => comb([Math.sin(p[0]) * Math.cos(p[1]), Math.sin(p[0]) * Math.sin(p[1]), Math.cos(p[0])], [ops.x, ops.y, ops.z]);
  function randHerm(n, r) { const re = [], im = []; for (let i = 0; i < n; i++) { re.push(new Float64Array(n)); im.push(new Float64Array(n)); }
    for (let i = 0; i < n; i++) for (let j = i; j < n; j++) { const a = C.gauss(r), b = i === j ? 0 : C.gauss(r); re[i][j] = re[j][i] = a; im[i][j] = b; im[j][i] = -b; } return { re, im }; }
  const gr = C.rng(404), G0 = randHerm(3, gr), G1 = randHerm(3, gr), G2 = randHerm(3, gr);
  const bandFamilies = {
    two: { label: 'two levels (d·σ; the polariton case)', H: dirH(SIG), point: r => [0.3 + 2.5 * r(), 2 * Math.PI * r()] },
    spin1: { label: 'spin 1, lowest band (d·S)', H: dirH(SPN), point: r => [0.3 + 2.5 * r(), 2 * Math.PI * r()] },
    generic: { label: 'generic three levels (H₀ + aH₁ + bH₂)', H: p => comb([1, p[0], p[1]], [G0, G1, G2]), point: r => [C.gauss(r) * 0.5, C.gauss(r) * 0.5] } };
  function bandCloud(fam, n, seed) { const r = C.rng(seed), F = bandFamilies[fam], pts = [];
    for (let t = 0; t < n; t++) { const Q = bandQ(F.H, F.point(r)); pts.push(pairFromQ(Q, [C.gauss(r), C.gauss(r)], [C.gauss(r), C.gauss(r)])); } return pts; }
  // Fisher speed of rho = Pi_k / k from the Bures fidelity of neighbouring projectors: F_Q = 8 (1 - (1/k) Σ cos θ_j) / ds²
  function fisherSpeed(F, inner, outer, s, k, ds = 1e-4) { const ids = [...Array(k).keys()];
    const a = C.eigh(C.feedOp(F, inner, outer, s - ds / 2)), b = C.eigh(C.feedOp(F, inner, outer, s + ds / 2));
    const th = C.principal(C.cols(a.V, ids), C.cols(b.V, ids)); const fq = 8 * (1 - th.reduce((x, t) => x + Math.cos(t), 0) / k) / (ds * ds);
    return { fq, v: Math.sqrt(k * fq / 4) }; }
  return { lowestVec, bandQ, pairFromQ, bandFamilies, bandCloud, fisherSpeed };
})());
if (typeof module !== 'undefined') module.exports = Core;
