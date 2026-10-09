// SCRIPT: GRASSMANN-MIXER-APP-V2
(() => {
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const COL = { ink: css('--ink'), blade: css('--blade'), stem: css('--stem'), warm: css('--warm'), vio: css('--vio'), rose: css('--rose'), rule: css('--rule'), bed: css('--bed'), ground: css('--ground') };
const CLUSTER = [COL.blade, COL.warm, COL.vio, COL.rose, '#6fe0b0', '#f2a65a'];
const $ = id => document.getElementById(id);
const S = { inner: 0.08, outer: 1.6, s: 0, noise: 0.004, thr: 0.06, k: 3, fam: 'polar', launch: 35, alpha: 90, mix: 0, tt: 'live', pp: 0, cone: 60, envfam: 'random', film: false, mode: 'mix', focus: 'levels', playing: false, t0: 0 };
const F = Core.makeFeed(11);
const fmt = (x, d = 3) => Number.isFinite(x) ? x.toFixed(d) : '–';

// ---------- cached computations (redone only when their inputs change) ----------
let cache = {};
function levelsTrack() {
  const key = [S.inner, S.outer].join();
  if (cache.levKey === key) return cache.lev;
  const ss = [], ws = [];
  for (let i = 0; i <= 150; i++) { const s = -1.5 + 3 * i / 150; ss.push(s); ws.push(Core.eigh(Core.feedOp(F, S.inner, S.outer, s)).w); }
  cache.levKey = key; cache.lev = { ss, ws }; return cache.lev;
}
function speedTrack() {
  const key = [S.inner, S.outer, S.k].join();
  if (cache.spKey === key) return cache.sp;
  const ss = [], v = [], c = [], fv = [];
  for (let i = 0; i <= 120; i++) { const s = -1.5 + 3 * i / 120; const r = Core.speed(F, S.inner, S.outer, s, S.k); ss.push(s); v.push(r.v); c.push(r.ceil); fv.push(Core.fisherSpeed(F, S.inner, S.outer, s, S.k).v); }
  cache.spKey = key; cache.sp = { ss, v, c, fv }; return cache.sp;
}
function tipTailLive() {
  const key = [S.inner, S.outer, S.s, S.noise].join();
  if (cache.ttKey === key) return cache.tt;
  const r = Core.rng(2024), A = Core.feedOp(F, S.inner, S.outer, S.s), runs = 9, acc = [];
  for (let t = 0; t < runs; t++) acc.push(Core.tipTail(A, S.noise, r));
  const rows = acc[0].map((row, i) => { const med = f => acc.map(a => a[i][f]).sort((x, y) => x - y)[runs >> 1]; return { k: row.k, vec: med('vec'), sub: med('sub'), gap: row.gap }; });
  cache.ttKey = key; cache.tt = rows; return rows;
}
function corridor() {
  const key = [S.outer, S.noise].join();
  if (cache.coKey === key) return cache.co;
  const gaps = [], tips = [], tails = []; const r = Core.rng(77);
  for (let i = 0; i <= 10; i++) {
    const g = Math.pow(10, -2.2 + 2.2 * i / 10); const a = [], b = [];
    for (let t = 0; t < 14; t++) { const rows = Core.tipTail(Core.feedOp(F, g, S.outer, 0), S.noise, r); a.push(rows[1].vec); b.push(rows[2].sub); }
    const med = x => x.sort((p, q) => p - q)[x.length >> 1];
    gaps.push(g); tips.push(med(a)); tails.push(med(b));
  }
  cache.coKey = key; cache.co = { gaps, tips, tails }; return cache.co;
}
function geodesics() {
  const key = [S.fam, S.launch].join();
  if (cache.geKey === key) return cache.ge;
  const f = Core.families[S.fam], x0 = S.fam === 'fold' ? -0.9 : 0;
  const fan = [20, 35, 50, 65, 80].map(th => Core.geodesic(f.a, f.b, x0, f.y0, th, S.fam === 'narrow' ? 4 : 2.2, 1500, S.fam === 'fold' ? -1 : -1));
  const main = Core.geodesic(f.a, f.b, x0, f.y0, S.launch, S.fam === 'narrow' ? 4 : 2.2, 3000);
  let ymin = Math.min(...main.map(p => p[1]));
  let pred = null;
  if (S.fam !== 'fold') { const C = f.a(0, f.y0) * Math.sin(S.launch * Math.PI / 180) ** 2; let lo = S.fam === 'polar' ? 0 : -3, hi = f.y0;
    if (S.fam === 'narrow') { lo = -3; hi = 0; }
    for (let i = 0; i < 80; i++) { const m = (lo + hi) / 2; (f.a(0, m) - C) * (f.a(0, hi) - C) > 0 ? hi = m : lo = m; } pred = (lo + hi) / 2; }
  cache.geKey = key; cache.ge = { fan, main, ymin, pred, x0 }; return cache.ge;
}
let envPts = [];
function envelopeCloud() {
  const key = [S.alpha, S.mix].join();
  if (cache.enKey === key) return cache.en;
  const r = Core.rng(5), pts = [];
  for (let t = 0; t < 160; t++) pts.push(Core.envelopePoint(r, 2, 4, r() * Math.PI, r() * 1.2));
  const cur = Core.envelopePoint(Core.rng(31), 2, 4, S.alpha * Math.PI / 180, S.mix);
  cache.enKey = key; cache.en = { pts, cur }; return cache.en;
}
function bandCloud(fam) { const key = 'band_' + fam; if (cache[key]) return cache[key]; cache[key] = Core.bandCloud(fam, 140, 17); return cache[key]; }
const PAIRS = [ { n: 20, A: [2, 2, 2, 2, 2], B: [10], la: 'five swaps', lb: 'one 10-cycle' },
  { n: 12, A: [2, 2, 2], B: [6], la: 'three swaps', lb: 'one 6-cycle' },
  { n: 12, A: [3, 3], B: [6], la: 'two 3-cycles', lb: 'one 6-cycle' },
  { n: 4, A: [2, 2], B: [2, 2], la: '(01)(23)', lb: '(02)(13)' } ];

// ---------- drawing helpers ----------
function setup(cv) {
  const dpr = Math.min(2, window.devicePixelRatio || 1), w = cv.clientWidth, h = cv.clientHeight;
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
  const c = cv.getContext('2d'); c.setTransform(dpr, 0, 0, dpr, 0, 0); c.clearRect(0, 0, w, h); c.fillStyle = COL.bed; c.fillRect(0, 0, w, h);
  return { c, w, h };
}
function text(c, s, x, y, col = COL.stem, size = 11, align = 'left', font = 'mono') {
  c.fillStyle = col; c.font = `${size}px ${font === 'mono' ? "'JetBrains Mono', monospace" : "'Newsreader', Georgia, serif"}`; c.textAlign = align; c.textBaseline = 'alphabetic'; c.fillText(s, x, y);
}
function axes(c, x0, y0, x1, y1) { c.strokeStyle = COL.rule; c.lineWidth = 1; c.beginPath(); c.moveTo(x0, y0); c.lineTo(x0, y1); c.moveTo(x0, y1); c.lineTo(x1, y1); c.stroke(); }
function line(c, pts, col, lw = 2, dash = null) { c.strokeStyle = col; c.lineWidth = lw; c.setLineDash(dash || []); c.beginPath(); pts.forEach((p, i) => i ? c.lineTo(p[0], p[1]) : c.moveTo(p[0], p[1])); c.stroke(); c.setLineDash([]); }
function dot(c, x, y, r, col) { c.fillStyle = col; c.beginPath(); c.arc(x, y, r, 0, 2 * Math.PI); c.fill(); }

// ---------- the eyes ----------
const EYES = {
  levels: { name: 'Levels and colours', sec: 'Sections 3, 10.2, 12.1',
    short: 'Eigenvalues of the feed as the sweep moves. Levels closer than the fuse threshold sit in one violet band: together they are one subspace.',
    long: 'Each curve is one eigenvalue of the six-level operator as s sweeps. Where two curves approach and turn away, the eigenvectors trade places. Levels that sit closer than the threshold fuse into one cluster, drawn as a violet band around them, because only the cluster, not its members, is held in place. Move the inner gap toward zero and watch the bottom three fuse.',
    draw(cv) { const { c, w, h } = setup(cv); const L = levelsTrack(); const pad = 34; const ymin = Math.min(...L.ws.flat()), ymax = Math.max(...L.ws.flat());
      const X = s => pad + (s + 1.5) / 3 * (w - pad - 10), Y = v => h - 22 - (v - ymin) / (ymax - ymin || 1) * (h - 40);
      axes(c, pad, 8, w - 8, h - 22);
      for (let i = 0; i < L.ss.length - 1; i++) { const cl = Core.clusters(L.ws[i], S.thr); c.fillStyle = 'rgba(169,155,255,0.22)';
        cl.filter(g => g.length > 1).forEach(g => { const lo = Y(L.ws[i][g[0]]) + 4, hi = Y(L.ws[i][g[g.length - 1]]) - 4; c.fillRect(X(L.ss[i]), hi, X(L.ss[i + 1]) - X(L.ss[i]) + 0.5, lo - hi); }); }
      for (let j = 0; j < L.ws[0].length; j++) line(c, L.ss.map((s, i) => [X(s), Y(L.ws[i][j])]), CLUSTER[j % 6], 2);
      c.strokeStyle = COL.ink; c.globalAlpha = 0.5; c.beginPath(); c.moveTo(X(S.s), 8); c.lineTo(X(S.s), h - 22); c.stroke(); c.globalAlpha = 1;
      text(c, 's', w - 14, h - 8);
      const w0 = Core.eigh(Core.feedOp(F, S.inner, S.outer, S.s)).w, cl = Core.clusters(w0, S.thr);
      text(c, `${cl.length} subspaces at s = ${fmt(S.s, 2)}: ` + cl.map(g => g.length).join(' + '), pad + 6, h - 8, COL.ink); } },
  grass: { name: 'The Grassmannian moving', sec: 'Sections 2, 5.2, 8.2',
    short: 'Principal angles between the rank-k subspace at s and at s = 0, and its speed under the gap ceiling.',
    long: 'Left: the moving subspace seen from its starting position. Each spoke is one principal angle between the rank-k eigenspace at s and the same eigenspace at s = 0; the Hilbert–Schmidt distance is the sum of their squared sines. Right: the speed of the subspace in the intrinsic metric along the whole sweep (blue) and the ceiling that the drive and the gap allow (amber, Theorem 5.2). The speed rises toward the ceiling where the gap is narrowest. The dashed line is the Fisher speed of ρ = Π/k computed independently, from the Bures fidelity of neighbouring projectors; it lies on the metric speed, which is Lemma 5.1, F_Q = 4g/k.',
    draw(cv) { const { c, w, h } = setup(cv); const k = S.k;
      const e0 = Core.eigh(Core.feedOp(F, S.inner, S.outer, 0)), e1 = Core.eigh(Core.feedOp(F, S.inner, S.outer, S.s));
      const ids = [...Array(k).keys()], th = Core.principal(Core.cols(e1.V, ids), Core.cols(e0.V, ids));
      const R = Math.min(w * 0.34, h - 76), cx = 18, cy = h - 30;
      c.strokeStyle = COL.rule; c.beginPath(); c.arc(cx, cy, R, -Math.PI / 2, 0); c.stroke();
      th.forEach((t, i) => line(c, [[cx, cy], [cx + R * Math.cos(t), cy - R * Math.sin(t)]], CLUSTER[i % 6], 2.4));
      text(c, 'angles to s = 0', cx, 16, COL.stem, 10); text(c, th.map(t => Core.deg(t).toFixed(1) + '°').join(' '), cx, 30, COL.ink, 10);
      text(c, 'Σ sin²θ = ' + fmt(th.reduce((a, t) => a + Math.sin(t) ** 2, 0)), cx, cy + 16, COL.ink);
      const sp = speedTrack(), x0 = w * 0.5, x1 = w - 10, y0 = 54, y1 = h - 26; axes(c, x0, y0, x1, y1);
      const vmax = Math.max(...sp.c) * 1.05, X = s => x0 + (s + 1.5) / 3 * (x1 - x0), Y = v => y1 - Math.log10(Math.max(v, vmax / 300) / (vmax / 300)) / Math.log10(300) * (y1 - y0);
      line(c, sp.ss.map((s, i) => [X(s), Y(sp.c[i])]), COL.warm, 2); line(c, sp.ss.map((s, i) => [X(s), Y(sp.v[i])]), COL.blade, 3.2); line(c, sp.ss.map((s, i) => [X(s), Y(sp.fv[i])]), COL.rose, 1.6, [5, 4]);
      const now = Core.speed(F, S.inner, S.outer, S.s, k); dot(c, X(S.s), Y(now.v), 4, COL.ink);
      text(c, 'ceiling √k‖∂Σ‖/gap', x0, 16, COL.warm, 10); text(c, 'speed √g', x0, 30, COL.blade, 10); text(c, 'Fisher speed √(kF_Q/4), dashed', x0, 44, COL.rose, 10); text(c, 's', x1 - 6, h - 10); } },
  bloch: { name: 'Bloch sphere: Gr(1, ℂ²)', sec: 'Sections 2, 6, Appendix B',
    short: 'The smallest Grassmannian is a sphere. A loop of the lower-band projector, its metric and the curvature flux it encloses.',
    long: 'For two levels the Grassmannian of lines is the Bloch sphere, and every projector is a point on it. The lower band of d·σ traces a loop of opening θ as s sweeps. The metric is a quarter of the round metric (g_θθ = 1/4) and the curvature is half the area element, so the Berry phase of the loop is half the solid angle it encloses, and the Chern number over the whole sphere is 1. Two levels always sit on the curvature bound: every pair of directions is extremal.',
    draw(cv, tsec) { const { c, w, h } = setup(cv); const R = Math.min(w * 0.3, h * 0.4), cx = w * 0.3, cy = h * 0.52;
      const rot = (tsec || 0) * 0.25, tilt = 0.45, th = S.cone * Math.PI / 180;
      const P = (t, p) => { const x = Math.sin(t) * Math.cos(p + rot), y = Math.sin(t) * Math.sin(p + rot), z = Math.cos(t); const y2 = y * Math.cos(tilt) - z * Math.sin(tilt), z2 = y * Math.sin(tilt) + z * Math.cos(tilt); return [cx + R * x, cy - R * z2, y2]; };
      const grad = c.createRadialGradient(cx - R * 0.3, cy - R * 0.3, R * 0.1, cx, cy, R); grad.addColorStop(0, '#13203a'); grad.addColorStop(1, '#0a1020');
      c.fillStyle = grad; c.beginPath(); c.arc(cx, cy, R, 0, 2 * Math.PI); c.fill(); c.strokeStyle = COL.rule; c.stroke();
      for (let la = 30; la < 180; la += 30) { const pts = []; for (let p = 0; p <= 64; p++) pts.push(P(la * Math.PI / 180, 2 * Math.PI * p / 64)); c.strokeStyle = COL.rule; c.lineWidth = 1; c.beginPath(); pts.forEach((q, i) => { if (q[2] < 0) { c.moveTo(q[0], q[1]); return; } i ? c.lineTo(q[0], q[1]) : c.moveTo(q[0], q[1]); }); c.stroke(); }
      c.fillStyle = 'rgba(169,155,255,0.18)'; c.beginPath(); for (let p = 0; p <= 64; p++) { const q = P(th, 2 * Math.PI * p / 64); p ? c.lineTo(q[0], q[1]) : c.moveTo(q[0], q[1]); } const np = P(0, 0); c.lineTo(np[0], np[1]); c.fill();
      const loop = []; for (let p = 0; p <= 96; p++) loop.push(P(th, 2 * Math.PI * p / 96)); line(c, loop.map(q => [q[0], q[1]]), COL.vio, 2);
      const phi = (S.s + 1.5) / 3 * 2 * Math.PI, q = P(th, phi); dot(c, q[0], q[1], 5, COL.blade);
      const solid = 2 * Math.PI * (1 - Math.cos(th));
      const tx = w * 0.64, rows = [['opening θ', S.cone + '°'], ['solid angle', fmt(solid, 3) + ' sr'], ['Berry phase', fmt(solid / 2, 3) + ' rad'], ['  = ½ solid angle', ''], ['loop length', fmt(Math.PI * Math.sin(th), 3)], ['Chern number', '1']];
      rows.forEach(([a, b], i) => { text(c, a, tx, 24 + i * 26, COL.stem, 10); if (b) text(c, b, tx, 37 + i * 26, COL.ink, 11); }); } },
  tiptail: { name: 'Tip and tail', sec: 'Sections 10.2, 12.3',
    short: 'Two independent halves vote on every direction. Single vectors (rose) against the rank-k subspace (blue), rank by rank.',
    long: 'Two independent noisy copies of the feed are diagonalised separately, and each rank is compared between them. A single eigenvector (rose) is held only by its nearest neighbour, so inside a cluster it swings. The subspace of the first k levels (blue) is held by the gap above it, so at the cut that closes a cluster it holds still. Choose the condensate in the settings to see the measured halves of the NIST images: the same law, in data.',
    draw(cv) { const { c, w, h } = setup(cv); let rows, label;
      if (S.tt === 'live') { rows = tipTailLive(); label = 'live feed, median of 9 noise draws'; } else { rows = BEC[S.tt].map(r => ({ k: r[0], vec: r[1], sub: r[2], gap: r[3] })); label = 'condensate: ' + S.tt; }
      const pad = 30, n = rows.length, bw = (w - pad - 10) / n, vmax = Math.max(10, ...rows.map(r => Math.max(r.vec, r.sub))) * 1.08, Y = v => h - 22 - v / vmax * (h - 44);
      axes(c, pad, 10, w - 8, h - 22);
      rows.forEach((r, i) => { const x = pad + i * bw + bw * 0.12;
        if (r.gap < 0.10 && S.tt !== 'live' || (S.tt === 'live' && r.gap < S.thr)) { c.fillStyle = 'rgba(169,155,255,0.10)'; c.fillRect(pad + i * bw, 10, bw, h - 32); }
        c.fillStyle = COL.rose; c.fillRect(x, Y(r.vec), bw * 0.36, h - 22 - Y(r.vec)); c.fillStyle = COL.blade; c.fillRect(x + bw * 0.38, Y(r.sub), bw * 0.36, h - 22 - Y(r.sub));
        text(c, String(r.k), pad + i * bw + bw / 2, h - 8, COL.stem, 10, 'center'); });
      text(c, label, pad + 6, 22, COL.ink); text(c, 'shaded: no resolved gap at that cut', pad + 6, 36); text(c, 'deg', 4, 18); } },
  corridor: { name: 'The corridor', sec: 'Section 12.1',
    short: 'Tip and tail against the inner gap. The vector follows the inner gap; the cluster stays inside the outer gap’s corridor.',
    long: 'The cluster of three levels is perturbed at fixed noise while its inner gap shrinks. One vector inside the cluster turns further and further, following noise over inner gap. The cluster subspace stays put, held by the outer gap. The marker is the current inner gap of the feed.',
    draw(cv) { const { c, w, h } = setup(cv); const co = corridor(); const pad = 34, x0 = pad, x1 = w - 10, y0 = 12, y1 = h - 24;
      axes(c, x0, y0, x1, y1); const all = co.tips.concat(co.tails).filter(v => v > 0); const lo = Math.log10(Math.min(...all) / 2), hi = Math.log10(Math.max(...all) * 2);
      const X = g => x0 + (Math.log10(g) + 2.2) / 2.2 * (x1 - x0), Y = v => y1 - (Math.log10(Math.max(v, 1e-6)) - lo) / (hi - lo) * (y1 - y0);
      line(c, co.gaps.map((g, i) => [X(g), Y(co.tips[i])]), COL.rose, 2.2); line(c, co.gaps.map((g, i) => [X(g), Y(co.tails[i])]), COL.blade, 2.2);
      co.gaps.forEach((g, i) => { dot(c, X(g), Y(co.tips[i]), 3, COL.rose); dot(c, X(g), Y(co.tails[i]), 3, COL.blade); });
      c.strokeStyle = COL.ink; c.globalAlpha = 0.5; c.beginPath(); c.moveTo(X(Math.max(0.0063, Math.min(1, S.inner))), y0); c.lineTo(X(Math.max(0.0063, Math.min(1, S.inner))), y1); c.stroke(); c.globalAlpha = 1;
      text(c, 'tip: one vector in the cluster', x1 - 4, y0 + 12, COL.rose, 10, 'right'); text(c, 'tail: the cluster subspace', x1 - 4, y0 + 26, COL.blade, 10, 'right');
      text(c, 'inner gap (log)', x1 - 4, h - 8, COL.stem, 11, 'right'); text(c, 'deg (log)', 4, y0 + 2); } },
  strata: { name: 'The hypersurface', sec: 'Section 8',
    short: 'A parameter plane coloured by the refractive index √g∥, with geodesics launched at the degenerate line.',
    long: 'Brightness is the refractive index √g∥, the rate at which the occupied subspace rotates along the interface. Geodesics of the intrinsic metric obey Clairaut’s law g∥ sin²θ = constant. Where the subspace stops rotating (polar stratum) they turn back by total internal reflection before the line. Where the gap narrows the index is largest and they cross and are guided. At a fold they cross freely. The highlighted path uses the launch angle in the settings, and the readout compares its turning point with the conserved-level prediction.',
    draw(cv) { const { c, w, h } = setup(cv); const f = Core.families[S.fam], [vx0, vx1, vy0, vy1] = f.view, G = geodesics();
      const X = x => (x - vx0) / (vx1 - vx0) * w, Y = y => h - (y - vy0) / (vy1 - vy0) * h;
      const nx = 64, ny = 40, img = c.createImageData(nx, ny); let amax = 0; const vals = [];
      for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const x = vx0 + (i + 0.5) / nx * (vx1 - vx0), y = vy1 - (j + 0.5) / ny * (vy1 - vy0); const v = Math.sqrt(Math.max(0, f.a(x, y))); vals.push(v); amax = Math.max(amax, v); }
      vals.forEach((v, idx) => { const t = Math.pow(v / amax, 0.8); img.data[idx * 4] = 12 + 40 * t; img.data[idx * 4 + 1] = 18 + 120 * t; img.data[idx * 4 + 2] = 32 + 190 * t; img.data[idx * 4 + 3] = 255; });
      const off = document.createElement('canvas'); off.width = nx; off.height = ny; off.getContext('2d').putImageData(img, 0, 0); c.imageSmoothingEnabled = true; c.drawImage(off, 0, 0, w, h);
      c.strokeStyle = COL.rose; c.lineWidth = 2; c.beginPath(); if (f.line === 'x0') { c.moveTo(X(0), 0); c.lineTo(X(0), h); } else { c.moveTo(0, Y(0)); c.lineTo(w, Y(0)); } c.stroke();
      G.fan.forEach(p => line(c, p.map(q => [X(q[0]), Y(q[1])]), 'rgba(224,184,92,0.55)', 1.3));
      line(c, G.main.map(q => [X(q[0]), Y(q[1])]), COL.warm, 2.6);
      if (G.pred !== null) { c.setLineDash([4, 4]); c.strokeStyle = COL.ink; c.globalAlpha = 0.6; c.beginPath(); c.moveTo(0, Y(G.pred)); c.lineTo(w, Y(G.pred)); c.stroke(); c.setLineDash([]); c.globalAlpha = 1; }
      c.font = "10px 'JetBrains Mono', monospace"; const words = f.note.split(' '), lines = [''], maxw = w - 28; words.forEach(wd => { const t = lines[lines.length - 1] ? lines[lines.length - 1] + ' ' + wd : wd; if (c.measureText(t).width > maxw) lines.push(wd); else lines[lines.length - 1] = t; });
      c.fillStyle = 'rgba(7,10,18,0.78)'; c.fillRect(6, 6, w - 12, 22 + 13 * lines.length); text(c, f.label, 12, 20, COL.ink, 11); lines.forEach((l, i) => text(c, l, 12, 34 + 13 * i, COL.stem, 10)); } },
  envelope: { name: 'Curvature inside the metric', sec: 'Section 6, Appendices B–C',
    short: 'Curvature per unit of capacity against the geometric angle. Nothing can rise above sin θG; the partner iX sits on it.',
    long: 'For a pair of tangent directions at a projector, the curvature divided by 2√(C(V)C(W)) can never exceed sin θG, the sine of their angle in the intrinsic metric: this is the sharp bound (6.1) in the reporting form (10.2). Grey points are random pairs. The bright point is the partner you set: turning the block X by α moves it along the envelope, and α = 90° (the partner iX) sits exactly on it; mixing in an unrelated block drops it inside. The band switch draws real bands: every two-level pair lands on the curve (the polariton measurement is this case), the lowest spin-1 band lands on it too because it is a spin coherent state with Kähler geometry, and a generic three-level band falls strictly inside (Appendices B and C).',
    draw(cv) { const { c, w, h } = setup(cv); const E = envelopeCloud(); const pad = 34, x0 = pad, x1 = w - 10, y0 = 12, y1 = h - 24; axes(c, x0, y0, x1, y1);
      const X = t => x0 + t / Math.PI * (x1 - x0), Y = v => y1 - v * (y1 - y0);
      const env = []; for (let i = 0; i <= 100; i++) { const t = Math.PI * i / 100; env.push([X(t), Y(Math.sin(t))]); } line(c, env, COL.warm, 2);
      E.pts.forEach(p => dot(c, X(p.thetaG), Y(p.ratio * Math.sin(p.thetaG)), 2.2, 'rgba(112,136,168,0.7)'));
      const cu = E.cur;
      if (S.envfam !== 'random') { const fc = { two: COL.blade, spin1: '#6fe0b0', generic: COL.rose }[S.envfam]; const B = bandCloud(S.envfam);
        B.forEach(p => dot(c, X(p.thetaG), Y(p.ratio * Math.sin(p.thetaG)), 3, fc)); const rs = B.map(p => p.ratio).sort((a, b) => a - b);
        text(c, Core.bandFamilies[S.envfam].label, x0 + 8, y0 + 38, fc, 10); text(c, `|Ω| / bound: ${fmt(rs[0] * 100, 1)}% to ${fmt(rs[rs.length - 1] * 100, 1)}%`, x0 + 8, y0 + 52, fc, 10); }
      dot(c, X(cu.thetaG), Y(cu.ratio * Math.sin(cu.thetaG)), 5.5, COL.vio);
      text(c, 'sin θG, the bound', x1 - 4, y0 + 10, COL.warm, 10, 'right'); text(c, 'θG', x1 - 14, h - 8); text(c, '|Ω| / 2√(C C)', x0 + 8, y0 + 10, COL.stem, 10);
      text(c, `partner: |Ω| = ${fmt(cu.ratio * 100, 1)}% of the bound`, x0 + 8, y0 + 24, COL.ink, 10); } },
  vote: { name: 'The vote', sec: 'Sections 10.2, 12.3',
    short: 'Every rank as a dial: the amber needle is one direction, half B against half A; the blue arc is the rank-k subspace.',
    long: 'The same two-halves comparison as the tip and tail bars, drawn as dials. At each rank the grey needle is half A, the amber needle is where half B puts the single direction, and the blue arc is how far the rank-k subspace moved. Inside a near-tie the needle swings wide while the arc stays short at the cut that closes the cluster. The highlighted dial follows the rank slider, and its numbers are written underneath. Uses the same data choice as the tip and tail eye: live from the feed, or the measured condensate halves.',
    draw(cv) { const { c, w, h } = setup(cv); const live = S.tt === 'live';
      const rows = live ? tipTailLive() : BEC[S.tt].map(r => ({ k: r[0], vec: r[1], sub: r[2], gap: r[3] }));
      const cols = live ? rows.length : 6, nr = Math.ceil(rows.length / cols), cw = w / cols, chh = (h - 46) / nr, R = Math.min(cw * 0.38, chh * 0.62);
      rows.forEach((r, i) => { const cx = (i % cols + 0.5) * cw, cy = 14 + Math.floor(i / cols) * chh + chh / 2 + R / 2 - 6, sel = r.k === S.k;
        c.strokeStyle = sel ? COL.ink : COL.rule; c.lineWidth = sel ? 1.6 : 1; c.beginPath(); c.arc(cx, cy, R, Math.PI, 2 * Math.PI); c.stroke();
        const ang = d => Math.PI + Math.min(90, d) / 90 * (Math.PI / 2);
        c.strokeStyle = COL.blade; c.lineWidth = 4; c.beginPath(); c.arc(cx, cy, R - 4, Math.PI, ang(r.sub)); c.stroke();
        line(c, [[cx, cy], [cx - R, cy]], COL.stem, 1.4); const a = ang(r.vec); line(c, [[cx, cy], [cx + R * Math.cos(a), cy + R * Math.sin(a)]], COL.warm, 2.2);
        if (!live && r.gap < 0.10 || live && r.gap < S.thr) { c.fillStyle = 'rgba(169,155,255,0.14)'; c.beginPath(); c.arc(cx, cy, R + 3, Math.PI, 2 * Math.PI); c.fill(); }
        text(c, 'k = ' + r.k, cx, cy + 13, sel ? COL.ink : COL.stem, 10, 'center'); });
      const r = rows.find(x => x.k === S.k) || rows[0];
      text(c, `rank ${r.k}: direction ${fmt(r.vec, 2)}°, subspace ${fmt(r.sub, 2)}°, relative gap ${fmt(r.gap, 3)}`, 8, h - 22, COL.ink, 11);
      text(c, (live ? 'live feed' : 'condensate, ' + S.tt) + '   ·   violet: no resolved gap at that cut', 8, h - 8, COL.stem, 10); } },
  pair: { name: 'The condensate pair', sec: 'Section 12.3',
    short: 'Measured eigen-images 2 and 3 of the condensate: the sine and cosine of one moving stripe pattern. The phase is the tip; the plane is the tail.',
    long: 'These are the real second and third eigen-images of the NIST dark-soliton absorption images (each carries about 6.7% of the variance; their overlap is 0.012 and their stripes are offset by a quarter period, 90.6°). The large image is cos φ · image 2 + sin φ · image 3: as φ turns the stripe slides, which is the soliton swinging in the trap (2.70 Hz in hold time, measured). Which single image comes out on top in one half of the data depends on where in the swing that half sits, so the phase, the tip, is free; the plane of the two images is the motion itself and holds. On the odd/even split the third eigen-image turns 40.7° between halves while the subspace that closes the pair turns 5.2°. The sweep slider sets φ; run the sweep to watch it slide.',
    draw(cv) { const { c, w, h } = setup(cv); const E = EIG, phi = (S.s + 1.5) / 3 * 2 * Math.PI, cs = Math.cos(phi), sn = Math.sin(phi);
      const img = (arr, x, y, sc) => { const im = c.createImageData(E.w, E.h); let m = 0; arr.forEach(v => m = Math.max(m, Math.abs(v)));
        arr.forEach((v, i) => { const t = v / (m || 1), p = Math.max(0, t), n = Math.max(0, -t); im.data[4 * i] = 12 + 64 * p + 243 * n; im.data[4 * i + 1] = 18 + 177 * p + 104 * n; im.data[4 * i + 2] = 32 + 223 * p + 101 * n; im.data[4 * i + 3] = 255; });
        const off = document.createElement('canvas'); off.width = E.w; off.height = E.h; off.getContext('2d').putImageData(im, 0, 0); c.imageSmoothingEnabled = true; c.drawImage(off, x, y, E.w * sc, E.h * sc); };
      const small = Math.min((w * 0.22) / E.w, (h * 0.36) / E.h), big = Math.min((w * 0.46) / E.w, (h - 40) / E.h);
      img(E.e2, 8, 22, small); img(E.e3, 8, 30 + E.h * small + 14, small);
      text(c, 'eigen-image 2', 8, 16, COL.stem, 10); text(c, 'eigen-image 3', 8, 30 + E.h * small + 10, COL.stem, 10);
      const mix = E.e2.map((v, i) => cs * v + sn * E.e3[i]); const bx = 16 + E.w * small; img(mix, bx, 22, big);
      text(c, 'cos φ · image 2 + sin φ · image 3', bx, 16, COL.ink, 10);
      const R = Math.min((w - (bx + E.w * big) - 30) / 2, h * 0.3), cx = bx + E.w * big + 14 + R, cy = h * 0.42;
      if (R > 18) { c.strokeStyle = COL.blade; c.lineWidth = 2; c.fillStyle = 'rgba(76,195,255,0.10)'; c.beginPath(); c.arc(cx, cy, R, 0, 2 * Math.PI); c.fill(); c.stroke();
        line(c, [[cx - R, cy], [cx + R, cy]], COL.rule, 1); line(c, [[cx, cy - R], [cx, cy + R]], COL.rule, 1); line(c, [[cx, cy], [cx + R * cs, cy - R * sn]], COL.warm, 2.6); dot(c, cx + R * cs, cy - R * sn, 4, COL.warm);
        text(c, 'the pair plane: tail', cx, cy + R + 16, COL.blade, 10, 'center'); text(c, 'phase φ = ' + (Core.deg(phi) % 360).toFixed(0) + '°: tip', cx, cy + R + 30, COL.warm, 10, 'center'); }
      text(c, 'odd/even, rank 3: single image 40.7°, pair-closing subspace 5.2°', 8, h - 20, COL.stem, 10); text(c, 'soliton swing 2.70 Hz in hold time (measured)', 8, h - 7, COL.stem, 10); } },
  perms: { name: 'Permutations: the card count', sec: 'Section 12.4',
    short: 'Two permutations with the same Hamming distance. Their graph subspaces have the same Σ sin²θ and different geodesic lengths.',
    long: 'The Hamming distance of a permutation is a squared chordal distance between its graph subspace and the graph of the identity: it reads Σ sin²θ over the principal angles. The intrinsic metric adds the geodesic length √Σθ². Equal counts can hide different arrangements: swaps put their angles at 90°, a long cycle spreads them in sine–cosine pairs. The bars are the displacement Gram spectrum, 2 sin²(πj/ℓ), whose equal pairs are the sine and cosine of one cyclic motion. The last pair in the list has identical spectra; there a fixed probe still tells them apart.',
    draw(cv) { const { c, w, h } = setup(cv); const P = PAIRS[S.pp], A = Core.permStats(P.A, P.n), B = Core.permStats(P.B, P.n);
      const half = w / 2, R = Math.min(half * 0.55, h * 0.3);
      [[A, P.la, COL.rose, 0], [B, P.lb, COL.blade, half]].forEach(([st, lab, col, ox]) => { const cx = ox + 14, cy = h * 0.44;
        c.strokeStyle = COL.rule; c.beginPath(); c.arc(cx, cy, R, -Math.PI / 2, 0); c.stroke();
        st.th.filter(t => t > 1e-9).forEach(t => line(c, [[cx, cy], [cx + R * Math.cos(t), cy - R * Math.sin(t)]], col, 2));
        text(c, lab, cx, 16, COL.ink, 12); text(c, `d_H = ${fmt(st.dH, 3)}`, cx, cy + 16, COL.stem, 10); text(c, `Σ sin²θ = ${fmt(st.s2, 3)}`, cx, cy + 29, COL.stem, 10); text(c, `√Σθ² = ${fmt(st.L, 3)}`, cx, cy + 43, col, 11);
        const bars = st.gram.filter(g => g.val > 1e-9), bx = cx, by = h - 8, bwid = Math.min(10, (half - 30) / Math.max(1, bars.length));
        bars.forEach((g, i) => { const hh = g.val / 2 * (h * 0.22); c.fillStyle = col; c.globalAlpha = 0.75; c.fillRect(bx + i * (bwid + 2), by - hh, bwid, hh); c.globalAlpha = 1; }); });
      if (S.pp === 3) text(c, 'same spectrum; probe (e0+e1)/√2 reads m_b = 1 vs 2/3', 14, 30, COL.vio, 10); } }
};
const ORDER = ['levels', 'grass', 'bloch', 'tiptail', 'vote', 'pair', 'corridor', 'strata', 'envelope', 'perms'];
const BEC = JSON.parse(document.getElementById('bec-data').textContent);
const EIG = JSON.parse(document.getElementById('eig-data').textContent);

// ---------- build the page ----------
const grid = $('grid'), list = $('eyelist'), panels = {};
ORDER.forEach(key => { const e = EYES[key];
  const b = document.createElement('button'); b.textContent = e.name; b.id = 'eye-' + key; b.setAttribute('aria-pressed', 'false');
  b.onclick = () => { S.focus = key; setMode('eye'); }; list.appendChild(b);
  const p = document.createElement('section'); p.className = 'panel';
  p.innerHTML = `<div class="head"><h2>${e.name}</h2><span class="sec">${e.sec}</span></div><canvas aria-label="${e.name}" role="img"></canvas><p>${e.short}</p><p class="long">${e.long}</p>`;
  grid.appendChild(p); panels[key] = p; });
function setMode(m) { S.mode = m; $('mode-mix').setAttribute('aria-pressed', m === 'mix'); $('mode-eye').setAttribute('aria-pressed', m === 'eye');
  grid.classList.toggle('single', m === 'eye'); ORDER.forEach(k => { panels[k].hidden = m === 'eye' && k !== S.focus; $('eye-' + k).setAttribute('aria-pressed', m === 'eye' && k === S.focus); }); redraw(); }
$('mode-mix').onclick = () => setMode('mix'); $('mode-eye').onclick = () => setMode('eye');

const sliders = { inner: v => v.toFixed(3), outer: v => v.toFixed(2), s: v => v.toFixed(2), noise: v => v.toFixed(4), thr: v => v.toFixed(2), k: v => String(v), launch: v => v + '°', alpha: v => v + '°', mix: v => v.toFixed(2), cone: v => v + '°' };
Object.entries(sliders).forEach(([id, f]) => { const el = $(id), out = $('o-' + id); const set = () => { S[id] = parseFloat(el.value); out.textContent = f(S[id]); }; set(); el.addEventListener('input', () => { set(); redraw(); }); });
['fam', 'tt', 'pp', 'envfam'].forEach(id => { const el = $(id); const set = () => { S[id] = id === 'pp' ? parseInt(el.value) : el.value; $('o-' + id).textContent = ''; }; set(); el.addEventListener('change', () => { set(); redraw(); }); });

function readout() {
  const e = Core.eigh(Core.feedOp(F, S.inner, S.outer, S.s)), cl = Core.clusters(e.w, S.thr), sp = Core.speed(F, S.inner, S.outer, S.s, S.k);
  const tt = tipTailLive(), G = geodesics(), P = PAIRS[S.pp], A = Core.permStats(P.A, P.n), B = Core.permStats(P.B, P.n), env = envelopeCloud().cur;
  const rows = [ ['eigenvalues', e.w.map(x => x.toFixed(2)).join(' ')], ['subspaces', cl.map(g => g.length).join(' + ')],
    ['gap above rank k', fmt(e.w[S.k] - e.w[S.k - 1], 3)], ['speed / ceiling', fmt(sp.v, 3) + ' / ' + fmt(sp.ceil, 3)],
    ['tip angle, rank 2', fmt(tt[1].vec, 2) + '°'], ['tail angle, rank 3', fmt(tt[2].sub, 2) + '°'],
    ['turning point', G.pred === null ? 'none (crosses)' : fmt(G.ymin, 4) + ' vs ' + fmt(G.pred, 4)],
    ['Fisher / metric speed', fmt(Core.fisherSpeed(F, S.inner, S.outer, S.s, S.k).v, 3) + ' / ' + fmt(sp.v, 3)],
    ['|Ω| of bound (partner)', fmt(env.ratio * 100, 1) + '%'], ['pair phase φ', ((S.s + 1.5) / 3 * 360).toFixed(0) + '°'], ['geodesic A / B', fmt(A.L, 3) + ' / ' + fmt(B.L, 3)] ];
  $('readout').innerHTML = rows.map(([a, b]) => `<dt>${a}</dt><dd>${b}</dd>`).join('');
}
function redraw(tsec) { ORDER.forEach(k => { if (!panels[k].hidden) EYES[k].draw(panels[k].querySelector('canvas'), tsec); }); readout(); }

// ---------- sweep animation ----------
const reduce = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
let last = 0;
function frame(ts) { if (!S.playing) return; const t = ts / 1000; if (t - last > 0.06) { last = t; S.s = Math.sin(t * 0.35) * 1.5; $('s').value = S.s; $('o-s').textContent = S.s.toFixed(2); redraw(t); } requestAnimationFrame(frame); }
$('play').onclick = () => { S.playing = !S.playing; $('play').setAttribute('aria-pressed', S.playing); $('play').textContent = S.playing ? 'Pause sweep' : 'Run sweep'; if (S.playing) requestAnimationFrame(frame); };
window.addEventListener('resize', () => redraw());
// ---------- film order: each eye in turn, sweep running, for recording ----------
let filmTimer = null;
function filmStep(i) { S.focus = ORDER[i % ORDER.length]; setMode('eye'); filmTimer = setTimeout(() => filmStep(i + 1), 12000); }
$('film').onclick = () => { S.film = !S.film; $('film').setAttribute('aria-pressed', S.film); $('film').textContent = S.film ? 'Stop film' : 'Run film';
  if (S.film) { if (!S.playing) $('play').click(); filmStep(ORDER.indexOf(S.focus) < 0 ? 0 : ORDER.indexOf(S.focus)); } else { clearTimeout(filmTimer); if (S.playing) $('play').click(); } };
// recording hook: show any eye at any sweep position, with optional settings, from the console or a capture script
window.MIX = { eyes: ORDER.slice(), state: S, show(k) { S.focus = k; setMode('eye'); },
  render(k, s, over) { if (over) Object.assign(S, over); if (s !== undefined) { S.s = s; $('s').value = s; $('o-s').textContent = (+s).toFixed(2); } S.focus = k; setMode('eye'); return true; },
  mix() { setMode('mix'); } };
setMode('mix');
if (!reduce) { /* start still; the reader starts the sweep */ }
})();
