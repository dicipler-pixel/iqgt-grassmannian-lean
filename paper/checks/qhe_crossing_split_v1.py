# SCRIPT: QHE-CROSSING-SPLIT-V1
# Dry run (synthetic files, same layout): planted structure LIVES in both eyes and both splits; pure shot
# noise is KILLED in both. Pure noise alone also gives mean steps of 17-22 deg and jumps up to 89 deg.
# Question: are the big jumps of the leading eigenvector in the heat-engine population trajectory
# real crossings, or ties decided by shot noise?
# Data: the authors' raw single-shot files single_shots_calibration.nc and single_shots_three_cycles.nc
# (Uusnakki et al., quantum heat engine), taken from the data_analysis zip in your Drive.
# Method, the same as the July notebook (Untitled21): 4-component Gaussian mixture on the calibration
# shots, the authors' overlap-correction matrix, 128-point g/e/f/h trajectory, leading eigenvector of
# the 4x4 population covariance over a sliding 5-point window. Two eyes: CORRECTED (correction matrix)
# and RAW (plain mixture labels). New here: every cell's shots are split into odd and even halves, and
# each half gets its own independent trajectory.
#
# Kill conditions, written BEFORE the run:
# K1 REPEAT  A big jump is a step of more than 45 deg. A real crossing is in the system, so it must jump
#            in BOTH halves at the same step (within 1 step). Kill (jumps are noise-decided ties) if fewer
#            than half of the odd-half jumps reappear in the even half, OR if the number of matches is not
#            above the 99th percentile of the shifted null (the even half slid along the trajectory).
# K2 GAP     Reported only, not a kill: the relative gap (l1 - l2)/l1 at the two windows around each repeated
#            jump. A crossing that falls BETWEEN samples need not show a small gap at either sampled window
#            (the dry run showed this: 13 of 15 planted crossings had large gaps on both sides).
# K3 SEAMS   The trajectory is stitched from strokes and cycles. A jump whose window steps across a seam
#            is reported separately: it may be the stitch, not physics.
# Opposite check: the same test on a trajectory with the shot order shuffled before splitting (a different
# but equally independent split). Real crossings must survive it too.
import os, sys, zipfile, json, subprocess
import numpy as np
if os.path.isdir('/content'):
    subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', 'xarray', 'netCDF4'])

OUT = 'QHE_crossing_split_v1'

def find_files():
    need = ['single_shots_calibration.nc', 'single_shots_three_cycles.nc']
    for base in ['.', '/content', '/content/qhe']:
        if all(os.path.exists(os.path.join(base, n)) for n in need):
            return [os.path.join(base, n) for n in need]
    root = '/content/drive/MyDrive'
    if os.path.isdir('/content') and not os.path.isdir(root):
        from google.colab import drive
        drive.mount('/content/drive')
    zips = []
    for dirpath, dirnames, filenames in os.walk(root):
        for f in filenames:
            if f.startswith('data_analysis_for_') and f.endswith('.zip'):
                zips.append(os.path.join(dirpath, f))
    zips.sort(key=os.path.getmtime, reverse=True)
    for zp in zips:
        with zipfile.ZipFile(zp) as z:
            names = z.namelist()
            hits = {n: [m for m in names if m.endswith(n)] for n in need}
            if all(hits[n] for n in need):
                os.makedirs('/content/qhe', exist_ok=True)
                out = []
                for n in need:
                    m = hits[n][0]
                    with z.open(m) as src, open(os.path.join('/content/qhe', n), 'wb') as dst:
                        dst.write(src.read())
                    out.append(os.path.join('/content/qhe', n))
                print('using', os.path.basename(zp))
                return out
    raise SystemExit('could not find the two .nc files or a data_analysis zip that holds them')

def ellipse_count(re, im, mean, cov):
    w, v = np.linalg.eigh(cov)
    u = v[:, 0] / np.linalg.norm(v[:, 0])
    ang = np.degrees(np.arctan(u[1] / u[0]))
    c, s = np.cos(np.radians(180.0 - ang)), np.sin(np.radians(180.0 - ang))
    xc, yc = re - mean[0], im - mean[1]
    xt, yt = xc * c - yc * s, xc * s + yc * c
    a, b = np.sqrt(w[0]), np.sqrt(w[1])
    return int(np.sum(xt ** 2 / a ** 2 + yt ** 2 / b ** 2 <= 1.0))

def main(paths=None, seed=0):
    import xarray as xr
    from sklearn.mixture import GaussianMixture
    rng = np.random.default_rng(seed)
    cal_p, cyc_p = paths or find_files()
    cal = xr.open_dataset(cal_p); cyc = xr.open_dataset(cyc_p)
    R, I = cal['real'].values, cal['imag'].values
    s0 = (R[0, 0, 0, 0, 0, 0, :], I[0, 0, 0, 0, 0, 0, :])
    s1 = (R[1, 0, 0, 0, 0, 0, :], I[1, 0, 0, 0, 0, 0, :])
    s2 = (R[2, 0, 0, 1, 0, 0, :], I[2, 0, 0, 1, 0, 0, :])
    X = np.column_stack([np.concatenate([s0[0], s1[0], s2[0]]), np.concatenate([s0[1], s1[1], s2[1]])])
    X = X[np.all(np.isfinite(X), axis=1)]
    gmm = GaussianMixture(n_components=4, random_state=0, n_init=3).fit(X)
    corr = np.zeros((4, 4))
    for i in range(4):
        pts = rng.multivariate_normal(gmm.means_[i], gmm.covariances_[i], 1000000)
        for j in range(4):
            corr[i, j] = ellipse_count(pts[:, 0], pts[:, 1], gmm.means_[j], gmm.covariances_[j]) / len(pts)

    def pop_corrected(re, im):
        cnt = np.array([ellipse_count(re, im, gmm.means_[j], 0.4 * gmm.covariances_[j]) for j in range(4)], float)
        v = corr @ cnt
        return v / v.sum()

    def pop_raw(re, im):
        lab = gmm.predict(np.column_stack([re, im]))
        v = np.bincount(lab, minlength=4).astype(float)
        return v / v.sum()

    shared = {2: [19, 20, 21, 22, 23, 24, 25, 26, 27, 28], 3: [1, 3, 6, 8, 10, 13, 15, 16, 17, 18],
              4: [0, 19, 20, 21, 22, 23, 24, 25, 26], 5: [1, 2, 4, 5, 7, 8, 9, 11, 12, 14]}
    first = {0: [0], 1: [1, 2, 4, 5, 7, 8, 9, 11, 12, 14]}
    CR, CI = cyc['real'].values, cyc['imag'].values
    cells = []
    for c in range(3):
        d = {**first, **shared} if c == 0 else dict(shared)
        for st, amps in d.items():
            for a in amps:
                cells.append((c, st, a))

    def trajectories(order_seed=None):
        out = {'corrected': {'all': [], 'odd': [], 'even': []}, 'raw': {'all': [], 'odd': [], 'even': []}}
        r2 = np.random.default_rng(order_seed) if order_seed is not None else None
        for (c, st, a) in cells:
            re, im = CR[c, st, a, :], CI[c, st, a, :]
            ok = np.isfinite(re) & np.isfinite(im)
            re, im = re[ok], im[ok]
            idx = np.arange(len(re))
            if r2 is not None:
                idx = r2.permutation(len(re))
            halves = {'all': idx, 'odd': idx[1::2], 'even': idx[0::2]}
            for h, sel in halves.items():
                out['corrected'][h].append(pop_corrected(re[sel], im[sel]))
                out['raw'][h].append(pop_raw(re[sel], im[sel]))
        return {m: {h: np.array(v) for h, v in d.items()} for m, d in out.items()}

    def path(P, window=5):
        vecs, gaps = [], []
        for i in range(len(P) - window):
            w, U = np.linalg.eigh(np.cov(P[i:i + window].T))
            o = np.argsort(w)[::-1]; w, U = w[o], U[:, o]
            vecs.append(U[:, 0]); gaps.append((w[0] - w[1]) / w[0] if w[0] > 0 else 0.0)
        steps = [float(np.degrees(np.arccos(min(1.0, abs(vecs[i] @ vecs[i + 1]))))) for i in range(len(vecs) - 1)]
        gstep = [min(gaps[i], gaps[i + 1]) for i in range(len(vecs) - 1)]
        return np.array(steps), np.array(gstep)

    seam = []
    for i in range(len(cells) - 6):
        seam.append(cells[i + 5][:2] != cells[i + 4][:2])
    seam = np.array(seam)

    def judge(T, tag):
        res = {}
        for m in ('corrected', 'raw'):
            sa, ga = path(T[m]['all']); so, go = path(T[m]['odd']); se, ge = path(T[m]['even'])
            jo = np.where(so > 45)[0]; je = np.where(se > 45)[0]
            match = [int(i) for i in jo if np.any(np.abs(je - i) <= 1)]
            n = len(so)
            null = []
            for sh in range(3, n - 3):
                jes = (je + sh) % n
                null.append(sum(1 for i in jo if np.any(np.abs(jes - i) <= 1)))
            q99 = float(np.percentile(null, 99)) if null else 0.0
            k1 = (len(jo) > 0) and (len(match) >= 0.5 * len(jo)) and (len(match) > q99)
            gap_ok = [i for i in match if go[i] <= 0.5 and ge[i] <= 0.5]
            seams = [i for i in match if seam[min(i, len(seam) - 1)]]
            res[m] = dict(steps_all_mean=float(sa.mean()), jumps_all=int(np.sum(sa > 45)), top_all=sorted(np.round(sa, 1).tolist())[-4:],
                          jumps_odd=len(jo), jumps_even=len(je), matched=len(match), matched_steps=match,
                          null_q99=q99, null_mean=float(np.mean(null)) if null else 0.0,
                          K1_repeat='LIVES' if k1 else 'KILLED', K2_matched_at_small_gap=len(gap_ok),
                          K3_matched_on_seams=len(seams), step_odd_vs_even_corr=float(np.corrcoef(so, se)[0, 1]))
            r = res[m]
            sm, ja, ta, nm, k1s, cc = r['steps_all_mean'], r['jumps_all'], r['top_all'], r['null_mean'], r['K1_repeat'], r['step_odd_vs_even_corr']
            print(f'[{tag}] {m:9s} full-data mean step {sm:.1f} deg, jumps over 45: {ja}, largest {ta}')
            print(f'[{tag}] {m:9s} odd jumps {len(jo)}, even jumps {len(je)}, repeated in both {len(match)} (shifted-null mean {nm:.2f}, 99th pct {q99:.1f})')
            print(f'[{tag}] {m:9s} K1 repeat: {k1s} | K2 repeated jumps at a small gap in both halves: {len(gap_ok)} of {len(match)} | K3 repeated jumps on a seam: {len(seams)} | step-size correlation odd vs even {cc:+.2f}')
        return res

    print('=== QHE-CROSSING-SPLIT-V1 ===')
    print('cells:', len(cells), '| shots per cell:', CR.shape[-1])
    R1 = judge(trajectories(None), 'odd/even')
    R2 = judge(trajectories(12345), 'shuffled split')
    print('')
    print('=== VERDICT ===')
    for m in ('corrected', 'raw'):
        both = R1[m]['K1_repeat'] == 'LIVES' and R2[m]['K1_repeat'] == 'LIVES'
        print(f'{m}: big jumps are real crossings (repeat in independent halves, both splits): {both}')
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, OUT + '_numbers.json'), 'w') as f:
        json.dump({'odd_even': R1, 'shuffled': R2}, f, indent=1)
    if os.path.isdir('/content/drive/MyDrive'):
        import shutil
        dest = os.path.join('/content/drive/MyDrive', OUT); os.makedirs(dest, exist_ok=True)
        shutil.copy(os.path.join(OUT, OUT + '_numbers.json'), dest)
        print('numbers saved to your Drive folder', OUT)
    return R1, R2

if __name__ == '__main__':
    main()
