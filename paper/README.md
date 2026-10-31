# Paper: Intrinsic Quantum Geometric Tensor on the Grassmannian and Spectral Lower Bounds on Dissipation

Jeromie Beasley

The paper (`iqgt_grassmannian_v7_3.html`, `iqgt_grassmannian_v7_3.pdf`), the Grassmannian Mixer
(`grassmannian_mixer_v8.html`, open it in any browser), and the code, checks and data behind them.
The Lean proofs are at the top of this repository.

## What is in this archive

| Folder | Contents |
| :--- | :--- |
| `checks/` | One numerical check script per section of the paper, each with its pass/fail conditions written at the top before it was run, and the log of its run. |
| `vault_tests/` | Tests of older ideas that were tried and recorded rather than discarded, including those that failed. `test_vault.html` lists every one with its result. |
| `figures/` | The figure scripts and the 22 figures as SVG. |
| `mixer/` | Source of the Grassmannian Mixer v3: maths core `core_v2.js`, its test `core_v2_test.js`, the bench script `app_v2.js`, the animated stage `stage.js` (seven eyes from the Tip and Tail Mixer), the page shell `shell_v4.html` and the condensate table `bec.json`. |
| `data/` | The condensate numbers used in Sections 10 and 12 (`bec_numbers.json`, `bec_cuts_by_split.json`) and the measured eigen-images 2 and 3 shown in the mixer (`bec_eigenimages_2_3.json`, 104 × 84; overlap 0.012, stripes offset 90.6°). |
| `build/` | The section sources of the paper and the scripts that assemble the HTML and print the PDF. |
| `test_vault.html` | Survived, failed-and-kept, and open tests. |

## Running the checks

Everything the checks and the figures need is in `../requirements.txt`, used with Python 3.12
(tested with 3.12.14):

```
python3 -m venv .venv
. .venv/bin/activate
pip install -r ../requirements.txt
cd checks
python3 s2_checks_v2.py
```

The first three entries of `requirements.txt` (NumPy, SciPy, Matplotlib) cover every script but
one. The last three (`xarray`, `netCDF4`, `scikit-learn`) are used only by
`qhe_crossing_split_v1.py`, which also needs the heat-engine data described below; without that
data it stops with the message `could not find the two .nc files …` rather than an import error.

The Lean audit (`../scripts/verify.py`) imports nothing outside the Python standard library, so it
runs with a bare Python 3.

One script writes into the folder it is run from: `ring_tiptail_check.py` saves `ring_rows.npy`,
the ring rows that `figures/fig_v7_3.py` reads. That file is a cache of the run, not an input, and
re-running the check will show it modified.

Expected results:

| Script | Paper | Result |
| :--- | :--- | :--- |
| `s2_checks_v2.py` | §2 framework | 23 / 23 |
| `s3_checks.py` | §3 spectral realization, gap bound | 11 / 11 |
| `prop35_check.py` | Proposition 3.5, Hilbert–Schmidt step (400 random cases, adjacent and non-adjacent spectral sets) | 2 / 2 |
| `s4_checks_v2.py` | §4 pullback, redistribution flow | 11 / 11 |
| `s5s6_checks_v2.py` | §5 dissipation, §6 curvature bound | 19 / 19 |
| `s7_checks_v2.py` | §7 transport tensors | 9 / 9 |
| `s8_checks.py`, `s8_fold_check.py` | §8 geodesic refraction across degenerate strata | 11 / 11, 1 / 1 |
| `s9_checks_v2.py` | §9 mixed states | 11 / 11 |
| `s10_checks.py` | §10 data | 8 / 9: the one FAIL is the tail-versus-core angle, which the noise reaches; it is reported as such in §10 and kept in the vault |
| `s10_swap_check_v2.py` | §10 edge balance | 2 / 2 |
| `s12_checks.py` | §12 corridor and cuts | 4 / 4 |
| `sofic_bridge.py` | §12.4 permutations and principal angles | 6 / 6 |
| `rigidity_checks.py` | §6.3 equality cases (gap bound, speed limit, constant speed, curvature deficit (6.3), trace condition), §8.4–8.5 (metric dipole, Christoffel divergence, silent second order, three-body fold), §12 (sandwich equality, Kuramoto threshold), Example 6.5 (holonomy) | 12 / 12 |
| `spinj_check.py` | Appendix C, Proposition 6.4: the lowest spin-j band (j = 1/2 … 4) meets the curvature bound and the trace condition with equality, Chern number 2j; broken control lifted | 6 / 6 spins + 5 controls |
| `silent_order_check.py` | Proposition 8.6: first nonvanishing order m = 2 … 5 is tangent, distance ε^m/m!, even m folds, odd m passes through | 4 / 4 |
| `eta_dipole_check.py` | §7.5: the regularised dipole weight vanishes at a closed gap, is the derivative of 1/(Δ²+η²), peaks at Δ = η/√3 with height (3√3/8)/η³, and returns 2/Δ³ as η → 0 | PASS |
| `ring_tiptail_check.py` | §12.1, Example 12.2: light both ways round a ring; tail under its Davis–Kahan bound, tip as 1/δ, mixing angle 22.5° at Ω = 2κ, Adler lock-in threshold | PASS |
| `app_checks_v2.py` | Appendices A–D | 13 / 13 |
| `spin1_kahler_check.py` | Appendix C, why the spin-1 band meets the bound: coherent-state (Veronese) form, saturation, generic three-band control | 3 / 3 |
| `polariton_check.py` | Appendix B, the polariton measurement: a two-angle extraction forces 4 det g = Bz² on any field | equality to 10⁻¹⁶, as stated |
| `qhe_crossing_split_v1.py` | §10 crossing test | needs the heat-engine data and three extra packages (below) |

`qhe_crossing_split_v1.py` reads the single-shot data of Uusnäkki et al., Zenodo
doi:10.5281/zenodo.20023151. Download the record's data files, install the last three packages of
`../requirements.txt`, and run the script in the same folder as the two `.nc` files.

## The Lean proofs

```
cd ..
lake exe cache get
lake build
python3 scripts/verify.py
```

The top-level `README.md` maps every theorem to the paper, and `LIMITATIONS.md` says exactly what is
and is not formalized. The proofs are also checked on every commit by the repository's
GitHub Actions.

## The mixer

`node mixer/core_v2_test.js` runs the twelve checks of the mixer's maths core against the paper's
numbers: turning points 0.32981248 (polar stratum, 35°) and −1.0969 (gap-narrowing stratum,
50°), the swap and cycle lengths 3.512 and 2.896, speed under the gap ceiling, curvature inside
the envelope, the tail holding while the tip turns, two-level and spin-1 bands on the curvature bound,
generic three-level bands strictly inside it, and the Fisher speed from fidelity equal to the metric
speed. The published mixer page is these files inlined into `shell_v4.html`.

## Rebuilding the paper

```
cd build
npm i mathjax-full@3.2.2
node assemble_rebuild.js iqgt_paper.html s1_body_v7.html s2_body_v1.html s3_body_v4.html s4_body_v2.html s5_body_v4.html s6_body_v4.html s7_body_v3.html s8_body_v5.html s9_body_v2.html s10_body_v7.html s11_body_v6.html s12_body_v8.html s13_body_v2.html app_body_v6.html refs_body_v7.html
```

This reproduces the published HTML byte for byte (the cover image `cover_v7.jpg` is read from the same folder). For the PDF, install
`@fontsource/newsreader@5.1.0 @fontsource/spectral@5.1.0 @fontsource/jetbrains-mono@5.1.0 @fontsource/archivo-black@5.1.0` and
Playwright, then run `node make_pdf.js iqgt_paper.html iqgt_paper.pdf`.

## Data sources

- Quantum heat engine: T. Uusnäkki et al., arXiv:2502.20143; data and analysis codes,
  doi:10.5281/zenodo.20023151.
- Exciton-polariton quantum geometric tensor (Appendix B, cited, not reanalysed): A. Gianfrate et al.,
  Nature 578, 381 (2020), arXiv:1901.03219.
- Dark-soliton Bose–Einstein condensate: A. Fritsch et al., Machine Learning: Science and
  Technology (2022), arXiv:2205.09114; data NIST mds2-2363.

## Licence

Code MIT, text and figures CC BY 4.0, as in `../LICENSING.md`.
