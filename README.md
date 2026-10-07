<div align="center">

# Intrinsic Quantum Geometric Tensor on the Grassmannian and Spectral Lower Bounds on Dissipation — Lean proofs

[![Lean proof check](https://github.com/dicipler-pixel/iqgt-grassmannian-lean/actions/workflows/build.yml/badge.svg)](https://github.com/dicipler-pixel/iqgt-grassmannian-lean/actions/workflows/build.yml)
![Lean](https://img.shields.io/badge/Lean-v4.34.1-blue)
![Theorems](https://img.shields.io/badge/theorems-11-2EA043)
![sorry](https://img.shields.io/badge/sorry-0-2EA043)
![Code: MIT](https://img.shields.io/badge/code-MIT-lightgrey)
![Text: CC BY 4.0](https://img.shields.io/badge/text-CC%20BY%204.0-lightgrey)
[![Paper DOI](https://img.shields.io/badge/paper-10.5281%2Fzenodo.20768258-blue)](https://doi.org/10.5281/zenodo.20768258)

Jeromie Beasley

</div>

---

## The idea in one line

The quantum geometric tensor `Q_Π(V, W) = Tr(ΠVW)` is a Frobenius inner product in disguise. Its
Gram matrix is therefore positive semidefinite, and that alone bounds the Berry curvature by the
quantum metric: `|Ω| ≤ 2√(g(V,V)g(W,W) − g(V,W)²)`. The inequality is checked here in Lean for
every finite dimension.

## What is proved

| Paper | Result | Theorem |
| :--- | :--- | :--- |
| Corollary 3.1 | A tangent vector at a projector (`V = ΠV + VΠ`) is purely off-diagonal: `ΠVΠ = 0`, `(1−Π)V(1−Π) = 0`, `V = ΠV(1−Π) + (1−Π)VΠ` | `tangent_offdiag` |
| Theorem 3.2 | `Tr V² = 2 Tr(XX†)` for `V = [[0, X], [X†, 0]]` | `trace_sq_offdiag` |
| Sec. 2 | `Q(V,W) = Tr((VΠ)†(WΠ))`; the pairing is Hermitian and equals a sum of entries | `Q_eq`, `frob_herm`, `frob_sum` |
| Theorem 6.1 | The Gram matrix of `Q` is positive semidefinite, so `(Im Q(V,W))² ≤ g(V,V)g(W,W) − g(V,W)²`, i.e. `\|Ω\| ≤ 2√(det g)` | `gram_entries`, `gram_det_nonneg`, `metric_control` |
| Corollary 6.1 | Metric rigidity `g(V,V) = 0` forces `Ω(V,W) = 0` | `rigidity_kills_curvature` |
| Sec. 6, two-level case | Lagrange's identity `\|A\|²\|B\|² − (A·B)² = \|A×B\|²`, the saturation step | `lagrange_identity` |
| Lemma 3.1 | The partial fraction whose residue gives the inverse-gap weight | `partial_fraction` |

Theorem 6.1 is bridge B1 of the Operator-First Atlas ("curvature cannot exceed what the metric
affords"). The file is [`IqgtGrassmannian/Basic.lean`](IqgtGrassmannian/Basic.lean). What is not
proved is in [`LIMITATIONS.md`](LIMITATIONS.md).

## How it is checked

Every push runs [the proof check](.github/workflows/build.yml): build against Lean v4.34.1 and
Mathlib v4.34.1, independent replay in Lean's kernel checker, an axiom audit (only `propext`,
`Classical.choice`, `Quot.sound`), and three deliberately false statements that must fail.

## The paper

*Intrinsic Quantum Geometric Tensor on the Grassmannian and Spectral Lower Bounds on Dissipation*, Jeromie Beasley. DOI
[10.5281/zenodo.20768258](https://doi.org/10.5281/zenodo.20768258) (always opens the newest version).

## Licence

Copyright (c) 2026 Jeromie Beasley. Code and proofs: [MIT](LICENSE). Written text:
[CC BY 4.0](LICENSE-CC-BY-4.0.md). See [`LICENSING.md`](LICENSING.md). Citation metadata is in
[`CITATION.cff`](CITATION.cff); how AI tools were used is stated in [`AI_USE.md`](AI_USE.md).
