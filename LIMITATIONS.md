# What is not proved here

Lean proves exactly the statements written, under exactly the hypotheses written.

* The results are proved for matrices, meaning finite-dimensional Hilbert spaces. The separable
  infinite-dimensional setting and the Riesz contour construction are not formalized. Lemma 3.1
  is proved only as the partial-fraction identity used in its residue step; the residue itself
  is not computed in Lean.
* The two-level saturation is proved only as Lagrange's identity in `ℝ³`. Writing `Q` of a
  two-band operator in terms of Bloch vectors, which turns the bound into an equality, is not
  formalized.
* `g` is taken as `Re Q` and `Ω` as `2 Im Q`. That `Q(V,V) = ½ Tr V²` for off-diagonal `V`
  follows from `tangent_offdiag` and `trace_sq_offdiag`, but it is not assembled as one
  statement.
* Theorem 5.1 (the thermodynamic-length bound on entropy production) rests on a cited
  quantum-Fisher-information bound and is not formalized. Neither is the pullback to parameter
  space (Sec. 4).
