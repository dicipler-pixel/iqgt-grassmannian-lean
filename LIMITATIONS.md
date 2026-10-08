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
* Section 5 of the rebuilt paper: only the algebra is formalized (the SLD solution `L = 2V`, the
  trace step `2 Tr(ΠVV) = Tr(VV)`, the finite-sum band and the finite-sum Cauchy–Schwarz length
  step). The friction law itself, its slow-driving derivation, the passage from sums to integrals,
  and the minimum-driving-time bound are not formalized; they are checked numerically.
* The redistribution flow `Π̇ = [[Ω,Π],Π]` (Brockett's double-bracket flow) is not formalized; only
  the operator identities are. Neither is the pullback to parameter space (Sec. 4).
* `Core.lean` is checked against a minimal ring declared in the file, with no library.
