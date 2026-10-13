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
* Proposition 3.5: only the scalar Hilbert–Schmidt step (a finite sum of squared entries divided by
  gaps of modulus at least `Δ`) is formalized; the identification of those entries with the
  resolvent formula (3.2) is not.
* Section 12.4: the sandwich between the Hamming count and the intrinsic length and the pairing of
  the cycle Gram eigenvalues are formalized as statements about real angles. The trace identity
  `Tr((P−Q)²) = Tr P + Tr Q − 2 Tr(PQ)` is formalized (`hs_dist_proj`); reading `Tr(PQ)` as
  `Σ cos²θ` over principal angles, and hence `d_H = (2/n) Σ sin²θ`, and the geodesic length
  `√Σθ²`, are checked numerically (`paper/checks/sofic_bridge.py`), not in Lean.
* Theorem 6.1 in block form (`curvature_bound`) is proved for vectors in `EuclideanSpace ℂ ι`; the
  identification of a tangent vector with its cross-gap block is stated in words, not formalized.
* Section 6.3 (equality cases, paper v7): only the scalar and algebraic cores are formalized: the
  termwise equality of the Hilbert–Schmidt step (`gap_bound_hs_step_eq_iff`), the constant-speed
  equality of the finite-sum length form (`constant_speed_iff`), equality in Cauchy–Schwarz over `ℂ`
  (`saturation_iff_complex_line`, a restatement of Mathlib's `norm_inner_eq_norm_iff`) and the
  two-number arithmetic–geometric mean step of the trace condition (`am_gm_eq_iff`). The singular-value
  step of the gap bound, the equality case of the minimum driving time, the deficit identity (6.3)
  as an identity of tangent vectors, and the pullback of the trace condition are checked numerically
  (`paper/checks/rigidity_checks.py`), not in Lean.
* Section 8.5: Lemma 8.5 is formalized as "a differentiable curve with zero derivative is constant"
  and Proposition 8.6 as a ring identity (`silent_second_order`), taking the twice-differentiated
  idempotence `Π″ = Π″Π + 2Π′Π′ + ΠΠ″` as a hypothesis; the differentiation itself, the Taylor
  estimate (8.3), the Christoffel divergence and the three-body fold are checked numerically.
* Section 12.4: the sandwich's equality case is formalized for real angles and for unit eigenvalues
  (`sandwich_upper_eq_iff`, `half_phase_right_iff`, `involution_eigenvalue`); the passage from a
  permutation matrix to its principal angles is not.
