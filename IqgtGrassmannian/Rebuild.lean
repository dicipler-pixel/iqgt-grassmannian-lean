/-
Intrinsic Quantum Geometric Tensor on the Grassmannian and Spectral Lower Bounds on Dissipation
(Jeromie Beasley): the statements added by the section-by-section rebuild (Sections 2, 4, 5),
in Mathlib's language. Numbering follows the rebuilt paper. `Basic.lean` is unchanged.

* §2 (converse of (2.1)) both off-diagonal blocks vanishing makes `V` tangent; one block is not enough.
* §2 the complement: `-V` is tangent at `1 - Π`, so the metric of `Π` and of `1 - Π` agree.
* §2 Remark 2.3: `Tr(ΠVW) + Tr(ΠWV) = Tr(VW)` as soon as `V` is tangent.
* §4.4 the redistribution operator `F = (1-Π)ΩΠ + Π Ω(1-Π) = [[Ω,Π],Π]`: it is tangent, and
  `F = 0 ↔ ΩΠ = ΠΩ`. Stated for an idempotent in any ring.
* §5 Lemma 5.1: `L = 2V` solves the SLD equation for `ρ = Π/k`, and `2 Tr(ΠVV) = Tr(VV)`,
  so `F_Q = 4 g / k`.
* §5 Theorem 5.3: the zero-temperature band (a convex reweighting of the metric) and the
  length form `(Σ √gᵢ)² ≤ N Σ gᵢ` (Cauchy–Schwarz).
* §3 Proposition 3.5: the Hilbert–Schmidt step of the gap bound.
* §12.4 the sofic bridge: the chordal distance `Tr((P-Q)²) = Tr P + Tr Q - 2 Tr(PQ)`, the sandwich `Σ sin²θ ≤ Σ θ² ≤ (π/2)² Σ sin²θ` and the paired
  cycle Gram eigenvalues.
* §6 Theorem 6.1 in block form: `Ω² ≤ 4 (g_VV g_WW − g_VW²)` for the cross-gap blocks, from
  Cauchy–Schwarz in `EuclideanSpace ℂ ι`.
-/
import Mathlib

namespace IqgtGrassmannian.Rebuild

open Matrix BigOperators

section Ring
variable {R : Type*} [Ring R]

/-- The linearised projector constraint `V = ΠV + VΠ` (§2). -/
def IsTangent (P V : R) : Prop := V = P * V + V * P

/-- §2, converse of (2.1): if both diagonal blocks vanish, `V` is tangent. -/
theorem tangent_of_blocks {P V : R} (h1 : P * V * P = 0) (h2 : (1 - P) * V * (1 - P) = 0) :
    IsTangent P V := by
  have e : (1 - P) * V * (1 - P) = V - (P * V + V * P) + P * V * P := by noncomm_ring
  rw [e, h1, add_zero, sub_eq_zero] at h2
  exact h2

/-- §2: `-V` moves the complement `1 - Π`. -/
theorem tangent_compl {P V : R} (hV : IsTangent P V) : IsTangent (1 - P) (-V) := by
  unfold IsTangent at *
  have e : (1 - P) * -V + -V * (1 - P) = -V - V + (P * V + V * P) := by noncomm_ring
  rw [e, ← hV]; abel

/-- The redistribution operator (§4.4). -/
def redist (Ω P : R) : R := (1 - P) * Ω * P + P * Ω * (1 - P)

theorem redist_double_bracket (Ω P : R) (hP : IsIdempotentElem P) :
    redist Ω P = (Ω * P - P * Ω) * P - P * (Ω * P - P * Ω) := by
  have h := hP.eq
  calc redist Ω P = Ω * P + P * Ω - P * Ω * P - P * Ω * P := by unfold redist; noncomm_ring
    _ = Ω * (P * P) + (P * P) * Ω - P * Ω * P - P * Ω * P := by rw [h]
    _ = (Ω * P - P * Ω) * P - P * (Ω * P - P * Ω) := by noncomm_ring

theorem proj_mul_compl {P : R} (hP : IsIdempotentElem P) : P * (1 - P) = 0 := by
  rw [mul_sub, mul_one, hP.eq, sub_self]

theorem compl_mul_proj {P : R} (hP : IsIdempotentElem P) : (1 - P) * P = 0 := by
  rw [sub_mul, one_mul, hP.eq, sub_self]

theorem proj_mul_redist (Ω P : R) (hP : IsIdempotentElem P) :
    P * redist Ω P = P * Ω * (1 - P) := by
  have a : P * redist Ω P = (P * (1 - P)) * Ω * P + (P * P) * Ω * (1 - P) := by
    unfold redist; noncomm_ring
  rw [a, proj_mul_compl hP, hP.eq, zero_mul, zero_mul, zero_add]

theorem redist_mul_proj (Ω P : R) (hP : IsIdempotentElem P) :
    redist Ω P * P = (1 - P) * Ω * P := by
  have b : redist Ω P * P = (1 - P) * Ω * (P * P) + P * Ω * ((1 - P) * P) := by
    unfold redist; noncomm_ring
  rw [b, compl_mul_proj hP, hP.eq, mul_zero, add_zero]

/-- §4.4: the redistribution operator is a tangent vector at `Π`. -/
theorem redist_tangent (Ω P : R) (hP : IsIdempotentElem P) : IsTangent P (redist Ω P) := by
  unfold IsTangent
  rw [proj_mul_redist Ω P hP, redist_mul_proj Ω P hP]
  unfold redist; exact add_comm _ _

/-- §4.4: `F(Ω, Π) = 0` exactly when `Ω` commutes with `Π`. -/
theorem redist_eq_zero_iff (Ω P : R) (hP : IsIdempotentElem P) :
    redist Ω P = 0 ↔ Ω * P = P * Ω := by
  constructor
  · intro hF
    have a := proj_mul_redist Ω P hP
    have b := redist_mul_proj Ω P hP
    rw [hF, mul_zero, mul_sub, mul_one] at a
    rw [hF, zero_mul, sub_mul, one_mul, sub_mul] at b
    have e1 : P * Ω = P * Ω * P := sub_eq_zero.mp a.symm
    have e2 : Ω * P = P * Ω * P := sub_eq_zero.mp b.symm
    rw [e2, ← e1]
  · intro hc
    rw [redist_double_bracket Ω P hP, hc, sub_self, zero_mul, mul_zero, sub_self]

/-- Lemma 5.1, algebraic core: `L = 2V` solves `Lρ + ρL = 2ρ'` for `ρ = Π/k` (the `1/k` cancels). -/
theorem sld_two_tangent {P V : R} (hV : IsTangent P V) : (2 * V) * P + P * (2 * V) = 2 * V := by
  have e : (2 * V) * P + P * (2 * V) = 2 * (P * V + V * P) := by noncomm_ring
  rw [e, ← hV]

end Ring

section Trace
variable {n 𝕜 : Type*} [Fintype n] [DecidableEq n] [CommRing 𝕜]

/-- Remark 2.3: `Tr(ΠVW) + Tr(ΠWV) = Tr(VW)` as soon as `V` is tangent at `Π`. -/
theorem trace_remark_2_3 {P V W : Matrix n n 𝕜} (hV : IsTangent P V) :
    trace (P * V * W) + trace (P * W * V) = trace (V * W) := by
  have h : V * W = (P * V + V * P) * W := by rw [← hV]
  rw [h, add_mul, trace_add, trace_mul_cycle V P W, trace_mul_cycle W V P]

/-- Lemma 5.1, trace step: `2 Tr(ΠVV) = Tr(VV)` for tangent `V`, hence `F_Q = 4 g / k`. -/
theorem two_trace_PVV {P V : Matrix n n 𝕜} (hV : IsTangent P V) :
    2 * trace (P * V * V) = trace (V * V) := by
  rw [two_mul]; exact trace_remark_2_3 (W := V) hV

/-- §12.4, the chordal distance between two projectors: for idempotent matrices,
`Tr((P - Q)²) = Tr P + Tr Q - 2 Tr(PQ)`. For orthogonal projectors of equal rank `k` this is
`2k - 2 Tr(PQ) = 2 Σ sin²θⱼ`, the identity behind the Hamming count of a permutation. -/
theorem hs_dist_proj {P Q : Matrix n n 𝕜} (hP : P * P = P) (hQ : Q * Q = Q) :
    trace ((P - Q) * (P - Q)) = trace P + trace Q - 2 * trace (P * Q) := by
  have e : (P - Q) * (P - Q) = P * P - P * Q - Q * P + Q * Q := by noncomm_ring
  rw [e, hP, hQ]
  first
  | (simp only [trace_add, trace_sub, trace_mul_comm Q P]; ring)
  | (simp [trace_add, trace_sub, trace_mul_comm Q P]; ring)

end Trace

section Friction
variable {ι : Type*}

/-- Theorem 5.3(ii), the zero-temperature band: a reweighting of the metric by weights in `[lo, hi]`. -/
theorem friction_band (s : Finset ι) (w g : ι → ℝ) (hg : ∀ i ∈ s, 0 ≤ g i)
    (lo hi : ℝ) (hlo : ∀ i ∈ s, lo ≤ w i) (hhi : ∀ i ∈ s, w i ≤ hi) :
    lo * ∑ i ∈ s, g i ≤ ∑ i ∈ s, w i * g i ∧ ∑ i ∈ s, w i * g i ≤ hi * ∑ i ∈ s, g i := by
  constructor
  · rw [Finset.mul_sum]
    exact Finset.sum_le_sum fun i hi' => mul_le_mul_of_nonneg_right (hlo i hi') (hg i hi')
  · rw [Finset.mul_sum]
    exact Finset.sum_le_sum fun i hi' => mul_le_mul_of_nonneg_right (hhi i hi') (hg i hi')

/-- Theorem 5.3, length form (discretised): `(Σ √gᵢ)² ≤ N · Σ gᵢ`, the step `ℒ² ≤ T ∫ g`. -/
theorem length_sq_le (s : Finset ι) (g : ι → ℝ) (hg : ∀ i ∈ s, 0 ≤ g i) :
    (∑ i ∈ s, Real.sqrt (g i)) ^ 2 ≤ (s.card : ℝ) * ∑ i ∈ s, g i := by
  have h := Finset.sum_mul_sq_le_sq_mul_sq s (fun _ => (1 : ℝ)) (fun i => Real.sqrt (g i))
  simp only [one_mul, one_pow, Finset.sum_const, nsmul_eq_mul, mul_one] at h
  calc (∑ i ∈ s, Real.sqrt (g i)) ^ 2 ≤ (s.card : ℝ) * ∑ i ∈ s, Real.sqrt (g i) ^ 2 := h
    _ = (s.card : ℝ) * ∑ i ∈ s, g i := by
        congr 1; exact Finset.sum_congr rfl fun i hi => Real.sq_sqrt (hg i hi)

end Friction
section Transport

/-- §7.2: friction read from the response function, one pair at a time:
`2Δg·τ/(1+τ²Δ²) = τ·[2Δg − 2Δ³g/(Δ² + 1/τ²)]`, i.e. `ζ_rot = τ[α(0) − α(i/τ)]`. -/
theorem friction_response_pair (Δ g τ : ℝ) (hτ : 0 < τ) :
    2 * Δ * g * τ / (1 + τ ^ 2 * Δ ^ 2) = τ * (2 * Δ * g - 2 * Δ ^ 3 * g / (Δ ^ 2 + 1 / τ ^ 2)) := by
  have h1 : (0 : ℝ) < 1 + τ ^ 2 * Δ ^ 2 := by positivity
  have h2 : (0 : ℝ) < Δ ^ 2 + 1 / τ ^ 2 := by positivity
  field_simp
  ring

/-- §7.5: the regularised pair metric is bounded at a closed gap: `a/(Δ² + η²) ≤ a/η²` for `a ≥ 0`. -/
theorem eta_regularised_bound (a Δ η : ℝ) (ha : 0 ≤ a) (hη : 0 < η) :
    a / (Δ ^ 2 + η ^ 2) ≤ a / η ^ 2 := by
  apply div_le_div_of_nonneg_left ha (by positivity)
  nlinarith [sq_nonneg Δ]

/-- §8.3, total internal reflection: if `C = g∥ sin²θ` with `sin²θ ≤ 1` and `g∥ ≥ 0`, then the
trajectory only visits `g∥ ≥ C`; at a stratum `g∥ = 0` it would need `C = 0`. -/
theorem reflection_bound (gpar s2 : ℝ) (hg : 0 ≤ gpar) (hs : 0 ≤ s2) (hs1 : s2 ≤ 1) :
    gpar * s2 ≤ gpar := by
  nlinarith

theorem stratum_needs_zero (s2 : ℝ) : (0 : ℝ) * s2 = 0 := zero_mul s2

end Transport

section Mixed
variable {R : Type*} [Ring R]

/-- §9, Remark 9.2: anything of the form `Lρ + ρL` with `ρ` supported on `Ran Π` has no
kernel–kernel block, so a variation with such a block can never solve the SLD equation. -/
theorem sld_kernel_block {P L : R} (hP : IsIdempotentElem P) :
    (1 - P) * (L * P + P * L) * (1 - P) = 0 := by
  have h1 : P * (1 - P) = 0 := proj_mul_compl hP
  have h2 : (1 - P) * P = 0 := compl_mul_proj hP
  have e : (1 - P) * (L * P + P * L) * (1 - P)
      = (1 - P) * L * (P * (1 - P)) + ((1 - P) * P) * L * (1 - P) := by noncomm_ring
  rw [e, h1, h2]; simp

end Mixed

section GapLaw
variable {R : Type*} [Ring R]

/-- §10.2, the gap identity behind the tip/tail law (finite Davis–Kahan). If `P` is a spectral
projector of `A` and `Q` the complementary spectral projector of `A' = A + E` (each commuting with
its operator), then `(QA'Q)(QP) − (QP)(PAP) = QEP`. The block `QP` measures how far the subspace
moved; this Sylvester equation is what turns a gap into a bound on it. -/
theorem gap_identity {P Q A A' : R} (hP : IsIdempotentElem P) (hQ : IsIdempotentElem Q)
    (hA : A * P = P * A) (hA' : A' * Q = Q * A') :
    (Q * A' * Q) * (Q * P) - (Q * P) * (P * A * P) = Q * (A' - A) * P := by
  have hPP : P * P = P := hP.eq
  have hQQ : Q * Q = Q := hQ.eq
  have e1 : Q * A' * Q * (Q * P) = Q * A' * P := by
    calc Q * A' * Q * (Q * P) = Q * A' * (Q * Q) * P := by noncomm_ring
      _ = Q * A' * Q * P := by rw [hQQ]
      _ = Q * (A' * Q) * P := by noncomm_ring
      _ = Q * (Q * A') * P := by rw [hA']
      _ = Q * Q * A' * P := by noncomm_ring
      _ = Q * A' * P := by rw [hQQ]
  have e2 : Q * P * (P * A * P) = Q * A * P := by
    calc Q * P * (P * A * P) = Q * (P * P) * A * P := by noncomm_ring
      _ = Q * P * A * P := by rw [hPP]
      _ = Q * (P * A) * P := by noncomm_ring
      _ = Q * (A * P) * P := by rw [hA]
      _ = Q * A * (P * P) := by noncomm_ring
      _ = Q * A * P := by rw [hPP]
  rw [e1, e2]; noncomm_ring

end GapLaw

/-- §10.2, the tip is free at a tie: inside an exactly degenerate eigenspace every combination is
again an eigenvector, so no single vector is selected. Only the subspace is determined. -/
theorem tip_free_at_tie {K V : Type*} [CommRing K] [AddCommGroup V] [Module K V]
    (A : V →ₗ[K] V) (c a b : K) {v w : V} (hv : A v = c • v) (hw : A w = c • w) :
    A (a • v + b • w) = c • (a • v + b • w) := by
  first
  | (rw [map_add, map_smul, map_smul, hv, hw, smul_add, smul_comm a c v, smul_comm b c w])
  | (simp [hv, hw, smul_smul, mul_comm])

/-- Appendix A, the second resolvent identity with its sign: for invertible `a = z - Σ` and
`b = z - Σ - δΣ`, `a⁻¹ - b⁻¹ = a⁻¹ (b - a) b⁻¹`, and `b - a = -δΣ`. Holds in any ring. -/
theorem second_resolvent {R : Type*} [Ring R] (a b : Rˣ) :
    ((a⁻¹ : Rˣ) : R) - ((b⁻¹ : Rˣ) : R) = ((a⁻¹ : Rˣ) : R) * ((b : R) - (a : R)) * ((b⁻¹ : Rˣ) : R) := by
  rw [mul_sub, sub_mul, mul_assoc, Units.mul_inv, mul_one, Units.inv_mul, one_mul]

section Bounds
variable {ι : Type*}

/-- §3, Proposition 3.5, the Hilbert–Schmidt step of the gap bound: dividing each entry of
`Π δΣ (1-Π)` by a gap of modulus at least `Δ` divides the sum of squared entries by at least `Δ²`.
This is the step that bounds the Hilbert–Schmidt norm (dividing entries does not, in general,
bound the operator norm). -/
theorem gap_bound_hs_step (s : Finset ι) (a d : ι → ℝ) (Δ : ℝ) (hΔ : 0 < Δ)
    (hd : ∀ i ∈ s, Δ ≤ |d i|) :
    ∑ i ∈ s, a i ^ 2 / d i ^ 2 ≤ (∑ i ∈ s, a i ^ 2) / Δ ^ 2 := by
  rw [Finset.sum_div]
  apply Finset.sum_le_sum
  intro i hi
  have hdi := hd i hi
  have hm := mul_le_mul hdi hdi hΔ.le (abs_nonneg (d i))
  have h1 : Δ ^ 2 ≤ d i ^ 2 := by nlinarith [sq_abs (d i)]
  have h2 : 0 < Δ ^ 2 := by positivity
  first
  | exact div_le_div_of_nonneg_left (sq_nonneg (a i)) h2 h1
  | (rw [div_eq_mul_inv, div_eq_mul_inv]
     exact mul_le_mul_of_nonneg_left ((inv_le_inv₀ (lt_of_lt_of_le h2 h1) h2).mpr h1) (sq_nonneg (a i)))

/-- §12.4, lower half of the sandwich, one principal angle: `sin² θ ≤ θ²` on `[0, π/2]`. -/
theorem sin_sq_le_angle_sq (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ Real.pi / 2) :
    Real.sin x ^ 2 ≤ x ^ 2 := by
  have hs : Real.sin x ≤ x := Real.sin_le h0
  have hpos : 0 ≤ Real.sin x :=
    Real.sin_nonneg_of_nonneg_of_le_pi h0 (by linarith [Real.pi_pos])
  nlinarith [mul_le_mul hs hs hpos h0]

/-- §12.4, upper half of the sandwich, one principal angle (Jordan's inequality):
`θ² ≤ (π/2)² sin² θ` on `[0, π/2]`. -/
theorem angle_sq_le_jordan (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ Real.pi / 2) :
    x ^ 2 ≤ (Real.pi / 2) ^ 2 * Real.sin x ^ 2 := by
  have hpi : 0 < Real.pi := Real.pi_pos
  have hj : 2 / Real.pi * x ≤ Real.sin x := by
    first
    | exact Real.mul_le_sin h0 h1
    | exact Real.two_div_pi_mul_le_sin h0 h1
  have hpi' : Real.pi ≠ 0 := hpi.ne'
  have he : Real.pi / 2 * (2 / Real.pi * x) = x := by
    first
    | (field_simp; done)
    | (field_simp; ring)
    | (rw [← mul_assoc, div_mul_div_comm, mul_comm Real.pi 2, div_self (by positivity), one_mul])
  have hx : x ≤ Real.pi / 2 * Real.sin x := by
    calc x = Real.pi / 2 * (2 / Real.pi * x) := he.symm
      _ ≤ Real.pi / 2 * Real.sin x := mul_le_mul_of_nonneg_left hj (by positivity)
  calc x ^ 2 = x * x := sq x
    _ ≤ (Real.pi / 2 * Real.sin x) * (Real.pi / 2 * Real.sin x) :=
        mul_le_mul hx hx h0 (le_trans h0 hx)
    _ = (Real.pi / 2) ^ 2 * Real.sin x ^ 2 := by ring

/-- §12.4, the sandwich between the Hamming count and the intrinsic length: over principal
angles in `[0, π/2]`, `Σ sin² θ ≤ Σ θ² ≤ (π/2)² Σ sin² θ`. With `d_H = (2/n) Σ sin² θ` this is
`d_H / 2 ≤ (1/n) Σ θ² ≤ (π²/8) d_H`. -/
theorem sofic_sandwich (s : Finset ι) (θ : ι → ℝ) (h0 : ∀ i ∈ s, 0 ≤ θ i)
    (h1 : ∀ i ∈ s, θ i ≤ Real.pi / 2) :
    ∑ i ∈ s, Real.sin (θ i) ^ 2 ≤ ∑ i ∈ s, θ i ^ 2 ∧
    ∑ i ∈ s, θ i ^ 2 ≤ (Real.pi / 2) ^ 2 * ∑ i ∈ s, Real.sin (θ i) ^ 2 := by
  refine ⟨Finset.sum_le_sum fun i hi => sin_sq_le_angle_sq _ (h0 i hi) (h1 i hi), ?_⟩
  rw [Finset.mul_sum]
  exact Finset.sum_le_sum fun i hi => angle_sq_le_jordan _ (h0 i hi) (h1 i hi)

/-- §12.4, the displacement Gram eigenvalues of an `ℓ`-cycle come in pairs:
`2 sin²(π(ℓ-j)/ℓ) = 2 sin²(πj/ℓ)`, the sine and cosine of one cyclic motion. -/
theorem cycle_gram_pair (l j : ℝ) (hl : l ≠ 0) :
    2 * Real.sin (Real.pi * (l - j) / l) ^ 2 = 2 * Real.sin (Real.pi * j / l) ^ 2 := by
  have h : Real.pi * (l - j) / l = Real.pi - Real.pi * j / l := by
    first
    | (field_simp; done)
    | (field_simp; ring)
    | (rw [mul_sub, sub_div, mul_div_assoc, div_self hl, mul_one])
  rw [h, Real.sin_pi_sub]

end Bounds

section Curvature
variable {ι : Type*} [Fintype ι]

/-- Theorem 6.1 in block form. For the cross-gap blocks `x, y` of two tangent vectors,
`Q(V,W) = ⟪y, x⟫`, `g_VV = ‖x‖²`, `g_WW = ‖y‖²`, `g_VW = Re Q` and `Ω = −2 Im Q`. Then
`Ω² ≤ 4 (g_VV g_WW − g_VW²)`, which is Cauchy–Schwarz for the blocks. -/
theorem curvature_bound (x y : EuclideanSpace ℂ ι) :
    (2 * (inner ℂ y x).im) ^ 2 ≤ 4 * (‖x‖ ^ 2 * ‖y‖ ^ 2 - (inner ℂ y x).re ^ 2) := by
  have h1 : ‖inner ℂ y x‖ ≤ ‖y‖ * ‖x‖ := norm_inner_le_norm y x
  have h2 : ‖inner ℂ y x‖ ^ 2 = (inner ℂ y x).re ^ 2 + (inner ℂ y x).im ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_apply]; ring
  have h3 : ‖inner ℂ y x‖ ^ 2 ≤ (‖y‖ * ‖x‖) ^ 2 :=
    pow_le_pow_left₀ (norm_nonneg _) h1 2
  nlinarith [h2, h3]

end Curvature

end IqgtGrassmannian.Rebuild
