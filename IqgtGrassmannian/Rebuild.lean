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
end IqgtGrassmannian.Rebuild
