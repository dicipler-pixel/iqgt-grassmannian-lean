/-
Intrinsic Quantum Geometric Tensor on the Grassmannian (Jeromie Beasley):
the two equality/defect statements of Section 6, v7.3.1 numbering.

* Proposition 6.4 (first step): for the positive semidefinite Gram matrix `[[a,b],[b,d]]` of two
  tangent deformations, `tr g = a + d ≥ 2√det g = 2√(ad − b²)`; equality holds exactly when
  `b = 0` and `a = d`, i.e. when `g` is conformal to the reference metric.
* Identity (6.3): the deficit of the curvature bound is an exact quadratic form whose kernel is
  the line of `v`; over the reals in which `g` is an inner product, the deficit of `w` equals the
  squared norm of the part of `w` orthogonal to `v`.
-/
import Mathlib

namespace IqgtGrassmannian

open Real
open scoped InnerProductSpace

/-- **Proposition 6.4, first step.** A positive semidefinite 2×2 Gram matrix satisfies
`a + d ≥ 2√(a d − b²)`, which is `tr g ≥ 2√det g`, i.e. `½ tr g ≥ √det g`. -/
theorem trace_ge_two_sqrt_det (a b d : ℝ) (ha : 0 ≤ a) (hd : 0 ≤ d) :
    a + d ≥ 2 * Real.sqrt (a * d - b ^ 2) := by
  have had : 0 ≤ a * d := mul_nonneg ha hd
  have hs : 0 ≤ a + d := by linarith
  have hmono : Real.sqrt (a * d - b ^ 2) ≤ Real.sqrt (a * d) := by
    apply Real.sqrt_le_sqrt
    nlinarith [sq_nonneg b]
  have h4 : 4 * (a * d) ≤ (a + d) ^ 2 := by nlinarith [sq_nonneg (a - d)]
  have hsq : 2 * Real.sqrt (a * d) ≤ a + d := by
    have h1 : Real.sqrt (4 * (a * d)) ≤ Real.sqrt ((a + d) ^ 2) := Real.sqrt_le_sqrt h4
    rw [Real.sqrt_sq hs] at h1
    rw [Real.sqrt_mul (by norm_num) (a * d)] at h1
    norm_num at h1
    linarith
  linarith

/-- **Proposition 6.4, equality case.** Equality in the trace condition forces `g` to be
conformal: `b = 0` and `a = d`. -/
theorem trace_eq_two_sqrt_det_iff (a b d : ℝ) (hpsd : 0 ≤ a * d - b ^ 2) (ha : 0 ≤ a) (hd : 0 ≤ d) :
    a + d = 2 * Real.sqrt (a * d - b ^ 2) ↔ b = 0 ∧ a = d := by
  have had : 0 ≤ a * d := mul_nonneg ha hd
  have hs : 0 ≤ a + d := by linarith
  have hmono : Real.sqrt (a * d - b ^ 2) ≤ Real.sqrt (a * d) := by
    apply Real.sqrt_le_sqrt
    nlinarith [sq_nonneg b]
  have h4 : 4 * (a * d) ≤ (a + d) ^ 2 := by nlinarith [sq_nonneg (a - d)]
  have hsq : 2 * Real.sqrt (a * d) ≤ a + d := by
    have h1 : Real.sqrt (4 * (a * d)) ≤ Real.sqrt ((a + d) ^ 2) := Real.sqrt_le_sqrt h4
    rw [Real.sqrt_sq hs] at h1
    rw [Real.sqrt_mul (by norm_num) (a * d)] at h1
    norm_num at h1
    linarith
  constructor
  · intro h
    -- both steps of the chain are equalities
    have he1 : Real.sqrt (a * d) = Real.sqrt (a * d - b ^ 2) := by linarith
    have he2 : 2 * Real.sqrt (a * d) = a + d := by linarith
    have hb : b ^ 2 = 0 := by
      have hsq : a * d - b ^ 2 = a * d := by
        have h := congrArg (fun x => x ^ 2) he1
        rw [Real.sq_sqrt had, Real.sq_sqrt hpsd] at h
        linarith
      linarith
    have hsq0 : (a - d) ^ 2 = 0 := by
      have h2 : (2 * Real.sqrt (a * d)) ^ 2 = (a + d) ^ 2 := by rw [he2]
      rw [mul_pow, Real.sq_sqrt had] at h2
      linarith
    exact ⟨pow_eq_zero_iff (two_ne_zero) |>.mp hb, sub_eq_zero.mp (pow_eq_zero_iff (two_ne_zero) |>.mp hsq0)⟩
  · rintro ⟨rfl, rfl⟩
    have hsq : a * a - 0 ^ 2 = a ^ 2 := by ring
    rw [hsq, Real.sqrt_sq ha]
    ring

/-- **Identity (6.3), abstract form.** In the real inner product space in which `g` takes its
values, the deficit `g(w,w) − g(v,w)²/g(v,v)` of the curvature bound is the squared norm of the
part of `w` orthogonal to `v`: a nonnegative quadratic form in `w` whose kernel is exactly the
line of `v`. -/
theorem deficit_eq_sq_norm_sub_projection {E : Type*} [NormedAddCommGroup E]
    [InnerProductSpace ℝ E] (v w : E) (hv : v ≠ 0) :
    ‖w‖ ^ 2 - (⟪v, w⟫_ℝ) ^ 2 / ‖v‖ ^ 2
      = ‖w - ((⟪v, w⟫_ℝ) / ‖v‖ ^ 2) • v‖ ^ 2 := by
  have hvv : ‖v‖ ^ 2 ≠ 0 := by
    have : 0 < ‖v‖ := norm_pos_iff.mpr hv
    positivity
  set c : ℝ := (⟪v, w⟫_ℝ) / ‖v‖ ^ 2 with hc
  have hinner : ⟪w, c • v⟫_ℝ = c * ⟪v, w⟫_ℝ := by
    rw [real_inner_smul_right, real_inner_comm v w]
  have hnorm : ‖c • v‖ ^ 2 = c ^ 2 * ‖v‖ ^ 2 := by
    rw [norm_smul, Real.norm_eq_abs, mul_pow, sq_abs]
  have hkey : ‖w - c • v‖ ^ 2 = ‖w‖ ^ 2 - 2 * c * ⟪v, w⟫_ℝ + c ^ 2 * ‖v‖ ^ 2 := by
    rw [norm_sub_sq_real, hinner, hnorm]
    ring
  rw [hkey, hc]
  field_simp
  ring

/-- **Identity (6.3): the kernel.** The deficit vanishes exactly on the line of `v`, so equality in
the curvature bound holds exactly when the two directions span a (complex) line. -/
theorem deficit_eq_zero_iff {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]
    (v w : E) (hv : v ≠ 0) :
    ‖w‖ ^ 2 - (⟪v, w⟫_ℝ) ^ 2 / ‖v‖ ^ 2 = 0 ↔ ∃ c : ℝ, w = c • v := by
  have hvv : ‖v‖ ^ 2 ≠ 0 := by
    have : 0 < ‖v‖ := norm_pos_iff.mpr hv
    positivity
  rw [deficit_eq_sq_norm_sub_projection v w hv, sq_eq_zero_iff, norm_eq_zero,
    sub_eq_zero]
  constructor
  · intro h
    exact ⟨(⟪v, w⟫_ℝ) / ‖v‖ ^ 2, h⟩
  · rintro ⟨c, rfl⟩
    have h1 : ⟪v, c • v⟫_ℝ = c * ‖v‖ ^ 2 := by
      rw [real_inner_smul_right, real_inner_self_eq_norm_sq]
    rw [h1]
    field_simp

end IqgtGrassmannian
