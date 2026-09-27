/-
Intrinsic Quantum Geometric Tensor on the Grassmannian and Spectral Lower Bounds on Dissipation
(Jeromie Beasley, DOI 10.5281/zenodo.20768258): the exact finite-dimensional results.

* Corollary 3.1: a tangent vector at a projector (`V = ΠV + VΠ`, the derivative of `Π² = Π`) is
  purely off-diagonal: `ΠVΠ = 0`, `(1−Π)V(1−Π) = 0`, `V = ΠV(1−Π) + (1−Π)VΠ`.
* Theorem 3.2 (metric part): for `V = [[0, X], [X†, 0]]`, `Tr V² = 2 Tr(XX†)`.
* Section 2: `Q(V, W) = Tr(ΠVW)` equals the Frobenius pairing `Tr((VΠ)†(WΠ))` for self-adjoint
  `V` and an orthogonal projector `Π`; hence `Q` is Hermitian and `Q(V, V)` is real.
* Theorem 6.1 (Metric control of Berry curvature): the Gram matrix of `Q` is positive
  semidefinite, so `(Im Q(V,W))² ≤ Re Q(V,V) · Re Q(W,W) − (Re Q(V,W))²`, i.e.
  `|Ω| ≤ 2√(g(V,V)g(W,W) − g(V,W)²)` with `g = Re Q`, `Ω = 2 Im Q`.
* Corollary 6.1: metric rigidity `g(V,V) = 0` forces `Ω(V,W) = 0`.
* The two-level saturation step: `|A|²|B|² − (A·B)² = |A × B|²` in `ℝ³`.
* The residue step of Lemma 3.1: `1/((z−a)(z−b)) = (1/(z−a) − 1/(z−b))/(a−b)`.
-/
import Mathlib

namespace IqgtGrassmannian

open Matrix ComplexConjugate
open scoped ComplexOrder

variable {n : Type*} [Fintype n] [DecidableEq n]

/-! ## Corollary 3.1: tangent vectors are off-diagonal -/

/-- **Corollary 3.1.** If `Π² = Π` and `V = ΠV + VΠ`, then `ΠVΠ = 0`, `(1−Π)V(1−Π) = 0` and
`V = ΠV(1−Π) + (1−Π)VΠ`. -/
theorem tangent_offdiag (P V : Matrix n n ℂ) (hP : P * P = P) (hV : V = P * V + V * P) :
    P * V * P = 0 ∧ (1 - P) * V * (1 - P) = 0 ∧ V = P * V * (1 - P) + (1 - P) * V * P := by
  have h1 : P * V * P = 0 := by
    have e : P * V * P = P * (P * V + V * P) * P := by rw [← hV]
    have e2 : P * (P * V + V * P) * P = P * V * P + P * V * P := by
      simp only [mul_add, add_mul, ← mul_assoc, hP]
      rw [mul_assoc (P * V) P P, hP]
    rw [e2] at e
    have : P * V * P + P * V * P - P * V * P = 0 := by rw [← e, sub_self]
    simpa using this
  refine ⟨h1, ?_, ?_⟩
  · have : (1 - P) * V * (1 - P) = V - (P * V + V * P) + P * V * P := by
      simp only [sub_mul, mul_sub, one_mul, mul_one]
      abel
    rw [this, ← hV, h1, sub_self, add_zero]
  · have : P * V * (1 - P) + (1 - P) * V * P = P * V + V * P - (P * V * P + P * V * P) := by
      simp only [sub_mul, mul_sub, one_mul, mul_one]
      abel
    rw [this, h1, add_zero, sub_zero, ← hV]

/-! ## Theorem 3.2: the metric of an off-diagonal tangent vector -/

/-- **Theorem 3.2, metric.** For `V = [[0, X], [X†, 0]]`, `Tr V² = 2 Tr(XX†)`. -/
theorem trace_sq_offdiag {m k : Type*} [Fintype m] [Fintype k] (X : Matrix m k ℂ) :
    trace (fromBlocks 0 X Xᴴ 0 * fromBlocks 0 X Xᴴ 0) = 2 * trace (X * Xᴴ) := by
  rw [fromBlocks_multiply]
  simp only [zero_mul, zero_add, mul_zero, add_zero]
  have tb : ∀ (A : Matrix m m ℂ) (B : Matrix m k ℂ) (C : Matrix k m ℂ) (D : Matrix k k ℂ),
      trace (fromBlocks A B C D) = trace A + trace D := by
    intro A B C D
    simp only [trace, diag, Fintype.sum_sum_type, fromBlocks_apply₁₁, fromBlocks_apply₂₂]
  rw [tb, trace_mul_comm Xᴴ X]
  ring

/-! ## The intrinsic tensor as a Frobenius pairing -/

/-- The intrinsic tensor `Q_Π(V, W) = Tr(ΠVW)`. -/
def Q (P V W : Matrix n n ℂ) : ℂ := trace (P * V * W)

/-- `Q(V, W) = Tr((VΠ)†(WΠ))` for self-adjoint `V` and an orthogonal projector `Π`. -/
theorem Q_eq (P V W : Matrix n n ℂ) (hP : Pᴴ = P) (hPP : P * P = P) (hV : Vᴴ = V) :
    Q P V W = trace ((V * P)ᴴ * (W * P)) := by
  unfold Q
  rw [conjTranspose_mul, hV, hP, ← mul_assoc, trace_mul_comm (P * V * W) P, ← mul_assoc,
    ← mul_assoc, hPP]

/-- The Frobenius pairing is Hermitian: `Tr(B†A) = conj Tr(A†B)`. -/
theorem frob_herm (A B : Matrix n n ℂ) : trace (Bᴴ * A) = conj (trace (Aᴴ * B)) := by
  rw [← Complex.star_def, ← trace_conjTranspose, conjTranspose_mul, conjTranspose_conjTranspose]

/-- The Frobenius pairing as a sum of entries. -/
theorem frob_sum (A B : Matrix n n ℂ) :
    ∑ p : n × n, conj (A p.1 p.2) * B p.1 p.2 = trace (Aᴴ * B) := by
  simp only [trace, diag, mul_apply, conjTranspose_apply, Fintype.sum_prod_type, Complex.star_def]
  exact Finset.sum_comm

/-- The two tangent vectors `VΠ`, `WΠ` as the columns of one matrix. -/
def cols (A B : Matrix n n ℂ) : Matrix (n × n) (Fin 2) ℂ := of fun p a => ![A, B] a p.1 p.2

theorem gram_entries (A B : Matrix n n ℂ) :
    ((cols A B)ᴴ * cols A B) 0 0 = trace (Aᴴ * A) ∧ ((cols A B)ᴴ * cols A B) 0 1 = trace (Aᴴ * B) ∧
      ((cols A B)ᴴ * cols A B) 1 0 = trace (Bᴴ * A) ∧
        ((cols A B)ᴴ * cols A B) 1 1 = trace (Bᴴ * B) := by
  simp only [mul_apply, conjTranspose_apply, cols, of_apply]
  refine ⟨?_, ?_, ?_, ?_⟩ <;> simp [frob_sum]

/-- **The Gram determinant is nonnegative** (positive semidefiniteness of `Q`):
`Q(V,V) Q(W,W) − Q(V,W) Q(W,V)` has nonnegative real part. -/
theorem gram_det_nonneg (A B : Matrix n n ℂ) :
    0 ≤ trace (Aᴴ * A) * trace (Bᴴ * B) - trace (Aᴴ * B) * trace (Bᴴ * A) := by
  have h := (posSemidef_conjTranspose_mul_self (cols A B)).det_nonneg
  rw [det_fin_two] at h
  obtain ⟨h00, h01, h10, h11⟩ := gram_entries A B
  rwa [h00, h01, h10, h11] at h

/-- **Theorem 6.1 (Metric control of Berry curvature).** With `g = Re Q` and `Ω = 2 Im Q`,
`(Im Q(V,W))² ≤ g(V,V) g(W,W) − g(V,W)²`, i.e. `|Ω| ≤ 2√(g(V,V)g(W,W) − g(V,W)²)`. -/
theorem metric_control (P V W : Matrix n n ℂ) (hP : Pᴴ = P) (hPP : P * P = P) (hV : Vᴴ = V)
    (hW : Wᴴ = W) :
    (Q P V W).im ^ 2 ≤ (Q P V V).re * (Q P W W).re - (Q P V W).re ^ 2 := by
  set A := V * P
  set B := W * P
  have eVV : Q P V V = trace (Aᴴ * A) := Q_eq P V V hP hPP hV
  have eWW : Q P W W = trace (Bᴴ * B) := Q_eq P W W hP hPP hW
  have eVW : Q P V W = trace (Aᴴ * B) := Q_eq P V W hP hPP hV
  have hBA : trace (Bᴴ * A) = conj (trace (Aᴴ * B)) := frob_herm A B
  have hAA : (trace (Aᴴ * A)).im = 0 := by
    have := frob_herm A A
    exact Complex.conj_eq_iff_im.mp this.symm
  have hBB : (trace (Bᴴ * B)).im = 0 := by
    have := frob_herm B B
    exact Complex.conj_eq_iff_im.mp this.symm
  have hd := gram_det_nonneg A B
  rw [hBA] at hd
  have hre := (Complex.nonneg_iff.mp hd).1
  rw [eVV, eWW, eVW]
  simp only [Complex.sub_re, Complex.mul_re, Complex.conj_re, Complex.conj_im, hAA, hBB,
    mul_zero, sub_zero] at hre
  nlinarith [hre]

/-- **Corollary 6.1.** Metric rigidity `g(V,V) = 0` forces `Ω(V,W) = 0`, given `g(W,W) ≥ 0`. -/
theorem rigidity_kills_curvature (P V W : Matrix n n ℂ) (hP : Pᴴ = P) (hPP : P * P = P)
    (hV : Vᴴ = V) (hW : Wᴴ = W) (h0 : (Q P V V).re = 0) : (Q P V W).im = 0 := by
  have h := metric_control P V W hP hPP hV hW
  rw [h0, zero_mul, zero_sub] at h
  nlinarith [sq_nonneg (Q P V W).re, sq_nonneg (Q P V W).im]

/-! ## Two-level saturation and the residue step -/

/-- **Two-level saturation.** Lagrange's identity `|A|²|B|² − (A·B)² = |A × B|²` in `ℝ³`, which
makes the bound an equality for a two-band operator. -/
theorem lagrange_identity (a1 a2 a3 b1 b2 b3 : ℝ) :
    (a1 ^ 2 + a2 ^ 2 + a3 ^ 2) * (b1 ^ 2 + b2 ^ 2 + b3 ^ 2) - (a1 * b1 + a2 * b2 + a3 * b3) ^ 2 =
      (a2 * b3 - a3 * b2) ^ 2 + (a3 * b1 - a1 * b3) ^ 2 + (a1 * b2 - a2 * b1) ^ 2 := by
  ring

/-- **Lemma 3.1, residue step.** `1/((z−a)(z−b)) = (1/(z−a) − 1/(z−b))/(a−b)`, whose residue
at `z = a` gives the inverse-gap weight `1/(λ_m − λ_n)`. -/
theorem partial_fraction (z a b : ℂ) (hza : z ≠ a) (hzb : z ≠ b) (hab : a ≠ b) :
    1 / ((z - a) * (z - b)) = (1 / (z - a) - 1 / (z - b)) / (a - b) := by
  have h1 : z - a ≠ 0 := sub_ne_zero.mpr hza
  have h2 : z - b ≠ 0 := sub_ne_zero.mpr hzb
  have h3 : a - b ≠ 0 := sub_ne_zero.mpr hab
  field_simp
  ring

end IqgtGrassmannian
