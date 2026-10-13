set_option linter.unusedSimpArgs false
/-
  SCRIPT: IQGT-LEAN-CORE (no Mathlib)
  Intrinsic Quantum Geometric Tensor — the ring-level identities of Sections 2, 4 and 5,
  checked in core Lean 4 with no library. A minimal noncommutative ring is declared here,
  so every statement holds for matrices, bounded operators and any other (unital) ring.
  Jeromie Beasley
-/
namespace IqgtGrassmannian.Core

class NCRing (R : Type) extends Add R, Mul R, Neg R where
  zero : R
  one : R
  add_assoc : ∀ a b c : R, a + b + c = a + (b + c)
  add_comm : ∀ a b : R, a + b = b + a
  zero_add : ∀ a : R, zero + a = a
  neg_add : ∀ a : R, -a + a = zero
  mul_assoc : ∀ a b c : R, a * b * c = a * (b * c)
  one_mul : ∀ a : R, one * a = a
  mul_one : ∀ a : R, a * one = a
  left_distrib : ∀ a b c : R, a * (b + c) = a * b + a * c
  right_distrib : ∀ a b c : R, (a + b) * c = a * c + b * c

variable {R : Type} [NCRing R]

instance : OfNat R 0 := ⟨NCRing.zero⟩
instance : OfNat R 1 := ⟨NCRing.one⟩
instance : Sub R := ⟨fun a b => a + -b⟩

open NCRing

theorem sub_def (a b : R) : a - b = a + -b := rfl
theorem zero_add' (a : R) : 0 + a = a := NCRing.zero_add a
theorem add_zero' (a : R) : a + 0 = a := by rw [NCRing.add_comm]; exact NCRing.zero_add a
theorem neg_add' (a : R) : -a + a = 0 := NCRing.neg_add a
theorem add_neg' (a : R) : a + -a = 0 := by rw [NCRing.add_comm]; exact NCRing.neg_add a
theorem one_mul' (a : R) : 1 * a = a := NCRing.one_mul a
theorem mul_one' (a : R) : a * 1 = a := NCRing.mul_one a

theorem add_left_cancel' {a b c : R} (h : a + b = a + c) : b = c := by
  have : -a + (a + b) = -a + (a + c) := by rw [h]
  rw [← NCRing.add_assoc, ← NCRing.add_assoc, neg_add', zero_add', zero_add'] at this
  exact this

theorem eq_zero_of_self_add {x : R} (h : x = x + x) : x = 0 := by
  have h2 : x + 0 = x + x := by rw [add_zero']; exact h
  exact (add_left_cancel' h2).symm

theorem mul_zero' (a : R) : a * 0 = 0 := by
  apply eq_zero_of_self_add
  rw [← NCRing.left_distrib, add_zero']

theorem zero_mul' (a : R) : 0 * a = 0 := by
  apply eq_zero_of_self_add
  rw [← NCRing.right_distrib, add_zero']

theorem neg_eq_of_add {a b : R} (h : a + b = 0) : -a = b := by
  have : -a + (a + b) = -a + 0 := by rw [h]
  rw [← NCRing.add_assoc, neg_add', zero_add', add_zero'] at this
  exact this.symm

theorem neg_mul' (a b : R) : -a * b = -(a * b) := by
  apply Eq.symm; apply neg_eq_of_add
  rw [← NCRing.right_distrib, add_neg', zero_mul']

theorem mul_neg' (a b : R) : a * -b = -(a * b) := by
  apply Eq.symm; apply neg_eq_of_add
  rw [← NCRing.left_distrib, add_neg', mul_zero']

theorem neg_neg' (a : R) : -(-a) = a := neg_eq_of_add (neg_add' a)

theorem neg_add_distrib (a b : R) : -(a + b) = -a + -b := by
  apply neg_eq_of_add
  calc a + b + (-a + -b) = a + (b + -a) + -b := by
        rw [NCRing.add_assoc, NCRing.add_assoc, NCRing.add_assoc]
    _ = a + (-a + b) + -b := by rw [NCRing.add_comm b (-a)]
    _ = 0 := by rw [← NCRing.add_assoc a (-a) b, add_neg', zero_add', add_neg']

theorem sub_mul' (a b c : R) : (a - b) * c = a * c - b * c := by
  rw [sub_def, sub_def, NCRing.right_distrib, neg_mul']

theorem mul_sub' (a b c : R) : a * (b - c) = a * b - a * c := by
  rw [sub_def, sub_def, NCRing.left_distrib, mul_neg']

theorem sub_self' (a : R) : a - a = 0 := add_neg' a

theorem eq_of_sub_eq_zero {a b : R} (h : a - b = 0) : a = b := by
  have : a + -b + b = 0 + b := by rw [← sub_def, h]
  rw [NCRing.add_assoc, neg_add', add_zero', zero_add'] at this
  exact this

/-! ## Section 2: the tangent law -/

/-- A projector is an idempotent. (For the paper it is also self-adjoint; nothing here needs it.) -/
def IsIdem (P : R) : Prop := P * P = P

/-- The linearised constraint obtained by differentiating `P * P = P`. -/
def IsTangent (P V : R) : Prop := V = P * V + V * P

theorem idem_left (P x : R) (hP : IsIdem P) : P * (P * x) = P * x := by
  rw [← NCRing.mul_assoc, hP]

theorem add_left_comm' (a b c : R) : a + (b + c) = b + (a + c) := by
  rw [← NCRing.add_assoc, NCRing.add_comm a b, NCRing.add_assoc]

/-- Equation (2.1), first half: `P V P = 0`. -/
theorem tangent_PVP (P V : R) (hP : IsIdem P) (hV : IsTangent P V) : P * V * P = 0 := by
  have hl : ∀ x : R, P * (P * x) = P * x := fun x => idem_left P x hP
  have hr : ∀ x : R, x * (P * P) = x * P := fun x => by rw [hP]
  apply eq_zero_of_self_add
  have e : P * V * P = P * (P * V + V * P) * P := congrArg (fun x => P * x * P) hV
  have e2 : P * (P * V + V * P) * P = P * V * P + P * V * P := by
    simp only [NCRing.left_distrib, NCRing.right_distrib, NCRing.mul_assoc, hl, hP, hr]
  exact e.trans e2

/-- Expansion used twice below: `(1-P) X (1-P) = X - (P X + X P) + P X P`. -/
theorem compl_sandwich (P X : R) : (1 - P) * X * (1 - P) = X + -(P * X + X * P) + P * X * P := by
  simp only [sub_def, NCRing.right_distrib, NCRing.left_distrib, neg_mul', mul_neg', one_mul', mul_one',
    neg_neg', neg_add_distrib, NCRing.mul_assoc, NCRing.add_assoc]

/-- Equation (2.1), second half: `(1 - P) V (1 - P) = 0`. -/
theorem tangent_QVQ (P V : R) (hP : IsIdem P) (hV : IsTangent P V) :
    (1 - P) * V * (1 - P) = 0 := by
  rw [compl_sandwich, tangent_PVP P V hP hV, add_zero']
  have : P * V + V * P = V := hV.symm
  rw [this, add_neg']

/-- Converse (Section 2): block off-diagonal elements satisfy the tangent law. -/
theorem tangent_of_blocks (P V : R) (h1 : P * V * P = 0) (h2 : (1 - P) * V * (1 - P) = 0) :
    IsTangent P V := by
  rw [compl_sandwich, h1, add_zero'] at h2
  have := eq_of_sub_eq_zero (show V - (P * V + V * P) = 0 from h2)
  exact this

/-- The tangent of the complement: if `V` moves `P`, then `-V` moves `1 - P` (so `g` agrees for `P` and `1-P`). -/
theorem tangent_compl (P V : R) (hV : IsTangent P V) : IsTangent (1 - P) (-V) := by
  unfold IsTangent at *
  have key : (1 - P) * (-V) + (-V) * (1 - P) = -V + -V + (P * V + V * P) := by
    simp only [sub_def, NCRing.right_distrib, NCRing.left_distrib, one_mul', mul_one', neg_mul', mul_neg',
      neg_neg', neg_add_distrib, NCRing.add_assoc]
    simp only [NCRing.add_assoc, NCRing.add_comm, add_left_comm']
  rw [key, ← hV, NCRing.add_assoc, neg_add', add_zero']

/-! ## Section 4.4: the redistribution operator -/

def redist (Ω P : R) : R := (1 - P) * Ω * P + P * Ω * (1 - P)

theorem redist_sum_form (Ω P : R) (hP : IsIdem P) :
    redist Ω P = Ω * P + P * Ω + -(P * Ω * P + P * Ω * P) := by
  have hl : ∀ x : R, P * (P * x) = P * x := fun x => idem_left P x hP
  unfold redist
  simp only [sub_def, NCRing.right_distrib, NCRing.left_distrib, neg_mul', mul_neg', one_mul', mul_one',
    neg_add_distrib, NCRing.mul_assoc, hl, hP, NCRing.add_assoc]
  simp only [NCRing.add_assoc, NCRing.add_comm, add_left_comm']

theorem redist_double_bracket (Ω P : R) (hP : IsIdem P) :
    redist Ω P = (Ω * P - P * Ω) * P - P * (Ω * P - P * Ω) := by
  have hl : ∀ x : R, P * (P * x) = P * x := fun x => idem_left P x hP
  have hr : ∀ x : R, x * (P * P) = x * P := fun x => by rw [hP]
  unfold redist
  simp only [sub_def, NCRing.right_distrib, NCRing.left_distrib, neg_mul', mul_neg', one_mul', mul_one',
    neg_neg', neg_add_distrib, NCRing.mul_assoc, hl, hP, hr, NCRing.add_assoc]
  simp only [NCRing.add_assoc, NCRing.add_comm, add_left_comm']

theorem proj_mul_compl (P : R) (hP : IsIdem P) : P * (1 - P) = 0 := by
  rw [mul_sub', mul_one', hP, sub_self']

theorem compl_mul_proj (P : R) (hP : IsIdem P) : (1 - P) * P = 0 := by
  rw [sub_mul', one_mul', hP, sub_self']

/-- The redistribution operator is a tangent vector at `P` (Section 4.4). -/
theorem redist_tangent (Ω P : R) (hP : IsIdem P) : IsTangent P (redist Ω P) := by
  unfold IsTangent redist
  have a1 : P * ((1 - P) * Ω * P) = 0 := by
    rw [← NCRing.mul_assoc, ← NCRing.mul_assoc, proj_mul_compl P hP, zero_mul', zero_mul']
  have a2 : P * (P * Ω * (1 - P)) = P * Ω * (1 - P) := by
    rw [← NCRing.mul_assoc, ← NCRing.mul_assoc, hP]
  have b1 : (1 - P) * Ω * P * P = (1 - P) * Ω * P := by
    rw [NCRing.mul_assoc _ P P, hP]
  have b2 : P * Ω * (1 - P) * P = 0 := by
    rw [NCRing.mul_assoc _ (1 - P) P, compl_mul_proj P hP, mul_zero']
  rw [NCRing.left_distrib, NCRing.right_distrib, a1, a2, b1, b2, zero_add', add_zero', NCRing.add_comm]

/-- `F(Ω,P) = 0` exactly when `Ω` commutes with `P` (Section 4.4). -/
theorem redist_eq_zero_iff (Ω P : R) (hP : IsIdem P) : redist Ω P = 0 ↔ Ω * P = P * Ω := by
  constructor
  · intro h
    unfold redist at h
    -- left-multiply by P:  P*F = P Ω (1-P)
    have hPF : P * ((1 - P) * Ω * P + P * Ω * (1 - P)) = P * Ω * (1 - P) := by
      rw [NCRing.left_distrib, ← NCRing.mul_assoc P ((1 - P) * Ω) P, ← NCRing.mul_assoc P (1 - P) Ω,
        proj_mul_compl P hP, zero_mul', zero_mul', zero_add', ← NCRing.mul_assoc P (P * Ω) (1 - P),
        ← NCRing.mul_assoc P P Ω, hP]
    have hFP : ((1 - P) * Ω * P + P * Ω * (1 - P)) * P = (1 - P) * Ω * P := by
      rw [NCRing.right_distrib, NCRing.mul_assoc ((1 - P) * Ω) P P, hP, NCRing.mul_assoc (P * Ω) (1 - P) P,
        compl_mul_proj P hP, mul_zero', add_zero']
    rw [h, mul_zero'] at hPF
    rw [h, zero_mul'] at hFP
    -- P Ω (1-P) = 0  gives  P Ω = P Ω P ;  (1-P) Ω P = 0  gives  Ω P = P Ω P
    have e1 : P * Ω - P * Ω * P = 0 := by
      rw [mul_sub', mul_one'] at hPF; exact hPF.symm
    have e2 : Ω * P - P * Ω * P = 0 := by
      rw [sub_mul', one_mul', sub_mul'] at hFP; exact hFP.symm
    calc Ω * P = P * Ω * P := eq_of_sub_eq_zero e2
      _ = P * Ω := (eq_of_sub_eq_zero e1).symm
  · intro h
    rw [redist_double_bracket Ω P hP, h, sub_self', zero_mul', mul_zero', sub_self']

/-- Lemma 5.1: `L = 2V` solves the symmetric-logarithmic-derivative equation for a moving projector:
    `(V + V) P + P (V + V) = V + V` (the factor 1/k of `ρ = P/k` cancels on both sides). -/
theorem sld_two_tangent (P V : R) (hV : IsTangent P V) : (V + V) * P + P * (V + V) = V + V := by
  have e : (V + V) * P + P * (V + V) = (P * V + V * P) + (P * V + V * P) := by
    simp only [NCRing.right_distrib, NCRing.left_distrib]
    simp only [NCRing.add_assoc, NCRing.add_comm, add_left_comm']
  rw [e, ← hV]

end IqgtGrassmannian.Core

#print axioms IqgtGrassmannian.Core.tangent_PVP
#print axioms IqgtGrassmannian.Core.tangent_QVQ
#print axioms IqgtGrassmannian.Core.tangent_of_blocks
#print axioms IqgtGrassmannian.Core.tangent_compl
#print axioms IqgtGrassmannian.Core.redist_sum_form
#print axioms IqgtGrassmannian.Core.redist_double_bracket
#print axioms IqgtGrassmannian.Core.redist_tangent
#print axioms IqgtGrassmannian.Core.sld_two_tangent
#print axioms IqgtGrassmannian.Core.redist_eq_zero_iff
