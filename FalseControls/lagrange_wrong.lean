import IqgtGrassmannian.Basic
-- Lagrange's identity has no cross term: for A = B = e₁ the left side is 0, not 1.
example : ((1 : ℝ) ^ 2 + 0 ^ 2 + 0 ^ 2) * (1 ^ 2 + 0 ^ 2 + 0 ^ 2) - (1 * 1 + 0 * 0 + 0 * 0) ^ 2 = 1 := by
  norm_num
