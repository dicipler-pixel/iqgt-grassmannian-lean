import IqgtGrassmannian.Basic
-- The off-diagonal block appears twice: Tr V² = 2 Tr(XX†), so for X = 1 (1×1) it is 2, not 1.
example : (2 : ℂ) * Matrix.trace ((1 : Matrix (Fin 1) (Fin 1) ℂ) * Matrix.conjTranspose 1) = 1 := by
  simp
