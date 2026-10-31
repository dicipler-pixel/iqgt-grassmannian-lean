import IqgtGrassmannian.Basic
-- Doubling on ℤ/9ℤ has a six-cycle, so three doublings (×8) do not return every residue: 8·1 = 8 ≢ 1.
example : ∀ r : Fin 9, (8 * r.val) % 9 = r.val % 9 := by decide
