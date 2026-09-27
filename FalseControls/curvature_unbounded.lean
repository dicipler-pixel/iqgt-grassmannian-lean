import IqgtGrassmannian.Basic
-- Curvature cannot exceed the metric budget: with g(V,V) = g(W,W) = 1, g(V,W) = 0,
-- (Im Q)² = 4 would violate (Im Q)² ≤ 1.
example : (2 : ℝ) ^ 2 ≤ 1 * 1 - 0 ^ 2 := by norm_num
