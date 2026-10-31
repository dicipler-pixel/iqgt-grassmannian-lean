import IqgtGrassmannian.Rebuild
-- Remark 2.3 needs V tangent. With P = V = W = diag(1,0) (not tangent), Tr(PVW) + Tr(PWV) = 2 but Tr(VW) = 1.
example : Matrix.trace (!![(1:ℚ),0;0,0] * !![1,0;0,0] * !![1,0;0,0]) +
    Matrix.trace (!![(1:ℚ),0;0,0] * !![1,0;0,0] * !![1,0;0,0]) =
    Matrix.trace (!![(1:ℚ),0;0,0] * !![1,0;0,0]) := by
  simp [Matrix.trace_fin_two]
