import IqgtGrassmannian.Rebuild
-- P V P = 0 alone does not make V tangent: P = diag(1,0), V = diag(0,1) gives P V + V P = 0 ≠ V.
example : IqgtGrassmannian.Rebuild.IsTangent (!![(1:ℚ),0;0,0]) (!![0,0;0,1]) := by
  unfold IqgtGrassmannian.Rebuild.IsTangent; simp
