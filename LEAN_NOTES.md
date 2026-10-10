# Which Lean version of each object to use — Mathlib or yours

In this repository the notes' `IQGT.*` statements live in `IqgtGrassmannian/Rebuild.lean` (namespace `IqgtGrassmannian.Rebuild`) and `IQGTCore` is `IqgtGrassmannian/Core.lean`.

Read against Mathlib master at commit 47d0979a (8 Oct 2026) and your public repos upg-lean (v4.33.0, 104 theorems) and offset-lean.
Rule used: take Mathlib's object when it is the standard one and your statement is a special case of it; keep yours when
it carries your other work (oblique idempotents, block couplings, the feedback/memory language) that Mathlib does not have.

| Object in the paper | Mathlib | Yours | Use for IQGT | Why |
|---|---|---|---|---|
| Projector Π | `IsStarProjection p` (self-adjoint idempotent, any star ring; `Algebra/Star/StarProjection`) and `IsIdempotentElem p` | `P * P = P` in any `Ring` (UPGAlgebra) | Mathlib `IsIdempotentElem` for the ring identities; `IsStarProjection` wherever self-adjointness is used | Your idempotent form is strictly more general (oblique idempotents, used in UPGDiracRefraction); keep it there. IQGT's tangent law needs only idempotence, so `IsIdempotentElem` is the exact fit and costs nothing. |
| Grassmannian Gr(k,ℋ) | `Module.Grassmannian` (Kenny Lau 2025): quotients of a module, functor of points, for algebraic geometry | none | Neither as a manifold; IQGT works pointwise with "projectors of trace k" | Mathlib's Grassmannian has no metric or tangent bundle and uses the quotient convention. A smooth/Riemannian Grassmannian is not in Mathlib. |
| Tangent law V = ΠV + VΠ, off-diagonal blocks | not present | implicit in block forms (`fromBlocks` with zero diagonal blocks, UPGBlockFeedback) | IQGT (`IsTangent`, both files) | New in this library, stated in any ring; the block form is your UPG picture. |
| Spectral projector (Riesz) | no Riesz/holomorphic projection; has `Matrix.IsHermitian.eigenvalues`, `eigenvectorBasis`, spectral theorem, and the continuous functional calculus `cfc` | none | numerics only for now; the Mathlib route is `cfc` of an indicator that is continuous on the spectrum (gap) | The contour integral (2.4) has no Mathlib API; `cfc` gives the same projector when the gap separates the spectrum. Not formalised yet. |
| Redistribution operator F | not present | `UPGAlgebra.redistributionOp`, sum and commutator forms, F=0 ⇔ [Ω,Π]=0, trace F = 0, complex block forms | yours for the UPG statements; IQGT restates it (`IQGT.redist`) and adds: F is tangent, double-bracket form | Your version is richer (trace, complex block energy identity) and ties to UPG; IQGT keeps its own copy so the paper is standalone, with the same definition byte-for-byte in meaning. |
| Trace identities | `Matrix.trace_mul_cycle`, `trace_mul_comm`, `trace_add` | `trace_fromBlocks_diag`, block energy identities (UPGAlgebra, UPGProjectorMetric) | Mathlib lemmas | Standard; IQGT's Remark 2.3 is three Mathlib rewrites. |
| Metric g = ½Tr(VW) as coupling energy | Frobenius norm (`Matrix.frobeniusNormedAddCommGroup`) | `half_projectorCrossEnergy_eq_couplingEnergy`: ½ cross energy = Σ B² (real), complex version in UPGAlgebra | yours | Your statement is exactly "g = ‖X‖²_HS" in block language, already kernel-checked. |
| Cauchy–Schwarz for Theorem 6.1 | `norm_inner_le_norm` (InnerProductSpace) on `EuclideanSpace ℂ ι` | none | Mathlib | Standard; the block of a tangent vector is a vector in `EuclideanSpace ℂ (r × h)`. |
| Friction band / length bound | `Finset.sum_le_sum`, `Finset.mul_sum`, `Finset.sum_mul_sq_le_sq_mul_sq` | none | Mathlib | Plain finite-sum inequalities. |
| Fubini–Study / quantum Fisher information / Kubo–Mori metric | not present | not present | IQGT (algebraic core only: SLD solution, F_Q = 4g/k trace step) | No Mathlib API; the analytic definitions would be new library work. |
| Double-bracket (Brockett) flow | not present | not present | not formalised | ODE on the Grassmannian; needs manifold/ODE infrastructure. |

Not to be sent to Mathlib: nothing in this repository is proposed for Mathlib (standing rule).
