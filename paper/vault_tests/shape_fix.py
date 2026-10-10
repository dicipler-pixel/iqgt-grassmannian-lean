# SCRIPT: IQGT-OLD-WORK-BATTERY-01b  (A2/A3 follow-up: which Gram matrix is the rotation-invariant shape)
# Kill: for G = X^T X / Tr, the kernel of X->G must be exactly scale + physical rotations X->RX (dim 4), rank 5.
import numpy as np
rng=np.random.default_rng(1729); h=1e-6
X=rng.normal(size=(3,3)); X/=np.linalg.norm(X)
sv=lambda W: np.array([W[0,0],W[1,1],W[2,2],W[0,1],W[0,2],W[1,2]])
def jac(f):
    J=[]
    for i in range(9):
        E=np.zeros(9);E[i]=1;E=E.reshape(3,3); J.append((sv(f(X+h*E))-sv(f(X-h*E)))/(2*h))
    return np.array(J).T
Gs=[np.array(g,float) for g in ([[0,-1,0],[1,0,0],[0,0,0]],[[0,0,-1],[0,0,0],[1,0,0]],[[0,0,0],[0,0,-1],[0,1,0]])]
for name,f in (('XX^T',lambda X:(X@X.T)/np.trace(X@X.T)),('X^T X',lambda X:(X.T@X)/np.trace(X.T@X))):
    J=jac(f); s=np.linalg.svd(J,compute_uv=False); rk=int((s>1e-6*s[0]).sum())
    left=max(np.linalg.norm(J@(G@X).ravel()) for G in Gs); right=max(np.linalg.norm(J@(X@G).ravel()) for G in Gs)
    print(f'{name}: rank {rk} | physical rotation X->RX moves it by {left:.2e} | Jacobi relabel X->XR moves it by {right:.2e} | scale {np.linalg.norm(J@X.ravel()):.1e}')
A=X@X.T; B=X.T@X
print('same spectrum:',np.allclose(np.linalg.eigvalsh(A),np.linalg.eigvalsh(B)))
