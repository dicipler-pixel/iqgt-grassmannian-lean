# SCRIPT: IQGT-SPIN1-TILT-CHECK
# Kill (before run): k=1 and k=2 metrics of spin-1 must stay equal once a term breaking the +-|d| symmetry
# is added (Sigma = d.S + c Sz^2), else the k=1 <-> k=2 equality needs that symmetry. Complement rule must still hold.
import numpy as np
rng=np.random.default_rng(8)
Sx=np.array([[0,1,0],[1,0,1],[0,1,0]])/np.sqrt(2); Sy=np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])/np.sqrt(2); Sz=np.diag([1,0,-1]).astype(complex)
def Pk(d,k,c):
    w,U=np.linalg.eigh(d[0]*Sx+d[1]*Sy+d[2]*Sz+c*Sz@Sz); return U[:,:k]@U[:,:k].conj().T
for c in (0.0,0.8):
    r=[];comp=0
    for t in range(50):
        d=rng.normal(size=3); v=rng.normal(size=3); h=1e-6
        g=[]
        for k in (1,2):
            V=(Pk(d+h*v,k,c)-Pk(d-h*v,k,c))/(2*h); g.append(0.5*np.trace(V@V).real)
            P=Pk(d,k,c); comp=max(comp,abs(np.trace(P@V@V).real-np.trace((np.eye(3)-P)@V@V).real))
        r.append(g[0]/g[1])
    print('c=',c,'g(k=1)/g(k=2): median %.4f, range %.3f-%.3f'%(np.median(r),min(r),max(r)),'| complement rule err %.1e'%comp)
