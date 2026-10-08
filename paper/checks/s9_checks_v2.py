# SCRIPT: IQGT-REBUILD-S9-CHECKS-02
# Section 9 (mixed states). Kill conditions per line.
import numpy as np
from scipy.linalg import eigh, solve_continuous_lyapunov, sqrtm
rng=np.random.default_rng(909)
def herm(n,s=1.0): A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return s*(A+A.conj().T)/2
def proj(H,k): w,U=np.linalg.eigh(H); return U[:,:k]@U[:,:k].conj().T
R=[]
def rep(tag,val,thr,note=''):
    ok=val<thr; R.append(ok); print(('PASS ' if ok else 'FAIL ')+f'{tag}: {val:.3e} (kill if >= {thr:g}) {note}')
def exact_tangent(H,D,k):
    w,U=eigh(H); A=U.conj().T@D@U; n=len(w); M=np.zeros((n,n),complex)
    for a in range(n):
        for b in range(n):
            if (a<k)!=(b<k): M[a,b]=A[a,b]/((w[a]-w[b]) if a<k else (w[b]-w[a]))
    return U@M@U.conj().T
n=6; H=herm(n); D1,D2=herm(n),herm(n)
# M1/M2: k=1, rho=Pi, L=2V
P=proj(H,1); V=exact_tangent(H,D1,1); W=exact_tangent(H,D2,1); LV,LW=2*V,2*W
QP=np.trace(P@V@W); gP=0.5*np.trace(V@W).real; OP=(1j*np.trace(P@(V@W-W@V))).real
Qr=0.25*np.trace(P@LV@LW); gr=(np.trace(P@(LV@LW+LW@LV))/8).real; Or=(0.25j*np.trace(P@(LV@LW-LW@LV))).real
rep('M1 k=1: Q_rho = 1/4 Tr(rho L L) reduces exactly to Q_Pi (and g, Omega)',abs(Qr-QP)+abs(gr-gP)+abs(Or-OP),1e-12)
rep('M2 without the 1/4, Tr(rho L_V L_W) = 4 Q_Pi (normalisation, kill if not exactly 4)',abs(np.trace(P@LV@LW)/QP-4),1e-12)
# M3 k>1: rho = Pi/k gives Q_Pi/k
k=3; P3=proj(H,k); V3=exact_tangent(H,D1,k); W3=exact_tangent(H,D2,k)
rep('M3 k=3: with rho = Pi/k the reduction carries 1/k',abs(0.25*np.trace((P3/k)@(2*V3)@(2*W3))-np.trace(P3@V3@W3)/k),1e-12)
# M4 full-rank state: g_rho = 1/8 Tr(rho{L,L}) = F_Q/4 = Bures metric from fidelity
A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); rho=A@A.conj().T; rho/=np.trace(rho).real
Dr=herm(n,0.3); Dr-=np.trace(Dr)/n*np.eye(n)                         # traceless direction
L=solve_continuous_lyapunov(rho,2*Dr)                                 # rho L + L rho = 2 drho
rep('M4a SLD solves drho = 1/2(L rho + rho L)',np.abs(0.5*(L@rho+rho@L)-Dr).max(),1e-10)
lam,U=eigh(rho); B=U.conj().T@Dr@U
FQ=sum(2*abs(B[i,j])**2/(lam[i]+lam[j]) for i in range(n) for j in range(n))
g8=(np.trace(rho@(L@L+L@L))/8).real
rep('M4b g_rho = 1/8 Tr(rho{L,L}) = F_Q/4',abs(g8-FQ/4)/g8,1e-10)
eps=1e-4; r1=rho-0.5*eps*Dr; r2=rho+0.5*eps*Dr   # symmetric points cancel the O(eps^3) term; v01 used one-sided points and tripped at 3.1e-3
sq=sqrtm(r1); fid=np.trace(sqrtm(sq@r2@sq)).real
dB2=2*(1-fid)
rep('M4c Bures distance^2 = g_rho * dt^2 to leading order',abs(dB2/eps**2-g8)/g8,1e-3,f'{dB2/eps**2:.6f} vs {g8:.6f}')
# M5 SLD existence for singular rho: requires (I-Pi) V (I-Pi) = 0
Pk=proj(H,2); rho_s=Pk/2
def sld_residual(Vv):
    # least-squares solve of 1/2(L rho + rho L) = Vv for L
    N=n*n; M=np.zeros((N,N),complex)
    for idx in range(N):
        E=np.zeros(N,complex); E[idx]=1; Y=E.reshape(n,n); M[:,idx]=(0.5*(Y@rho_s+rho_s@Y)).reshape(-1)
    sol,res,rk,sv=np.linalg.lstsq(M,Vv.reshape(-1),rcond=None); return np.linalg.norm(M@sol-Vv.reshape(-1))
Vt=exact_tangent(H,D1,2)/2
Q=np.eye(n)-Pk; Vbad=Vt+Q@herm(n)@Q
rep('M5a tangent variation: SLD equation solvable (residual 0)',sld_residual(Vt),1e-10)
rep('M5b opposite: a kernel-kernel block makes it unsolvable (kill if residual ~ 0)',1e-6/sld_residual(Vbad),1.0,f'residual {sld_residual(Vbad):.3f}')
# M6 Kronecker form of the regularised SLD solve: A = 1/2(I (x) rho + rho^T (x) I) + gamma I, column stacking
gam=0.05; X=herm(n); I=np.eye(n)
Amat=0.5*(np.kron(I,rho_s)+np.kron(rho_s.T,I))+gam*np.eye(n*n)
Y=np.linalg.solve(Amat,X.reshape(-1,order='F')).reshape(n,n,order='F')
rep('M6 Kronecker solve returns Y with 1/2(rho Y + Y rho) + gamma Y = X (works for singular rho)',np.abs(0.5*(rho_s@Y+Y@rho_s)+gam*Y-X).max(),1e-10)
# M7 the i*eta prescription is not a tangent vector (not self-adjoint)
w,U=eigh(H); Aeig=U.conj().T@D1@U; eta=0.1; M=np.zeros((n,n),complex)
for a in range(n):
    for b in range(n):
        if (a<1)!=(b<1): M[a,b]=Aeig[a,b]/(w[a]-w[b]+1j*eta)
Veta=U@M@U.conj().T
rep('M7 V_eta is not self-adjoint (kill if it is)',1e-6/np.abs(Veta-Veta.conj().T).max(),1.0,f'|V_eta - V_eta^dag| = {np.abs(Veta-Veta.conj().T).max():.3f}')
print('SUMMARY',sum(R),'/',len(R))
# M8 (added): gamma -> 0 limit of (9.4) returns an SLD for a tangent variation of a singular rho
Vt2=exact_tangent(H,D1,2)/2
for gm in (1e-3,1e-6,1e-9):
    Am=0.5*(np.kron(I,rho_s)+np.kron(rho_s.T,I))+gm*np.eye(n*n)
    Yg=np.linalg.solve(Am,Vt2.reshape(-1,order='F')).reshape(n,n,order='F')
    print(f'     gamma={gm:.0e}: residual of (9.1) = {np.abs(0.5*(rho_s@Yg+Yg@rho_s)-Vt2).max():.2e}, |Y| = {np.linalg.norm(Yg):.3f}')
rep('M8 gamma=1e-9 solve satisfies (9.1) for a tangent variation, with bounded Y',np.abs(0.5*(rho_s@Yg+Yg@rho_s)-Vt2).max(),1e-6)
print('SUMMARY (with M8)',sum(R),'/',len(R))
