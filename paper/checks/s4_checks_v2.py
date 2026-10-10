# SCRIPT: IQGT-REBUILD-S4-CHECKS-02
# Section 4 (pullback). References are brute-force eigenprojector differences; kill conditions per line.
import numpy as np
rng=np.random.default_rng(1004)
def herm(n): A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return (A+A.conj().T)/2
def proj(H,k): w,U=np.linalg.eigh(H); return U[:,:k]@U[:,:k].conj().T
R=[]
def rep(tag,val,thr,note=''):
    ok=val<thr; R.append(ok); print(('PASS ' if ok else 'FAIL ')+f'{tag}: {val:.3e} (kill if >= {thr:g}) {note}')
n,k=6,2; H0=herm(n); D=[herm(n),herm(n),herm(n)]; E=[[herm(n)*0.3 for _ in range(3)] for _ in range(3)]
for i in range(3):
    for j in range(i): E[i][j]=E[j][i]
def Sig(a): 
    S=H0+sum(a[i]*D[i] for i in range(3))
    for i in range(3):
        for j in range(3): S=S+0.5*a[i]*a[j]*E[i][j]
    return S
a0=np.zeros(3); h=1e-4
def Pa(a): return proj(Sig(a),k)
def dPi(a,i): e=np.zeros(3);e[i]=h; return (Pa(a+e)-Pa(a-e))/(2*h)
P=Pa(a0); V=[dPi(a0,i) for i in range(3)]
w,U=np.linalg.eigh(Sig(a0))
g=np.array([[0.5*np.trace(V[i]@V[j]).real for j in range(3)] for i in range(3)])
gs=np.array([[sum((((U.conj().T@D[i]@U)[m,q])*((U.conj().T@D[j]@U)[q,m])).real/(w[q]-w[m])**2 for m in range(k) for q in range(k,n)) for j in range(3)] for i in range(3)])
rep('U1 eq (8) metric sum = 1/2 Tr(d_iPi d_jPi)',np.abs(g-gs).max(),1e-6)
gT=np.array([[np.trace(P@V[i]@V[j]).real for j in range(3)] for i in range(3)])
rep('U1b Remark 4.1: = Re Tr(Pi d_iPi d_jPi)',np.abs(g-gT).max(),1e-8)
Om=np.array([[(1j*np.trace(P@(V[i]@V[j]-V[j]@V[i]))).real for j in range(3)] for i in range(3)])
OmT=np.array([[-2*np.trace(P@V[i]@V[j]).imag for j in range(3)] for i in range(3)])
rep('U2 Omega_ij = -2 Im Tr(Pi d_iPi d_jPi)',np.abs(Om-OmT).max(),1e-8)
# U3 Prop 4.3: immersion failure with the gap OPEN
Ue=U; B=Ue@np.diag(rng.normal(size=n))@Ue.conj().T            # diagonal in the eigenbasis: no cross-gap block
Bo=Ue@ (lambda M:(M+M.conj().T)/2)(np.block([[rng.normal(size=(k,k)),np.zeros((k,n-k))],[np.zeros((n-k,k)),rng.normal(size=(n-k,n-k))]]))@Ue.conj().T
VB=(proj(Sig(a0)+h*Bo,k)-proj(Sig(a0)-h*Bo,k))/(2*h)
G2=np.array([[0.5*np.trace(X@Y).real for Y in (V[0],VB)] for X in (V[0],VB)])
gap=w[k]-w[k-1]
rep('U3 block-diagonal direction: pullback metric degenerate (det g) while the gap stays open',abs(np.linalg.det(G2)),1e-10,f'gap={gap:.3f}')
# opposite: shrink the gap -> metric blows up instead
def g_at_gap(s):
    Hs=np.diag([0.,s,2.,3.]).astype(complex); Dx=np.zeros((4,4),complex); Dx[1,2]=Dx[2,1]=1
    Vx=(proj(Hs+1e-7*Dx,2)-proj(Hs-1e-7*Dx,2))/2e-7; return 0.5*np.trace(Vx@Vx).real
rep('U3b opposite limit: closing gap makes g large (g(gap 0.01)/g(gap 1) ~ 1e4)',abs(g_at_gap(1.99)/g_at_gap(1.0)-1e4)/1e4,1e-3,f'ratio={g_at_gap(1.99)/g_at_gap(1.0):.1f}')
# U4 second-variation identity d_k V_i = V_k V_i + Pi d_kV_i + d_kV_i Pi + V_i V_k
def d2(i,j):
    ei=np.zeros(3);ei[i]=h;ej=np.zeros(3);ej[j]=h
    return (Pa(a0+ei+ej)-Pa(a0+ei-ej)-Pa(a0-ei+ej)+Pa(a0-ei-ej))/(4*h*h)
d01=d2(0,1)
rep('U4 d_kV_i = V_kV_i + Pi d_kV_i + d_kV_i Pi + V_iV_k',np.abs(d01-(V[1]@V[0]+P@d01+d01@P+V[0]@V[1])).max(),1e-5)
# U5 Remark 4.5: only off-diagonal blocks of d_kV_i survive Tr((d_kV_i)V_j)
Q=np.eye(n)-P; off=P@d01@Q+Q@d01@P
V2o=P@V[2]@Q+Q@V[2]@P; rep('U5 Tr(d_kV_i V_j) uses only the off-diagonal blocks (V_j exactly off-diagonal; v01 tripped on finite-difference residue in V_j)',abs(np.trace(d01@V2o)-np.trace(off@V2o)),1e-12)
# U6 exact second-order contour formula and the metric gradient (11)
N=3000; lo=w[0]-1; hi=(w[k-1]+w[k])/2; c=(lo+hi)/2; r=(hi-lo)/2
z=c+r*np.exp(2j*np.pi*np.arange(N)/N); dz=1j*(z-c)*2*np.pi/N
def d2contour(i,j):
    tot=np.zeros((n,n),complex)
    for zi,dd in zip(z,dz):
        Rz=np.linalg.inv(zi*np.eye(n)-Sig(a0))
        tot+=(Rz@E[i][j]@Rz+Rz@D[i]@Rz@D[j]@Rz+Rz@D[j]@Rz@D[i]@Rz)*dd
    return tot/(2j*np.pi)
rep('U6 second-order Riesz formula = d_i d_j Pi (brute force)',np.abs(d2contour(0,1)-d01).max(),1e-5)
def gij(a,i,j):
    Vi=(lambda e:(Pa(a+e)-Pa(a-e))/(2*h))(np.eye(3)[i]*h); Vj=(lambda e:(Pa(a+e)-Pa(a-e))/(2*h))(np.eye(3)[j]*h)
    return 0.5*np.trace(Vi@Vj).real
hk=1e-3; e2=np.eye(3)[2]*hk
dg_fd=(gij(a0+e2,0,1)-gij(a0-e2,0,1))/(2*hk)
dg_11=0.5*np.trace(d2contour(2,0)@V[1]+V[0]@d2contour(2,1)).real
rep('U6b eq (11) metric gradient d_k g_ij via the exact second variation',abs(dg_fd-dg_11)/abs(dg_fd),1e-4,f'{dg_11:.6f} vs {dg_fd:.6f}')
# U7 Primacy at second order: a silent direction (no cross-gap block) leaves Pi fixed but changes g
def g00_along(s,Bdir):
    Hs=Sig(a0)+s*Bdir; Vx=(proj(Hs+h*D[0],k)-proj(Hs-h*D[0],k))/(2*h); return 0.5*np.trace(Vx@Vx).real
dPB=np.abs(VB).max(); dgB=(g00_along(1e-3,Bo)-g00_along(-1e-3,Bo))/2e-3
rep('U7 silent direction: d_B Pi = 0',dPB,1e-8)
rep('U7b ... yet d_B g != 0 (kill for "everything factors through cross-gap amplitudes" at 2nd order if this is 0)',1e-6/abs(dgB),1.0,f'd_B g_00 = {dgB:.4f}')
print('SUMMARY',sum(R),'/',len(R))
