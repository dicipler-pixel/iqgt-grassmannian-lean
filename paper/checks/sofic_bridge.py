# SCRIPT: IQGT-SOFIC-BRIDGE-V1
# Kill conditions (before run):
# S1 his identity d_H(U,V) = (1/n)||P_U-P_V||_F^2 = (2/n) sum sin^2 theta on 300 random permutation pairs: kill if > 1e-12
# S2 along U(t)=exp(t log U) the intrinsic metric g = (1/2)Tr(dP/dt^2) is constant and the IQGT length equals
#    sqrt(sum theta_j^2) (geodesic): kill if rel err > 1e-6
# S3 arrangement: five disjoint swaps vs one 10-cycle on n=20 points have EQUAL Hamming (sum sin^2) and must give
#    DIFFERENT intrinsic lengths; kill if equal within 1e-9 (then the intrinsic metric forgets arrangement too)
# S4 sandwich d_H/2 <= (1/n) sum theta^2 <= pi^2/8 d_H holds on 2000 random pairs: kill on violation
# S5 conjugate permutations (same cycle type) have equal intrinsic distance from the identity (isometry): kill if not
# S6 sign: on a 2-parameter unitary family, Omega_IQGT = i Tr(P[dP_a,dP_b]) vs article form -i Tr(P[.,.]) : report which
#    matches Berry curvature F = dA, A = i<u|du> (two-level reference), so the two papers can be put on one convention.
import numpy as np
from scipy.linalg import expm, logm
rng=np.random.default_rng(1); res={}
def perm(p):
    n=len(p); U=np.zeros((n,n)); U[p,np.arange(n)]=1; return U
def P_of(U):
    n=U.shape[0]; Q=np.vstack([np.eye(n),U])/np.sqrt(2); return Q@Q.conj().T
def angles(U,V):
    s=np.linalg.svd((np.eye(len(U))+U.conj().T@V)/2,compute_uv=False); return np.arccos(np.clip(s,0,1))
w1=0
for t in range(300):
    n=rng.integers(2,13); U=perm(rng.permutation(n)); V=perm(rng.permutation(n))
    dH=np.mean(np.any(U!=V,axis=0)); th=angles(U,V)
    w1=max(w1,abs(dH-np.linalg.norm(P_of(U)-P_of(V))**2/n),abs(dH-2*np.sum(np.sin(th)**2)/n))
res['S1']=(w1,w1<1e-12)
def length(U,steps=400):
    L=logm(U); ts=np.linspace(0,1,steps+1); gs=[]
    for t in ts[:-1]:
        h=1e-6; dP=(P_of(expm((t+h)*L))-P_of(expm((t-h)*L)))/(2*h); gs.append(0.5*np.trace(dP@dP).real)
    return np.sqrt(np.mean(gs)), np.ptp(gs)/np.mean(gs)
w2=0; 
for t in range(20):
    n=rng.integers(3,9); U=perm(rng.permutation(n)); 
    if np.allclose(U,np.eye(n)): continue
    Lg,var=length(U); th=angles(np.eye(n),U); w2=max(w2,abs(Lg-np.sqrt(np.sum(th**2)))/np.sqrt(np.sum(th**2)),var)
res['S2']=(w2,w2<1e-6)
n=20; swaps=list(range(n)); 
for i in range(0,10,2): swaps[i],swaps[i+1]=swaps[i+1],swaps[i]
cyc=list(range(n)); cyc[:10]=list(range(1,10))+[0]
A,B=perm(swaps),perm(cyc); I=np.eye(n)
hA,hB=np.mean(np.any(A!=I,axis=0)),np.mean(np.any(B!=I,axis=0))
tA,tB=angles(I,A),angles(I,B)
res['S3']=((hA,hB,round(np.sum(np.sin(tA)**2),6),round(np.sum(np.sin(tB)**2),6),round(np.sqrt(np.sum(tA**2)),6),round(np.sqrt(np.sum(tB**2)),6)),abs(np.sum(tA**2)-np.sum(tB**2))>1e-9)
w4=0
for t in range(2000):
    n=rng.integers(2,13); U=perm(rng.permutation(n)); V=perm(rng.permutation(n))
    dH=np.mean(np.any(U!=V,axis=0)); th=angles(U,V); m=np.sum(th**2)/n
    if not (dH/2-1e-12<=m<=np.pi**2/8*dH+1e-12): w4+=1
res['S4']=(w4,w4==0)
p=rng.permutation(8); s=rng.permutation(8); U=perm(p); S=perm(s); V=S@U@S.T
res['S5']=(abs(np.sum(angles(I[:8,:8],U)**2)-np.sum(angles(I[:8,:8],V)**2)),abs(np.sum(angles(I[:8,:8],U)**2)-np.sum(angles(I[:8,:8],V)**2))<1e-10)
# S6 sign on the Bloch sphere: rank-1 P(theta,phi), Berry curvature F_theta_phi = (1/2) sin theta for the state |n> (upper)
def Pn(th,ph):
    u=np.array([np.cos(th/2),np.exp(1j*ph)*np.sin(th/2)]); return np.outer(u,u.conj()),u
th,ph,h=0.9,0.4,1e-5
P,u=Pn(th,ph); Pa=(Pn(th+h,ph)[0]-Pn(th-h,ph)[0])/(2*h); Pb=(Pn(th,ph+h)[0]-Pn(th,ph-h)[0])/(2*h)
plus=(1j*np.trace(P@(Pa@Pb-Pb@Pa))).real; minus=-plus
# reference: A = i<u|du>, F = d_th A_ph - d_ph A_th
def A(th,ph,which):
    u=Pn(th,ph)[1]
    if which=='th': du=(Pn(th+h,ph)[1]-Pn(th-h,ph)[1])/(2*h)
    else: du=(Pn(th,ph+h)[1]-Pn(th,ph-h)[1])/(2*h)
    return (1j*np.vdot(u,du)).real
F=(A(th+h,ph,'ph')-A(th-h,ph,'ph'))/(2*h)-(A(th,ph+h,'th')-A(th,ph-h,'th'))/(2*h)
res['S6']=((round(plus,6),round(minus,6),round(F,6)), True)
for k,(v,o) in res.items(): print(k,'PASS' if o else 'FAIL',v)
