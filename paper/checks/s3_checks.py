# SCRIPT: IQGT-REBUILD-S3-CHECKS-01
# Section 3 (spectral realization). Kill conditions stated per line; derivatives by central differences of eigenprojectors (no paper formula used as reference).
import numpy as np
rng=np.random.default_rng(613)
def herm(n): A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return (A+A.conj().T)/2
def proj(H,k): w,U=np.linalg.eigh(H); return U[:,:k]@U[:,:k].conj().T
R=[]
def rep(tag,val,thr,note=''):
    ok=val<thr; R.append(ok); print(('PASS ' if ok else 'FAIL ')+f'{tag}: {val:.3e} (kill if >= {thr:g}) {note}')
n,k=7,2; H=herm(n); w,U=np.linalg.eigh(H); P=proj(H,k); D1,D2=herm(n),herm(n); h=1e-5
dP=lambda D:(proj(H+h*D,k)-proj(H-h*D,k))/(2*h)
V,W=dP(D1),dP(D2)
# T1 lemma, both cross cases, and zero blocks, in the eigenbasis
Ve=U.conj().T@V@U; A=U.conj().T@D1@U
pred=np.zeros((n,n),complex)
for m in range(n):
    for q in range(n):
        if (m<k)!=(q<k):
            lin,lout=(w[m],w[q]) if m<k else (w[q],w[m])
            pred[m,q]=A[m,q]/(lin-lout)
rep('T1 Lemma 3.1 entries <m|V|n> = <m|dS|n>/(l_in - l_out), zero otherwise',np.abs(Ve-pred).max(),1e-7)
# T2 contour first-order formula
N=4000; c=(w[0]+w[1])/2; r=(w[k]+w[k-1])/2-c+ ( (w[k]-w[k-1])/2*0 )
r=( (w[k-1]+w[k])/2 - (w[0]-1) )/2; c=(w[0]-1)+r
z=c+r*np.exp(2j*np.pi*np.arange(N)/N); dz=1j*(z-c)*2*np.pi/N
Vc=sum(np.linalg.inv(zi*np.eye(n)-H)@D1@np.linalg.inv(zi*np.eye(n)-H)*d for zi,d in zip(z,dz))/(2j*np.pi)
rep('T2 first-order resolvent integral = dPi',np.abs(Vc-V).max(),1e-7)
# T3 degeneracy inside the occupied block: basis choice irrelevant
Hd=np.diag([-1.,-1.,0.5,1.0,1.7,2.2,3.0]).astype(complex); Q=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))[0]; Hd=Q@Hd@Q.conj().T
wd,Ud=np.linalg.eigh(Hd); Ud2=Ud.copy(); G=np.linalg.qr(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))[0]; Ud2[:,:2]=Ud[:,:2]@G
def lemV(Ub):
    Ae=Ub.conj().T@D1@Ub; M=np.zeros((n,n),complex)
    for m in range(n):
        for q in range(n):
            if (m<k)!=(q<k):
                lin,lout=(wd[m],wd[q]) if m<k else (wd[q],wd[m]); M[m,q]=Ae[m,q]/(lin-lout)
    return Ub@M@Ub.conj().T
rep('T3 degenerate occupied pair: V independent of eigenbasis choice',np.abs(lemV(Ud)-lemV(Ud2)).max(),1e-12)
Vd=(proj(Hd+h*D1,k)-proj(Hd-h*D1,k))/(2*h)
rep('T3b and equals the true derivative at the degenerate point',np.abs(lemV(Ud)-Vd).max(),1e-7)
rep('T4 V self-adjoint',np.abs(V-V.conj().T).max(),1e-9)
# T5 gap bound ||V||_HS <= sqrt(2k)||dS||/Delta over many trials
worst=0
for _ in range(3000):
    Hs=herm(6); ks=rng.integers(1,5); ws=np.linalg.eigvalsh(Hs); Dl=ws[ks]-ws[ks-1]; Ds=herm(6)
    Vs=(proj(Hs+1e-6*Ds,ks)-proj(Hs-1e-6*Ds,ks))/2e-6
    worst=max(worst,np.linalg.norm(Vs)/(np.sqrt(2*ks)*np.linalg.norm(Ds,2)/Dl))
rep('T5 Prop 3.6 gap bound, worst ratio over 3000 trials (must stay <= 1)',worst,1.0+1e-6,f'ratio={worst:.4f}')
# opposite: saturating example k=1, 2 levels, dS purely off-diagonal
H2=np.diag([0.,1.]).astype(complex); D=np.array([[0,1],[1,0]],complex)
V2=(proj(H2+1e-6*D,1)-proj(H2-1e-6*D,1))/2e-6
rep('T5b opposite: two-level off-diagonal push saturates the bound (ratio 1)',abs(np.linalg.norm(V2)/(np.sqrt(2)*1/1)-1),1e-6)
# T6/T7 sums vs brute force
B=U.conj().T@D2@U
gs=sum((A[m,q]*B[q,m]).real/(w[q]-w[m])**2 for m in range(k) for q in range(k,n))
Os=1j*sum((A[m,q]*B[q,m]-B[m,q]*A[q,m])/(w[q]-w[m])**2 for m in range(k) for q in range(k,n))
g_true=0.5*np.trace(V@W).real; O_aug=(1j*np.trace(P@(V@W-W@V))).real
rep('T6 metric sum (general V,W) = 1/2 Tr(VW)',abs(gs-g_true),1e-7)
rep('T7 curvature sum = Omega with Omega = i Tr(Pi[V,W])',abs(Os.real-O_aug)+abs(Os.imag),1e-7)
print('     with the opposite sign definition the sum would be off by', abs(Os.real+O_aug),'(sign check)')
# T8 the numerator is gap-weighted, not the projected commutator
unw=1j*sum((A[m,q]*B[q,m]-B[m,q]*A[q,m]) for m in range(k) for q in range(k,n))
comm=1j*np.trace(P@(D1@D2-D2@D1)@ (np.eye(n)) )  # for reference only
rep('T8 weighted and unweighted numerators differ (kill if equal)',1e-3/abs(Os-unw/((w[k]-w[k-1])**2)),1.0)
# T9 onto: dS_mn=(l_m-l_n)X_mn generates the target X
X=np.zeros((n,n),complex); X[:k,k:]=rng.normal(size=(k,n-k))+1j*rng.normal(size=(k,n-k))
Vt=U@(X+X.conj().T)@U.conj().T
Dg=np.zeros((n,n),complex)
for m in range(k):
    for q in range(k,n): Dg[m,q]=(w[m]-w[q])*X[m,q]; Dg[q,m]=np.conj(Dg[m,q])
Dg=U@Dg@U.conj().T
rep('T9 Remark 3.10: dS_mn=(l_m-l_n)X_mn produces exactly V(X)',np.abs(dP(Dg)-Vt).max(),1e-7)
print('SUMMARY',sum(R),'/',len(R))
