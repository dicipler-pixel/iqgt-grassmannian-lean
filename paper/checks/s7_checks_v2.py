# SCRIPT: IQGT-REBUILD-S7-CHECKS-02
# Section 7 (transport tensors). Brute-force references by finite differences of eigenprojectors / energies.
import numpy as np
from scipy.linalg import eigh
rng=np.random.default_rng(707)
def herm(n,s=1.0): A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return s*(A+A.conj().T)/2
def proj(H,k): w,U=np.linalg.eigh(H); return U[:,:k]@U[:,:k].conj().T
R=[]
def rep(tag,val,thr,note=''):
    ok=val<thr; R.append(ok); print(('PASS ' if ok else 'FAIL ')+f'{tag}: {val:.3e} (kill if >= {thr:g}) {note}')
# ---------- dipole: his two-index form (17) + rotation remainder = exact ----------
worst_full=0; worst_17=0
for trial in range(40):
    n=6; k=int(rng.integers(1,4)); H0=herm(n); D=[herm(n,0.6) for _ in range(3)]; E2=[[None]*3 for _ in range(3)]
    for a in range(3):
        for b in range(a,3): E2[a][b]=E2[b][a]=herm(n,0.3)
    S=lambda x: H0+sum(x[a]*D[a] for a in range(3))+0.5*sum(x[a]*x[b]*E2[a][b] for a in range(3) for b in range(3))
    def g(x,i,j,h=1e-5):
        e=np.eye(3)*h; Vi=(proj(S(x+e[i]),k)-proj(S(x-e[i]),k))/(2*h); Vj=(proj(S(x+e[j]),k)-proj(S(x-e[j]),k))/(2*h)
        return 0.5*np.trace(Vi@Vj).real
    i,j,kk=0,1,2; x0=np.zeros(3); hk=1e-3
    # reference: exact second variation by the Riesz contour (4.6), metric gradient (4.5); v01 used nested finite
    # differences and tripped at 4e-4 on step error
    ww_=np.linalg.eigvalsh(S(x0)); lo=ww_[0]-1; hi=(ww_[k-1]+ww_[k])/2; c=(lo+hi)/2; r=(hi-lo)/2; N=2000
    zs=c+r*np.exp(2j*np.pi*np.arange(N)/N); dzs=1j*(zs-c)*2*np.pi/N
    def d2c(a,b):
        tot=np.zeros((n,n),complex)
        for z_,dz_ in zip(zs,dzs):
            Rz=np.linalg.inv(z_*np.eye(n)-S(x0)); tot+=(Rz@E2[a][b]@Rz+Rz@D[a]@Rz@D[b]@Rz+Rz@D[b]@Rz@D[a]@Rz)*dz_
        return tot/(2j*np.pi)
    def d1c(a):
        tot=np.zeros((n,n),complex)
        for z_,dz_ in zip(zs,dzs):
            Rz=np.linalg.inv(z_*np.eye(n)-S(x0)); tot+=Rz@D[a]@Rz*dz_
        return tot/(2j*np.pi)
    G_fd=0.5*np.trace(d2c(kk,i)@d1c(j)+d1c(i)@d2c(kk,j)).real
    w,U=eigh(S(x0)); A=[U.conj().T@D[a]@U for a in range(3)]; B=lambda a,b:U.conj().T@E2[a][b]@U
    I=range(k); Ic=range(k,n)
    G17=0; Rem=0
    for m in I:
        for q in Ic:
            Dl=w[q]-w[m]; dDl=(A[kk][q,q]-A[kk][m,m]).real
            G17+=(( B(kk,i)[m,q]*A[j][q,m] + A[i][m,q]*B(kk,j)[q,m]).real)/Dl**2 - 2*((A[i][m,q]*A[j][q,m]).real)*dDl/Dl**3
            def rot(Ai,mm,qq):
                s=sum(A[kk][mm,l]*Ai[l,qq]/(w[mm]-w[l]) for l in range(n) if l!=mm)
                s+=sum(Ai[mm,l]*A[kk][l,qq]/(w[qq]-w[l]) for l in range(n) if l!=qq)
                return s
            Rem+=((rot(A[i],m,q)*A[j][q,m] + A[i][m,q]*rot(A[j],q,m)).real)/Dl**2
    worst_full=max(worst_full,abs(G17+Rem-G_fd)/abs(G_fd)); worst_17=max(worst_17,abs(G17-G_fd)/abs(G_fd))
rep('D1 exact dipole = two-index form (17) + eigenvector-rotation remainder (40 random families)',worst_full,1e-4)
rep('D2 opposite: (17) alone is NOT exact (kill if within 1e-3)',1e-3/max(worst_17,1e-30),1.0,f'worst rel error of (17) alone = {worst_17:.2f}')
# ---------- weights family ----------
wa=[0,0,0,0,0]
for t in range(200):
    n=int(rng.integers(3,7)); H0=herm(n); Hx=herm(n,0.5); Hy=herm(n,0.5); tau=10**rng.uniform(-1,1)
    E,U=eigh(H0); A=U.conj().T@Hx@U; B=U.conj().T@Hy@U; Dn=E[1:]-E[0]; gn=np.abs(A[0,1:])**2/Dn**2
    Om=-2*np.imag(A[0,1:]*B[1:,0])/Dn**2
    alpha=lambda x: np.sum(2*Dn**3/(Dn**2-x**2)*gn).real
    h=1e-4; e0=lambda a: eigh(H0+a*Hx)[0][0]
    wa[0]=max(wa[0],abs(alpha(0)+(e0(h)-2*e0(0)+e0(-h))/h**2)/alpha(0))
    def u0(a,b):
        ww,V=eigh(H0+a*Hx+b*Hy); v=V[:,0]; return v*np.exp(-1j*np.angle(U[:,0].conj()@v))
    def um(a,b,m):
        ww,V=eigh(H0+a*Hx+b*Hy); v=V[:,m]; return v*np.exp(-1j*np.angle(U[:,m].conj()@v))
    h2=1e-5  # v01 used 1e-4 here and tripped at 1.5e-6 on step error
    dx=(u0(h2,0)-u0(-h2,0))/(2*h2); dy=(u0(0,h2)-u0(0,-h2))/(2*h2)
    lhs=np.imag(dx.conj()@(H0-E[0]*np.eye(n))@dy); rhs=-0.5*np.sum(Dn*Om)
    wa[1]=max(wa[1],abs(lhs-rhs)/abs(rhs))
    # Berry connection polarizability, Gao-Yang-Niu / Kaplan convention: G = 2 Re sum_m A_0m A_m0 /(e0-em), A_nm = i<n|d u_m>
    Gx=0
    for m in range(1,n):
        dum=(um(h,0,m)-um(-h,0,m))/(2*h); A0m=1j*(U[:,0].conj()@dum)          # A_0m = i<0|d_x m>
        Gx+=2*np.real(A0m*np.conj(A0m))/(E[0]-E[m])                          # A_m0 = conj(A_0m)
    wa[2]=max(wa[2],abs(Gx-(-np.sum(2*gn/Dn)))/abs(Gx))
    z1=np.sum(2*tau/(1+tau**2*Dn**2)*Dn*gn); wa[3]=max(wa[3],abs(z1-tau*(alpha(0)-alpha(1j/tau)))/z1)
rep('P1 static response: alpha(0)=sum 2 D g_n equals -E0'' (finite differences)',wa[0],1e-5)
rep('P2 orbital-moment kernel Im<dx u|(H-E0)|dy u> = -1/2 sum D_n Omega_n',wa[1],1e-6)
rep('P3 Berry connection polarizability G_xx (Berry connections) = -sum 2 g_n/D_n',wa[2],1e-6)
rep('P4 friction = tau[alpha(0)-alpha(i/tau)] (precession bath)',wa[3],1e-10)
# ---------- gap-closing exponents ----------
def fam(Dl):
    H=np.diag([0,Dl,2.0+Dl,3.1]).astype(complex); X=np.zeros((4,4),complex); X[0,1]=X[1,0]=0.3; X[0,2]=X[2,0]=0.2; X[1,3]=X[3,1]=0.25
    Y=np.zeros((4,4),complex); Y[0,1]=-0.3j; Y[1,0]=0.3j; Y[0,3]=Y[3,0]=0.1
    Z=np.diag([0.,0.4,0.1,0.]).astype(complex); Z[1,2]=Z[2,1]=0.2
    return H,X,Y,Z
def quantities(Dl,beta=1.0,tau=1.0):
    H,X,Y,Z=fam(Dl); h=1e-6*Dl
    P=lambda M:proj(M,1)
    Vx=(P(H+h*X)-P(H-h*X))/(2*h); Vy=(P(H+h*Y)-P(H-h*Y))/(2*h)
    gg=0.5*np.trace(Vx@Vx).real; Om=(1j*np.trace(P(H)@(Vx@Vy-Vy@Vx))).real
    hz=1e-3*Dl
    def gxx(s):
        M=H+s*Z; Vv=(P(M+h*X)-P(M-h*X))/(2*h); return 0.5*np.trace(Vv@Vv).real
    dip=(gxx(hz)-gxx(-hz))/(2*hz)
    E,Uu=eigh(H); Ax=Uu.conj().T@X@Uu
    z0=sum(2*tau/(1+0*tau**2)*(E[n_]-E[0])*abs(Ax[0,n_])**2/(E[n_]-E[0])**2 for n_ in range(1,4))
    p=np.exp(-beta*E); p/=p.sum()
    zT=sum(tau*abs(Ax[m,n_])**2*(p[m]-p[n_])/(E[n_]-E[m]) for m in range(4) for n_ in range(4) if m!=n_ and (m==0)!=(n_==0))
    return gg,abs(Om),abs(dip),z0,zT
Ds=np.array([0.04,0.02,0.01,0.005])
Q=np.array([quantities(d) for d in Ds])
sl=[np.polyfit(np.log(Ds),np.log(Q[:,c]),1)[0] for c in range(5)]
print('     slopes (g, Omega, dipole, zeta_T0, zeta_T>0):',np.round(sl,3))
target=[-2,-2,-3,-1,0]
rep('X1 gap-closing exponents -2,-2,-3,-1,0',max(abs(a-b) for a,b in zip(sl,target)),0.05)
# ---------- eta regularization ----------
H,X,Y,Z=fam(0.0+1e-9); E,Uu=eigh(H); Ax=Uu.conj().T@X@Uu
def g_eta(eta): return sum(abs(Ax[0,n_])**2/((E[n_]-E[0])**2+eta**2) for n_ in range(1,4))
bound=lambda eta: sum(abs(Ax[0,n_])**2 for n_ in range(1,4))/eta**2
rep('E1 at a closed gap the regularized metric stays below sum|A|^2/eta^2 (eta=0.1)',max(0,g_eta(0.1)-bound(0.1)),1e-12,f'g_eta={g_eta(0.1):.3f} bound={bound(0.1):.3f}')
H2,X2,_,_=fam(0.5); E2_,U2=eigh(H2); A2=U2.conj().T@X2@U2
g0=sum(abs(A2[0,n_])**2/(E2_[n_]-E2_[0])**2 for n_ in range(1,4)); ge=sum(abs(A2[0,n_])**2/((E2_[n_]-E2_[0])**2+1e-8**2) for n_ in range(1,4))
rep('E2 eta -> 0 recovers the metric on an open gap',abs(ge-g0)/g0,1e-12)
print('SUMMARY',sum(R),'/',len(R))
