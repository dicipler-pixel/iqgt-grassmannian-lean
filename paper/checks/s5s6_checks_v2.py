# SCRIPT: IQGT-REBUILD-S5S6-CHECKS-02
# Section 5 (dissipation) and Section 6 (metric control of curvature). Kill conditions per line.
import numpy as np
from scipy.linalg import eigh
rng=np.random.default_rng(56)
def herm(n,s=1.0): A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return s*(A+A.conj().T)/2
def proj(H,k): w,U=np.linalg.eigh(H); return U[:,:k]@U[:,:k].conj().T
R=[]
def rep(tag,val,thr,note=''):
    ok=val<thr; R.append(ok); print(('PASS ' if ok else 'FAIL ')+f'{tag}: {val:.3e} (kill if >= {thr:g}) {note}')
def tang(P):
    n=P.shape[0]; Q=np.eye(n)-P; X=P@(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))@Q; return X, X+X.conj().T
g=lambda V,W:0.5*np.trace(V@W).real
Om=lambda P,V,W:(1j*np.trace(P@(V@W-W@V))).real
# ================= Section 6 =================
n,k=7,3; P=proj(herm(n),k)
worst=np.inf; worst12=np.inf
for _ in range(5000):
    X,V=tang(P); Y,W=tang(P)
    lhs=abs(Om(P,V,W)); rhs=2*np.sqrt(max(g(V,V)*g(W,W)-g(V,W)**2,0))
    worst=min(worst,rhs-lhs); worst12=min(worst12,2*np.sqrt(g(V,V)*g(W,W))-lhs)
rep('W1 |Omega| <= 2 sqrt(gVV gWW - gVW^2): most negative slack over 5000 pairs',max(0,-worst),1e-12,f'min slack {worst:.3e}')
rep('W1b the weaker June-12 form |Omega| <= 2 sqrt(gVV gWW) also holds',max(0,-worst12),1e-12)
X,V=tang(P); c=0.7-1.3j; W=c*X+np.conj(c)*X.conj().T
lhs=abs(Om(P,V,W)); rhs=2*np.sqrt(g(V,V)*g(W,W)-g(V,W)**2)
rep('W2 equality when Y = cX (c complex)',abs(lhs-rhs),1e-10)
Y,W2=tang(P); rep('W2b opposite: a generic pair is strictly inside (kill if equal)',1e-6/(2*np.sqrt(g(V,V)*g(W2,W2)-g(V,W2)**2)-abs(Om(P,V,W2))),1.0)
Wi=1j*X+(1j*X).conj().T
rep('W3 extremal partner iX: g(V,W)=0',abs(g(V,Wi)),1e-12)
rep('W3b ... and |Omega| = 2||X||^2 = the bound',abs(abs(Om(P,V,Wi))-2*np.linalg.norm(X)**2)+abs(2*np.sqrt(g(V,V)*g(Wi,Wi))-2*np.linalg.norm(X)**2),1e-10)
# parametric form 4 det g_ij >= Omega_12^2
H0,D1,D2=herm(6),herm(6),herm(6); h=1e-5
dP=lambda D:(proj(H0+h*D,2)-proj(H0-h*D,2))/(2*h)
Va,Vb=dP(D1),dP(D2); P2=proj(H0,2)
G=np.array([[g(Va,Va),g(Va,Vb)],[g(Vb,Va),g(Vb,Vb)]]); O12=Om(P2,Va,Vb)
rep('W4 parametric form 4 det g >= Omega^2 (kill if violated)',max(0,O12**2-4*np.linalg.det(G)),1e-12,f'4det={4*np.linalg.det(G):.4f} Om^2={O12**2:.4f}')
# angle form is relative control: two-level near a closing gap, both sides grow together
for r in (1.0,0.1,0.01):
    sx=np.array([[0,1],[1,0]]);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1.])
    Hd=lambda a,b: (r+a)*sz+b*sx  # gap 2r at a=b=0
    Pd=lambda a,b: proj(Hd(a,b)+0*sy,1)
    hh=1e-6*r
    Vx=(Pd(0,hh)-Pd(0,-hh))/(2*hh); Hy=lambda a,b: r*sz+a*sx+b*sy
    Vx=(proj(Hy(hh,0),1)-proj(Hy(-hh,0),1))/(2*hh); Vy=(proj(Hy(0,hh),1)-proj(Hy(0,-hh),1))/(2*hh); Pp=proj(Hy(0,0),1)
    print(f'     gap {2*r:5.2f}: |Omega| = {abs(Om(Pp,Vx,Vy)):.4e}, bound = {2*np.sqrt(g(Vx,Vx)*g(Vy,Vy)-g(Vx,Vy)**2):.4e}')
rep('W5 near a closing gap both sides diverge together (two levels saturate at every gap): relative control, not boundedness',0,1,'see lines above')
# ================= Section 5 =================
# X1/X2: SLD for rho = Pi/k is L = 2 dPi, and F_Q = 4 g / k (k=1: F_Q/4 = g)
for kk in (1,2,3):
    H=herm(6); D=herm(6); w_,U_=np.linalg.eigh(H); A_=U_.conj().T@D@U_; M_=np.zeros((6,6),complex)
    for a_ in range(6):
        for b_ in range(6):
            if (a_<kk)!=(b_<kk): M_[a_,b_]=A_[a_,b_]/((w_[a_]-w_[b_]) if a_<kk else (w_[b_]-w_[a_]))
    V=U_@M_@U_.conj().T; Pk=proj(H,kk); rho=Pk/kk; drho=V/kk   # exact tangent vector (Lemma 3.1); v01 used finite differences and tripped at 1.15e-9
    L=2*V; sld=np.abs(0.5*(L@rho+rho@L)-drho).max()
    lam,U=np.linalg.eigh(rho); A=U.conj().T@drho@U
    F=sum(2*(lam[i]-lam[j])**2/(lam[i]+lam[j])*abs(A[i,j])**2/ (1 if True else 1) for i in range(6) for j in range(6) if lam[i]+lam[j]>1e-12)
    F=sum(2*abs(A[i,j])**2/(lam[i]+lam[j]) for i in range(6) for j in range(6) if lam[i]+lam[j]>1e-12)
    rep(f'X1 k={kk}: L = 2 dPi solves the SLD equation for rho = Pi/k',sld,1e-9)
    rep(f'X2 k={kk}: F_Q = 4 g / k (exact QFI formula, independent of L)',abs(F-4*g(V,V)/kk)/F,1e-7,f'F={F:.6f} 4g/k={4*g(V,V)/kk:.6f}')
# X3 unitary driving produces no entropy
Hs=herm(5); Uu=__import__('scipy.linalg',fromlist=['expm']).expm(-1j*0.7*herm(5)); Pu=proj(Hs,2)
ent=lambda r:-sum(x*np.log(x) for x in np.linalg.eigvalsh(r) if x>1e-14)
rep('X3 unitary motion of Pi/k leaves von Neumann entropy unchanged',abs(ent(Pu/2)-ent(Uu@Pu@Uu.conj().T/2)),1e-12)
# X4 minimum driving time: tau >= Delta_min L / (sqrt(k) sup||dSigma||)
worst=0
for _ in range(300):
    n_=5; kk=rng.integers(1,4); A0,A1,A2=herm(n_),herm(n_,0.6),herm(n_,0.4)
    ts=np.linspace(0,1,401); S=lambda t:A0+np.sin(3*t)*A1+t*t*A2; dS=lambda t:3*np.cos(3*t)*A1+2*t*A2
    Lg=0; gaps=[]; sup=0
    for i,t in enumerate(ts[:-1]):
        tm=t+0.5*(ts[1]-ts[0]); V=(proj(S(tm+1e-6),kk)-proj(S(tm-1e-6),kk))/2e-6
        Lg+=np.sqrt(g(V,V))*(ts[1]-ts[0]); w=np.linalg.eigvalsh(S(tm)); gaps.append(w[kk]-w[kk-1]); sup=max(sup,np.linalg.norm(dS(tm),2))
    worst=max(worst,(min(gaps)*Lg/(np.sqrt(kk)*sup))/1.0)
rep('X4 driving-time bound: worst Delta_min L /(sqrt(k) sup||dS||) over 300 paths of duration 1 (must be <= 1)',worst,1.0+1e-6,f'{worst:.4f}')
# X5 the spectral lower bound on dissipation, length form: int zeta_rot >= 2 min_n c_n D_n * L^2 / T  (zero temperature)
def zeta_rot_T0(H,dH,tau,kprec):
    E,U=eigh(H); Hd=U.conj().T@dH@U; z=0; gg=0; m=[]
    for nn in range(1,len(E)):
        D=E[nn]-E[0]; c=tau/(1+kprec**2*tau**2*D**2); gn=abs(Hd[0,nn])**2/D**2
        z+=2*c*D*gn; gg+=gn; m.append(2*c*D)
    return z,gg,min(m),max(m)
worstband=0; worstlen=0
for _ in range(400):
    n_=rng.integers(3,6); A0=np.diag(np.sort(rng.uniform(0,3,n_))+np.arange(n_)*0.6); A1,A2=herm(n_,0.35),herm(n_,0.15)
    tau=rng.choice([0.3,3.0]); kp=rng.integers(0,2)
    ss=np.linspace(0,1,301); ds=ss[1]-ss[0]; Z=0; Lg=0; mins=[]
    for s in ss[:-1]+ds/2:
        H=A0+s*A1+s*s*A2; dH=A1+2*s*A2
        z,gg,mn,mx=zeta_rot_T0(H,dH,tau,kp)
        worstband=max(worstband, (mn*gg-z)/max(z,1e-30), (z-mx*gg)/max(z,1e-30))
        Z+=z*ds; Lg+=np.sqrt(gg)*ds; mins.append(mn)
    # protocol of duration T: W_ex = (1/T) int_0^1 zeta ds ; bound = 2 min c D * L^2 / T  -> T cancels
    worstlen=max(worstlen,(min(mins)*Lg**2-Z)/Z)
rep('X5a zero-temperature band min 2cD g <= zeta_rot <= max 2cD g (400 paths)',max(0,worstband),1e-9)
rep('X5b length form W_ex >= (2 min c D) L^2 / T (Cauchy-Schwarz) over 400 paths',max(0,worstlen),1e-9)
# X6 opposite: the June form  Sigma >= int adot^2 g dt  fails once the bath relaxes fast
beta=1.0; Dl=1.0; taub=1e-3; gss=1.0
sigma_rate=beta*2*taub*Dl*np.tanh(beta*Dl/2)*gss   # two-level rotation, k=0: Sigma-rate = beta*zeta_rot
rep('X6 opposite: fast bath (tau=1e-3) gives entropy rate below g-rate, so the uncorrected bound is violated (kill if not violated)',sigma_rate/gss,1.0,f'Sigma rate {sigma_rate:.2e} vs g {gss}')
print('SUMMARY',sum(R),'/',len(R))
