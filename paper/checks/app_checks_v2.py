# SCRIPT: IQGT-APPENDIX-CHECKS-V2
# v2: D2 inner-product order fixed (test bug); C3 split: spin-1 lowest band SATURATES (v1 kill fired), opposite run on a generic 3-band family must be strict
# Kill conditions written before the run.
# A1 second resolvent identity with the MINUS sign: R0 - R1 = -R0 dS R1; kill if residual > 1e-12 (and the + sign must fail)
# A2 Neumann condition sup_z ||dS R0(z)|| < 1 on the gap contour -> perturbed Riesz projector keeps rank k; kill if rank changes
# A3 unbounded-case gap bound ||V||_HS <= sqrt(2k)(a max|l_m| + b)/Delta for dS = a*Sigma + b*B (||B||<=1): kill on violation
# B1 two-band: g = |d1_perp . d2_perp|/(4|d|^2) form and Omega = d.(d1 x d2)/(2|d|^3): kill if residual > 1e-10
# B2 two-band saturates (6.1) for every pair: kill if |Omega| / bound differs from 1 by > 1e-10
# B3 Chern number of the lower band over the sphere = 1 (our sign): kill if |C - 1| > 1e-3
# C1 spin-1 lowest band: g = (1/2)|dd_perp|^2/|d|^2 for k=1 and k=2: kill if off by > 1e-8 (9/16 must fail)
# C2 g(Pi) = g(I - Pi) for random Pi (no flat band): kill if > 1e-12
# C3 three-band strict: some pair has |Omega| < bound by > 1%  (saturation is not forced beyond two bands)
# D1 Krylov: CG on [(H-l)^2 + eta^2] X = (I-Pi) dH|m> gives g_ij -> exact pair sum as eta -> 0; kill if rel err at eta=1e-6 > 1e-6
# D2 the same solves give Omega_ij = -2 Im sum_m <psi_i|X_j>; kill if rel err > 1e-6
# D3 cost: CG iterations per solve, k solves; reported
# G1 gap-closing refraction: on a family where only the gap varies across the interface, g_par sin^2(theta) is conserved
#    along a geodesic (Clairaut) and lambda_min-type candidates are not; kill if drift of g_par sin^2 > 1e-8 relative
import numpy as np
from scipy.sparse.linalg import cg, LinearOperator
from scipy.integrate import solve_ivp
rng=np.random.default_rng(3); res={}
def herm(n):
    a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return (a+a.conj().T)/2
def riesz(S,c,r,npts=400):
    n=S.shape[0]; P=np.zeros((n,n),complex)
    for t in np.linspace(0,2*np.pi,npts,endpoint=False):
        z=c+r*np.exp(1j*t); dz=1j*r*np.exp(1j*t)*(2*np.pi/npts)
        P+=np.linalg.solve(z*np.eye(n)-S,np.eye(n))*dz
    return P/(2j*np.pi)
# A1
n=6; S=herm(n); dS=herm(n)*0.1; z=0.3+0.7j; I=np.eye(n)
R0=np.linalg.inv(z*I-S); R1=np.linalg.inv(z*I-S-dS)
a1m=np.abs(R0-R1+R0@dS@R1).max(); a1p=np.abs(R0-R1-R0@dS@R1).max()
res['A1']=((a1m,a1p),a1m<1e-12 and a1p>1e-3)
# A2
w=np.linalg.eigvalsh(S); k=2; c=(w[0]+w[1])/2; r=(w[2]-w[0])/2*0.99+ (w[1]-w[0])/4
c=(w[0]+w[1])/2; r=(w[2]-c)*0.5+(c-w[0])*0.5
zs=[c+r*np.exp(1j*t) for t in np.linspace(0,2*np.pi,200,endpoint=False)]
sup=max(np.linalg.norm(dS@np.linalg.inv(zz*I-S),2) for zz in zs)
P1=riesz(S+dS,c,r); res['A2']=((round(sup,3),round(np.trace(P1).real,8)),sup<1 and abs(np.trace(P1).real-k)<1e-6)
# A3
worst=0
for t in range(500):
    n=7; k=rng.integers(1,4); S=herm(n); w,U=np.linalg.eigh(S); a=rng.uniform(0,0.1); b=rng.uniform(0,0.1)
    B=herm(n); B/=np.linalg.norm(B,2); dS=a*S+b*B
    P=U[:,:k]@U[:,:k].conj().T; X=U.conj().T@dS@U
    Xb=np.zeros((n,n),complex)
    for i in range(k):
        for j in range(k,n): Xb[i,j]=X[i,j]/(w[i]-w[j])
    V=U@(Xb+Xb.conj().T)@U.conj().T
    bound=np.sqrt(2*k)*(a*max(abs(w[:k]))+b)/(w[k]-w[k-1])
    worst=max(worst,np.linalg.norm(V)/bound)
res['A3']=(round(worst,4),worst<=1+1e-12)
# B
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.diag([1,-1]).astype(complex)
def Pm(d):
    dn=d/np.linalg.norm(d); return 0.5*(np.eye(2)-(dn[0]*sx+dn[1]*sy+dn[2]*sz))
def dP(f,d,v,h=1e-6): return (f(d+h*v)-f(d-h*v))/(2*h)
eb1=eb2=0
for t in range(200):
    d=rng.normal(size=3); v1=rng.normal(size=3); v2=rng.normal(size=3); P=Pm(d)
    V=dP(Pm,d,v1); W=dP(Pm,d,v2); dh=d/np.linalg.norm(d)
    perp=lambda v: v-(v@dh)*dh
    g12=0.5*np.trace(V@W).real; Om=(1j*np.trace(P@(V@W-W@V))).real
    eb1=max(eb1,abs(g12-perp(v1)@perp(v2)/(4*d@d)),abs(Om-d@np.cross(v1,v2)/(2*np.linalg.norm(d)**3)))
    g11=0.5*np.trace(V@V).real; g22=0.5*np.trace(W@W).real
    eb2=max(eb2,abs(abs(Om)/(2*np.sqrt(g11*g22-g12**2))-1))
res['B1']=(eb1,eb1<1e-8); res['B2']=(eb2,eb2<1e-7)
th=np.linspace(0,np.pi,401); ph=np.linspace(0,2*np.pi,401); tot=0
for i in range(len(th)-1):
    for j in range(len(ph)-1):
        tm,pm=(th[i]+th[i+1])/2,(ph[j]+ph[j+1])/2
        d=np.array([np.sin(tm)*np.cos(pm),np.sin(tm)*np.sin(pm),np.cos(tm)])
        vt=np.array([np.cos(tm)*np.cos(pm),np.cos(tm)*np.sin(pm),-np.sin(tm)]); vp=np.array([-np.sin(tm)*np.sin(pm),np.sin(tm)*np.cos(pm),0])
        tot+=d@np.cross(vt,vp)/2*(th[i+1]-th[i])*(ph[j+1]-ph[j])
res['B3']=(tot/(2*np.pi),abs(tot/(2*np.pi)-1)<1e-3)
# C spin-1
Sx=np.array([[0,1,0],[1,0,1],[0,1,0]])/np.sqrt(2); Sy=np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])/np.sqrt(2); Sz=np.diag([1,0,-1]).astype(complex)
def Pk(d,k):
    w,U=np.linalg.eigh(d[0]*Sx+d[1]*Sy+d[2]*Sz); return U[:,:k]@U[:,:k].conj().T
ec1=0; ratio916=[]
for t in range(100):
    d=rng.normal(size=3); v=rng.normal(size=3); dh=d/np.linalg.norm(d); vp=v-(v@dh)*dh
    for k in (1,2):
        V=dP(lambda x: Pk(x,k),d,v); g=0.5*np.trace(V@V).real
        ec1=max(ec1,abs(g-0.5*(vp@vp)/(d@d))/g); ratio916.append(g/((9/16)*(vp@vp)/(d@d)))
res['C1']=((ec1,round(np.mean(ratio916),4)),ec1<1e-6)
ec2=0
for t in range(100):
    n=6; S=herm(n); dS=herm(n); f=lambda s: (lambda w,U: U[:,:2]@U[:,:2].conj().T)(*np.linalg.eigh(S+s*dS))
    V=(f(1e-6)-f(-1e-6))/2e-6; g1=0.5*np.trace(V@V).real; g2=0.5*np.trace((-V)@(-V)).real
    Q=np.eye(n)-f(0); ec2=max(ec2,abs(np.trace(f(0)@V@V).real-np.trace(Q@V@V).real))
res['C2']=(ec2,ec2<1e-8)
sat=0
for t in range(200):
    d=rng.normal(size=3); v1=rng.normal(size=3); v2=rng.normal(size=3)
    P=Pk(d,1); V=dP(lambda x: Pk(x,1),d,v1); W=dP(lambda x: Pk(x,1),d,v2)
    g11=0.5*np.trace(V@V).real; g22=0.5*np.trace(W@W).real; g12=0.5*np.trace(V@W).real
    Om=(1j*np.trace(P@(V@W-W@V))).real; bd=2*np.sqrt(max(g11*g22-g12**2,0)); sat=max(sat,abs(1-abs(Om)/bd))
res['C3a_spin1_saturates']=(sat,sat<1e-6)
H0,H1,H2=herm(3),herm(3),herm(3)
def Pg(p): return (lambda w,U: U[:,:1]@U[:,:1].conj().T)(*np.linalg.eigh(H0+p[0]*H1+p[1]*H2))
strict=0
for t in range(200):
    p=rng.normal(size=2)*0.5; e1=np.array([1.,0]); e2=np.array([0,1.])
    P=Pg(p); V=dP(Pg,p,e1); W=dP(Pg,p,e2)
    g11=0.5*np.trace(V@V).real; g22=0.5*np.trace(W@W).real; g12=0.5*np.trace(V@W).real
    Om=(1j*np.trace(P@(V@W-W@V))).real; bd=2*np.sqrt(max(g11*g22-g12**2,0)); strict=max(strict,1-abs(Om)/bd)
res['C3b_generic3band_strict']=(round(strict,4),strict>0.01)
# D Krylov on a 2D tight-binding lattice with a parameter
L=14; N=L*L
def lattice(a1,a2):
    H=np.zeros((N,N),complex)
    for x in range(L):
        for y in range(L):
            i=x*L+y; j=((x+1)%L)*L+y; l=x*L+(y+1)%L
            H[i,j]+=-1*np.exp(1j*a1); H[i,l]+=-1*np.exp(1j*a2*x/L)
            H[i,i]+=0.6*np.cos(2*np.pi*x/L)+0.4*np.sin(2*np.pi*y/L+a1)
    H=H+H.conj().T-np.diag(np.diag(H)); return H
a=(0.31,0.47); H=lattice(*a); h=1e-6
D1=(lattice(a[0]+h,a[1])-lattice(a[0]-h,a[1]))/(2*h); D2=(lattice(a[0],a[1]+h)-lattice(a[0],a[1]-h))/(2*h)
w,U=np.linalg.eigh(H); k=3
gex=np.zeros((2,2)); Oex=0
Ds=[D1,D2]
for m in range(k):
    for n2 in range(k,N):
        A=[U[:,m].conj()@Dd@U[:,n2] for Dd in Ds]
        for i in range(2):
            for j in range(2): gex[i,j]+=(A[i]*np.conj(A[j])).real/(w[n2]-w[m])**2
        Oex+=-2*(A[0]*np.conj(A[1])).imag/(w[n2]-w[m])**2
P=U[:,:k]@U[:,:k].conj().T; Q=np.eye(N)-P
def krylov(eta):
    g=np.zeros((2,2)); O=0; its=[]
    for m in range(k):
        lm=w[m]; M=LinearOperator((N,N),matvec=lambda x,lm=lm:(H-lm*np.eye(N))@((H-lm*np.eye(N))@x)+eta**2*x,dtype=complex)
        psi=[Dd@U[:,m] for Dd in Ds]; X=[]
        for p in psi:
            cnt=[0]; sol,info=cg(M,Q@p,rtol=1e-13,maxiter=5000,callback=lambda xk:cnt.__setitem__(0,cnt[0]+1)); X.append(Q@sol); its.append(cnt[0])
        for i in range(2):
            for j in range(2): g[i,j]+=(np.vdot(psi[i],X[j])).real
        O+=-2*(np.vdot(psi[0],X[1])).imag
    return g,O,its
out=[]
for eta in (1e-1,1e-2,1e-3,1e-6):
    g,O,its=krylov(eta); out.append((eta,np.abs(g-gex).max()/np.abs(gex).max(),abs(O-Oex)/abs(Oex),int(np.mean(its))))
res['D1']=([(e,f'{r:.1e}',it) for e,r,_,it in out],out[-1][1]<1e-6)
res['D2']=([(e,f'{r:.1e}') for e,_,r,_ in out],out[-1][2]<1e-6)
# G1 gap-closing refraction: H = cos(Kx) sx + sin(Kx) sy + m(y) sz, gap 2 sqrt(1+m^2) varies with y only
K=3.0; mfun=lambda y: 2.5*y; dm=2.5
gpar=lambda y: K*K/(4*(1+mfun(y)**2)); dgpar=lambda y: -K*K*2*mfun(y)*dm/(4*(1+mfun(y)**2)**2)
gyy=lambda y: dm**2/(4*(1+mfun(y)**2)**2); dgyy=lambda y: -dm**2*4*mfun(y)*dm/(4*(1+mfun(y)**2)**3)
def rhs(t,s):
    x,y,vx,vy=s
    ax=-(dgpar(y)/gpar(y))*vx*vy
    ay=(dgpar(y)*vx**2-dgyy(y)*vy**2)/(2*gyy(y))
    return [vx,vy,ax,ay]
y0=0.8; th0=np.deg2rad(50); s0=[0,y0,np.sin(th0)/np.sqrt(gpar(y0)),-np.cos(th0)/np.sqrt(gyy(y0))]
sol=solve_ivp(rhs,(0,3),s0,rtol=1e-12,atol=1e-14,dense_output=True); tt=np.linspace(0,3,600); x,y,vx,vy=sol.sol(tt)
sp=np.sqrt(gpar(y)*vx**2+gyy(y)*vy**2); sn=np.sqrt(gpar(y))*np.abs(vx)/sp
c1=gpar(y)*sn**2; lm=np.minimum(gpar(y),gyy(y)); c2=lm*sn**2
res['G1']=((float(np.ptp(c1)/c1[0]),float(np.ptp(c2)/c2[0]),round(float(y.min()),4)),np.ptp(c1)/c1[0]<1e-8)
for kk,(v,ok) in res.items(): print(kk,'PASS' if ok else 'FAIL',v)
print(sum(ok for v,ok in res.values()),'/',len(res))
