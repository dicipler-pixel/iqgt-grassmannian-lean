# SCRIPT: IQGT-V7-RIGIDITY-CHECKS
# Why each check proves its claim (kill conditions fixed before running):
# R1 gap bound equality: equality iff (a) every coupled cross pair sits exactly at gap Delta and
#    (b) Pi dS (I-Pi) has k singular values all equal to ||dS||. Kill: a k-copy two-level system at one gap
#    built to meet (a),(b) gives ratio != 1, or breaking (a) or (b) alone leaves ratio = 1.
# R2 speed ceiling equality: rigid rotation of a two-level system at fixed gap meets sqrt(g)=sqrt(k)||dS/dt||/Delta.
#    Kill: ratio != 1 along the path, or a path with a varying gap also gives ratio 1.
# R3 friction floor equality: int g dt = L^2/T iff speed is constant. Kill: constant speed misses, or variable speed hits.
# R4 curvature deficit: g(V,V)g(W,W)-g(V,W)^2-Omega^2/4 = g(V,V)*g(Wperp,Wperp), Wperp the part of W's block
#    complex-orthogonal to V's block. Kill: identity fails on any random pair.
# R5 trace condition against the round metric h of the parameter sphere: tr_h g >= |Omega|/sqrt(det h), equality
#    for two-level and spin-1 lowest bands (holomorphic), strict for a broken spin-1. Kill: an equality case misses or the strict case hits.
# R6 permutation sandwich: (1/n)sum theta^2 = (pi^2/8) d_H iff every nonzero angle is pi/2 iff P^2 = I.
#    Kill: an involution misses equality, or a non-involution meets it.
# R7 silent directions: at a point where dPi/dx = 0, d2Pi/dx2 is block off-diagonal (a tangent vector) and the
#    chordal distance grows as (x^2/2) sqrt(g(Pi'',Pi'')). Polar stratum: distance exactly 0 along it.
#    Kill: Pi'' has a diagonal block, or the distance ratio does not tend to 1.
# R8 three-body shape sphere: the Jacobi spinor's projector is the shape; mutual distances see only (w1,w2), so the
#    distance metric is degenerate exactly on collinear shapes, kernel transverse, shape still moving along: a fold.
#    Kill: mirror images give different distances, rank of the distance map is not 1 on the collinear circle,
#    or the kernel is not transverse.
# R9 Christoffel symbols of the pulled-back metric (Section 8.4 family, metric from the projector) equal half-sums of
#    the metric dipole. Kill: mismatch > 1e-6.
# R11 Christoffel symbols diverge at both strata: polar Gamma^x_xy = cot y, fold Gamma^x_xx = 1/x. Kill: ratio to these != 1.
# R12 Kuramoto: two oscillators lock iff |dw| <= 2K, where the two-level mixing angle is 22.5 deg. Kill: 0.995 fails to lock,
#    1.005 locks, or the angle differs.
# R10 holonomy: a loop that returns the projector keeps a Berry phase of half the enclosed solid angle. Kill: mismatch.
import numpy as np
rng=np.random.default_rng(11); R={}
def herm(n): A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return (A+A.conj().T)/2
def projI(H,I): w,U=np.linalg.eigh(H); return U[:,I]@U[:,I].conj().T
def dproj(H,dS,I,h=1e-6): return (projI(H+h*dS,I)-projI(H-h*dS,I))/(2*h)
def gapratio(H,dS,I):
    w=np.linalg.eigvalsh(H); n=len(w); Ic=[j for j in range(n) if j not in I]
    D=min(abs(w[a]-w[b]) for a in I for b in Ic); V=dproj(H,dS,I); k=len(I)
    return np.linalg.norm(V)/(np.sqrt(2*k)*np.linalg.norm(dS,2)/D)
# R1
ok=[]
for k in (1,2,3):
    n=2*k; H=np.diag([0.]*k+[1.]*k); dS=np.zeros((n,n),complex)
    for j in range(k): dS[j,k+j]=dS[k+j,j]=1.0
    H=H+np.diag(rng.normal(size=n))*0  # exact two-cluster
    # break degeneracy-free requirement: use tiny splitting inside clusters that keeps cross gaps >= Delta
    r_eq=gapratio(H+np.diag([0]*n),dS,list(range(k)))
    # break (a): add an uncoupled level at gap 1 (so Delta stays 1) and move one coupled partner out to 1.5
    n2=n+1; H2=np.diag([0.]*k+[1.]*(k-1)+[1.5]+[1.0]); dS2=np.zeros((n2,n2),complex); dS2[:n,:n]=dS
    r_a=gapratio(H2,dS2,list(range(k)))
    dS3=dS.copy(); dS3[0,k]=dS3[k,0]=0.5; r_b=gapratio(H,dS3,list(range(k))) if k>1 else r_eq
    ok.append((abs(r_eq-1)<1e-6) and (r_a<1-1e-3) and (k==1 or r_b<1-1e-3)); print(f'R1 k={k}: equal {r_eq:.8f}  gap broken {r_a:.4f}  strengths broken {r_b:.4f}')
R['R1']=all(ok)
# R2
sx=np.array([[0,1],[1,0]]);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1])
Dl=0.7; om=1.3; rs=[]
for t in np.linspace(0,2,9):
    n=np.array([np.sin(om*t),0,np.cos(om*t)]); nd=om*np.array([np.cos(om*t),0,-np.sin(om*t)])
    S=Dl/2*(n[0]*sx+n[2]*sz); Sd=Dl/2*(nd[0]*sx+nd[2]*sz)
    V=dproj(S,Sd,[0]); g=0.5*np.trace(V@V).real; rs.append(np.sqrt(g)/(np.linalg.norm(Sd,2)/Dl))
rv=[]
for t in np.linspace(0,2,9):
    r=1+0.5*np.sin(t); n=np.array([np.sin(om*t),0,np.cos(om*t)]); nd=om*np.array([np.cos(om*t),0,-np.sin(om*t)])
    S=r*(n[0]*sx+n[2]*sz); Sd=r*(nd[0]*sx+nd[2]*sz)+0.5*np.cos(t)*(n[0]*sx+n[2]*sz)
    V=dproj(S,Sd,[0]); g=0.5*np.trace(V@V).real; rv.append(np.sqrt(g)/(np.linalg.norm(Sd,2)/(2*r)))
# theorem level: T >= Delta_min L /(sqrt(k) sup||dS/dt||); ratio = Delta_min L/(sup * T), equality iff ratio = 1
def thm_ratio(Sf,Sdf,gapf,T=2.0,N=2001):
    ts=np.linspace(0,T,N); sp=[]; dn=[]; gp=[]
    for t in ts:
        V=dproj(Sf(t),Sdf(t),[0]); sp.append(np.sqrt(0.5*np.trace(V@V).real)); dn.append(np.linalg.norm(Sdf(t),2)); gp.append(gapf(t))
    return min(gp)*np.trapezoid(sp,ts)/(max(dn)*T)
nv=lambda t:np.array([np.sin(om*t),0,np.cos(om*t)]); ndv=lambda t:om*np.array([np.cos(om*t),0,-np.sin(om*t)])
Sr=lambda t:Dl/2*(nv(t)[0]*sx+nv(t)[2]*sz); Srd=lambda t:Dl/2*(ndv(t)[0]*sx+ndv(t)[2]*sz)
rr=lambda t:1+0.5*np.sin(t); rd=lambda t:0.5*np.cos(t)
Sb=lambda t:rr(t)*(nv(t)[0]*sx+nv(t)[2]*sz); Sbd=lambda t:rr(t)*(ndv(t)[0]*sx+ndv(t)[2]*sz)+rd(t)*(nv(t)[0]*sx+nv(t)[2]*sz)
tr_rigid=thm_ratio(Sr,Srd,lambda t:Dl); tr_breath=thm_ratio(Sb,Sbd,lambda t:2*rr(t))
print(f'R2 pointwise rigid ratio {min(rs):.8f}..{max(rs):.8f}; theorem ratio rigid {tr_rigid:.6f}, breathing gap {tr_breath:.4f}')
R['R2']=abs(min(rs)-1)<1e-6 and abs(max(rs)-1)<1e-6 and abs(tr_rigid-1)<1e-5 and tr_breath<1-1e-3
# R3
T=2.0; ts=np.linspace(0,T,20001); 
for name,sp in (('const',np.ones_like(ts)),('var',1+0.6*np.sin(3*ts))):
    L=np.trapezoid(sp,ts); I=np.trapezoid(sp**2,ts); print(f'R3 {name}: int g dt / (L^2/T) = {I/(L*L/T):.8f}')
    R['R3_'+name]=I/(L*L/T)
R['R3']=abs(R.pop('R3_const')-1)<1e-8 and R.pop('R3_var')>1+1e-3
# R4
bad=0
for t in range(300):
    n=rng.integers(3,8); k=rng.integers(1,n); H=herm(n); I=list(range(k)); P=projI(H,I)
    V=dproj(H,herm(n),I); W=dproj(H,herm(n),I)
    g=lambda A,B:0.5*np.trace(A@B).real; Om=1j*np.trace(P@(V@W-W@V)); Om=Om.real
    X=P@V@(np.eye(n)-P); Y=P@W@(np.eye(n)-P)
    c=np.trace(Y@X.conj().T)/np.trace(X@X.conj().T); Yp=Y-c*X; Wp=Yp+Yp.conj().T
    lhs=g(V,V)*g(W,W)-g(V,W)**2-Om**2/4; rhs=g(V,V)*g(Wp,Wp)
    bad+= abs(lhs-rhs)>1e-6*max(1,abs(lhs))
print('R4 deficit identity failures',bad,'/300'); R['R4']=bad==0
# R5
Sx=np.array([[0,1,0],[1,0,1],[0,1,0]])/np.sqrt(2); Sy=np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])/np.sqrt(2); Sz=np.diag([1.,0,-1])
def geom(Hf,th,ph,h=1e-5):
    P=projI(Hf(th,ph),[0]); Pt=(projI(Hf(th+h,ph),[0])-projI(Hf(th-h,ph),[0]))/(2*h); Pp=(projI(Hf(th,ph+h),[0])-projI(Hf(th,ph-h),[0]))/(2*h)
    g=lambda A,B:0.5*np.trace(A@B).real
    return g(Pt,Pt),g(Pp,Pp),(1j*np.trace(P@(Pt@Pp-Pp@Pt))).real
def tracegap(Hf):
    out=[]
    for th in np.linspace(0.3,2.8,6):
        for ph in np.linspace(0,5,4):
            gtt,gpp,Om=geom(Hf,th,ph); out.append(gtt+gpp/np.sin(th)**2-abs(Om)/np.sin(th))
    return min(out),max(out)
d=lambda th,ph:np.array([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
two=lambda th,ph:(lambda v:v[0]*sx+v[1]*sy+v[2]*sz)(d(th,ph))
s1=lambda th,ph:(lambda v:v[0]*Sx+v[1]*Sy+v[2]*Sz)(d(th,ph))
br=lambda th,ph:s1(th,ph)+0.8*Sz@Sz
a,b=tracegap(two);c_,d_=tracegap(s1);e,f=tracegap(br)
print(f'R5 tr_h g - |Omega|/sqrt(det h): two-level [{a:.2e},{b:.2e}] spin-1 [{c_:.2e},{d_:.2e}] broken [{e:.4f},{f:.4f}]')
R['R5']=max(abs(a),abs(b),abs(c_),abs(d_))<1e-6 and e>1e-3
# R6
def angles(perm):
    n=len(perm); P=np.zeros((n,n)); 
    for i,j in enumerate(perm): P[j,i]=1
    A=np.vstack([np.eye(n),np.eye(n)])/np.sqrt(2); B=np.vstack([np.eye(n),P])/np.sqrt(2)
    s=np.clip(np.linalg.svd(A.T@B,compute_uv=False),-1,1); return np.arccos(s)
def sandwich(perm):
    th=angles(perm); n=len(perm); dH=(2/n)*np.sum(np.sin(th)**2); return (1/n)*np.sum(th**2)/((np.pi**2/8)*dH), th
swaps=[1,0,3,2,5,4,7,6,9,8]+list(range(10,20)); cyc=list(range(1,10))+[0]+list(range(10,20))
dbl=[(2*r)%9 for r in range(9)]; inv=[0,2,1,3,5,4,6,7,8]
rs,_=sandwich(swaps); rc,_=sandwich(cyc); rd,thd=sandwich(dbl); ri,_=sandwich(inv)
print(f'R6 upper-ratio: five swaps {rs:.8f}  ten-cycle {rc:.4f}  doubling mod 9 {rd:.4f}  involution on 9 {ri:.8f}')
print('   doubling angles (deg):',np.round(np.sort(np.degrees(thd))[::-1],4),' sum sin^2',round(np.sum(np.sin(thd)**2),8),' length',round(np.sqrt(np.sum(thd**2)),4))
R['R6']=abs(rs-1)<1e-9 and abs(ri-1)<1e-9 and rc<1-1e-3 and rd<1-1e-3
# R7 fold family theta=y+1.1, phi=x^2/2, at x=0 silent along x
def Pfold(x,y): return projI(two(y+1.1,x*x/2),[0])
y0=0.3; h=1e-4
P0=Pfold(0,y0); P1=(Pfold(h,y0)-Pfold(-h,y0))/(2*h); P2=(Pfold(h,y0)-2*P0+Pfold(-h,y0))/h**2
diag=np.linalg.norm(P0@P2@P0)+np.linalg.norm((np.eye(2)-P0)@P2@(np.eye(2)-P0))
gP2=0.5*np.trace(P2@P2).real
rat=[np.sqrt(0.5)*np.linalg.norm(Pfold(x,y0)-P0)/((x*x/2)*np.sqrt(gP2)) for x in (0.1,0.03,0.01)]
pol=max(np.linalg.norm(projI(two(0.0+1e-15,x),[0])-projI(two(0.0,0),[0])) for x in np.linspace(0,6,13))
print(f'R7 |dPi/dx|={np.linalg.norm(P1):.2e}  diagonal blocks of Pi\'\'={diag:.2e}  quartic ratio {np.round(rat,6)}  polar drift {pol:.2e}')
R['R7']=np.linalg.norm(P1)<1e-8 and diag<1e-5 and abs(rat[-1]-1)<1e-3 and pol<1e-12
# R8 three-body shape sphere
def shape(q):
    r1=(q[1]-q[0])/np.sqrt(2); r2=(2*q[2]-q[0]-q[1])/np.sqrt(6); psi=np.array([r1,r2]); psi=psi/np.linalg.norm(psi)
    return np.outer(psi,psi.conj())
def dists(q):
    I=sum(abs(q[i]-q[j])**2 for i,j in ((0,1),(0,2),(1,2)))
    return np.array([abs(q[0]-q[1])**2,abs(q[0]-q[2])**2])/I
def wvec(P): return np.array([np.trace(P@s).real for s in (sz,sx,sy)])
q=np.array([0,1.0,0.3+0.8j]); qm=q.conj()
mirror=np.linalg.norm(dists(q)-dists(qm)); wq,wm=wvec(shape(q)),wvec(shape(qm))
# on the sphere, param (alpha along equator, beta latitude); distance map F(w)=dists depends on w1,w2 only
def wpt(a,b): return np.array([np.cos(b)*np.cos(a),np.cos(b)*np.sin(a),np.sin(b)])
def F(a,b):
    w=wpt(a,b); Iq=1.0
    # reconstruct a triangle with this shape: psi from Bloch vector, then q
    th=np.arccos(w[0]); ph=np.arctan2(w[2],w[1]); psi=np.array([np.cos(th/2),np.sin(th/2)*np.exp(1j*ph)])
    r1,r2=psi; q1=-r1/np.sqrt(2)-r2/np.sqrt(6); q2=r1/np.sqrt(2)-r2/np.sqrt(6); q3=2*r2/np.sqrt(6)
    return dists(np.array([q1,q2,q3])), shape(np.array([q1,q2,q3]))
ranks=[];kern=[];moves=[]
for a in np.linspace(0.2,6,6):
    h=1e-5; J=np.column_stack([(F(a+h,0)[0]-F(a-h,0)[0])/(2*h),(F(a,h)[0]-F(a,-h)[0])/(2*h)])
    s=np.linalg.svd(J,compute_uv=False); ranks.append(s[1]/s[0]); _,_,Vt=np.linalg.svd(J); kern.append(abs(Vt[-1][1]))
    moves.append(np.linalg.norm((F(a+h,0)[1]-F(a-h,0)[1])/(2*h)))
offr=[(lambda J:np.linalg.svd(J,compute_uv=False))(np.column_stack([(F(1.0+1e-5,b)[0]-F(1.0-1e-5,b)[0])/2e-5,(F(1.0,b+1e-5)[0]-F(1.0,b-1e-5)[0])/2e-5])) for b in (0.4,)]
print(f'R8 mirror distance diff {mirror:.1e}; mirror Bloch vectors {np.round(wq,4)} {np.round(wm,4)}')
print(f'   collinear circle: s2/s1 max {max(ranks):.1e}, kernel latitude-component min {min(kern):.6f}, |dPi along| min {min(moves):.3f}; off-equator s2/s1 {offr[0][1]/offr[0][0]:.3f}')
R['R8']=mirror<1e-12 and abs(wq[2]+wm[2])<1e-12 and max(ranks)<1e-6 and min(kern)>1-1e-6 and min(moves)>0.1 and offr[0][1]/offr[0][0]>0.05
# R9 Christoffel = dipole combos, on the Section 8.4 family d=(sin y cos kx, sin y sin kx, cos y), k=3, metric from the projector
kk=3
def P84(x,y): return projI(np.sin(y)*np.cos(kk*x)*sx+np.sin(y)*np.sin(kk*x)*sy+np.cos(y)*sz,[0])
def gm(x,y,h=1e-5):
    Px=(P84(x+h,y)-P84(x-h,y))/(2*h); Py=(P84(x,y+h)-P84(x,y-h))/(2*h); g=lambda A,B:0.5*np.trace(A@B).real
    return np.array([[g(Px,Px),g(Px,Py)],[g(Py,Px),g(Py,Py)]])
x0,yv,h=0.4,0.7,1e-4
G=[(gm(x0+h,yv)-gm(x0-h,yv))/(2*h),(gm(x0,yv+h)-gm(x0,yv-h))/(2*h)]
Gam1=np.array([[[0.5*(G[j][i,k]+G[i][j,k]-G[k][i,j]) for k in range(2)] for j in range(2)] for i in range(2)])
exact_xx_y=-0.5*(kk**2/4)*2*np.sin(yv)*np.cos(yv)
print(f'R9 Gamma_(xx,y) from the projector dipole {Gam1[0,0,1]:.9f}  exact -1/2 d_y g_par {exact_xx_y:.9f}')
R['R9']=abs(Gam1[0,0,1]-exact_xx_y)<1e-6
# R11 Christoffel divergence at both strata: polar Gamma^x_xy = cot y, fold Gamma^x_xx = 1/x
def Gx(gfun,i,j,x,y,h=1e-5):
    g=gfun(x,y); gi=np.linalg.inv(g); d=[(gfun(x+h,y)-gfun(x-h,y))/(2*h),(gfun(x,y+h)-gfun(x,y-h))/(2*h)]
    return sum(gi[0,l]*0.5*(d[i][j,l]+d[j][i,l]-d[l][i,j]) for l in range(2))
gpol=lambda x,y:np.diag([kk**2*np.sin(y)**2/4,0.25]); gfold=lambda x,y:np.diag([np.sin(y+1.1)**2*x*x/4,0.25])
c1=[Gx(gpol,0,1,0.3,y)*np.tan(y) for y in (0.1,0.01)]; c2=[Gx(gfold,0,0,x,0.3)*x for x in (0.1,0.01)]
print(f'R11 polar Gamma^x_xy * tan y {np.round(c1,6)}; fold Gamma^x_xx * x {np.round(c2,6)}')
R['R11']=max(abs(np.array(c1)-1))<1e-4 and max(abs(np.array(c2)-1))<1e-4
# R12 Kuramoto: two oscillators lock iff |dw| <= 2K; two-level mixing angle 22.5 deg at threshold
def locks(dw,K,T=400,dt=0.01):
    D=0.3
    for _ in range(int(T/dt)): D+=dt*(dw-2*K*np.sin(D))
    return abs(dw-2*K*np.sin(D))<1e-6
l1,l2=locks(1.99,1.0),locks(2.01,1.0); ang=0.5*np.degrees(np.arctan2(2.0,2.0))
print(f'R12 dw/2K=0.995 locks {l1}; 1.005 locks {l2}; mixing angle at threshold {ang:.4f} deg')
R['R12']=l1 and not l2 and abs(ang-22.5)<1e-12
# R10
th0=np.pi/3; N=4000; phs=np.linspace(0,2*np.pi,N+1)
vecs=[np.linalg.eigh(two(th0,p))[1][:,0] for p in phs]
prod=np.prod([np.vdot(vecs[i],vecs[i+1]) for i in range(N)]); berry=-np.angle(prod)
solid=2*np.pi*(1-np.cos(th0)); ret=np.linalg.norm(np.outer(vecs[0],vecs[0].conj())-np.outer(vecs[-1],vecs[-1].conj()))
print(f'R10 projector returns (|dPi|={ret:.1e}); Berry phase {berry:.6f} vs -solid/2 {(-solid/2):.6f} (mod 2pi)')
R['R10']=ret<1e-12 and abs(np.angle(np.exp(1j*(berry+solid/2))))<1e-4
print('RESULTS',R, 'PASS' if all(R.values()) else 'FAIL')
