# SCRIPT: IQGT-REBUILD-S2-CHECKS-02
# Kill conditions are stated before each test; a test FAILS if its kill condition trips.
import numpy as np, scipy.linalg as sl
rng=np.random.default_rng(20260612)
def herm(n): A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return (A+A.conj().T)/2
def proj(H,k): w,U=np.linalg.eigh(H); return U[:,:k]@U[:,:k].conj().T
def tangent(P):
    n=P.shape[0]; X=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); Q=np.eye(n)-P
    V=P@X@Q; return V+V.conj().T
R=[]
def rep(name,val,thr,kill):
    ok=val<thr; R.append((name,val,thr,ok)); print(f"{'PASS' if ok else 'FAIL'} {name}: {val:.3e} (kill if >= {thr:g}) -- {kill}")
n,k=7,2; H=herm(n); P=proj(H,k); Qc=np.eye(n)-P
# K1 tangent law both directions
Ps=[proj(H+t*herm(n)*0+t*np.diag(np.arange(n)),k) for t in (0,)]
D=herm(n); h=1e-6
Pt=lambda t: proj(H+t*D,k)
V=(Pt(h)-Pt(-h))/(2*h)
rep("K1a dPi satisfies V=PV+VP",np.linalg.norm(V-P@V-V@P),1e-6,"a true projector velocity violating the linearised law")
rep("K1b diagonal blocks vanish",np.linalg.norm(P@V@P)+np.linalg.norm(Qc@V@Qc),1e-6,"nonzero intra-subspace block")
W=tangent(P); rep("K1c every off-diagonal V obeys the law",np.linalg.norm(W-P@W-W@P),1e-12,"converse fails")
# K2 surjectivity (Lemma 2.1): Pi(t)=e^{itA}Pi e^{-itA}, A=[[0,iX],[-iX*,0]]
Xb=P@(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))@Qc; Vt=Xb+Xb.conj().T
A=1j*Xb-1j*Xb.conj().T
Pc=lambda t: sl.expm(1j*t*A)@P@sl.expm(-1j*t*A)
Vd=(Pc(h)-Pc(-h))/(2*h)
rep("K2a curve velocity equals V",np.linalg.norm(Vd-Vt),1e-7,"constructed curve does not hit the target V")
Pq=Pc(0.7); rep("K2b curve stays rank-k projector",np.linalg.norm(Pq@Pq-Pq)+abs(np.trace(Pq).real-k)+np.linalg.norm(Pq-Pq.conj().T),1e-10,"curve leaves Gr(k)")
# K3 dimension 2k(n-k)
basis=[]
U=np.linalg.eigh(H)[1]
for i in range(k):
    for j in range(k,n):
        for c in (1,1j):
            Y=c*np.outer(U[:,i],U[:,j].conj()); basis.append(np.concatenate([(Y+Y.conj().T).real.ravel(),(Y+Y.conj().T).imag.ravel()]))
rk=np.linalg.matrix_rank(np.array(basis)); rep("K3 real dim = 2k(n-k)",abs(rk-2*k*(n-k)),0.5,f"rank {rk} vs {2*k*(n-k)}")
# K4 positivity: Q(A,A)=||A P||^2 >=0 on complexification; definite on real; NOT definite on complexification
Qf=lambda A,B: np.trace(P@A.conj().T@B)
worst=min((Qf(A_,A_).real for A_ in [tangent(P)+1j*tangent(P) for _ in range(500)]))
rep("K4a Q(A,A) >= 0 on complexification",max(0,-worst),1e-12,"negative value")
mins=min(Qf(V_,V_).real/np.linalg.norm(V_)**2 for V_ in [tangent(P) for _ in range(500)])
rep("K4b definite on real tangent (min Q(V,V)/|V|^2 = 1/2)",abs(mins-0.5),1e-12,"degenerate real direction")
Aup=P@(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))@Qc   # upper block only: = (V+iW)/... complex tangent with A P = 0
rep("K4c complexified form only semi-definite (A!=0, Q(A,A)=0)",abs(Qf(Aup,Aup)),1e-12,"if this were nonzero the 'semi' would be unnecessary")
# K5 block formula and both sign conventions
Xv=P@tangent(P)@Qc; V1=Xv+Xv.conj().T; Yv=P@tangent(P)@Qc+P@(rng.normal(size=(n,n))+0j)@Qc; W1=Yv+Yv.conj().T
q=Qf(V1,W1); txy=np.trace(Xv@Yv.conj().T)
g=0.5*np.trace(V1@W1).real; Om_aug=(1j*np.trace(P@(V1@W1-W1@V1))).real; Om_jun=(-1j*np.trace(P@(V1@W1-W1@V1))).real
rep("K5a Q(V,W)=Tr(XY*)",abs(q-txy),1e-12,"block formula wrong")
rep("K5b Aug convention Q=g-(i/2)Omega, Omega=iTr(P[V,W])",abs(q-(g-0.5j*Om_aug)),1e-12,"")
rep("K5c June convention Q=g+(i/2)Omega, Omega=-iTr(P[V,W])",abs(q-(g+0.5j*Om_jun)),1e-12,"June definition self-inconsistent")
rep("K5d Omega_Aug = -2 Im Tr(XY*)",abs(Om_aug+2*txy.imag),1e-12,"")
# K6 which sign is the standard Berry curvature F=dA, A=i<u|du>, Q_std=g-(i/2)F ?
def Hab(a,b): return H+a*D+b*D2
D2=herm(n)
def u0(a,b): w,U=np.linalg.eigh(Hab(a,b)); return U[:,0]
# small Wilson loop: Berry phase = -Im log prod <u_i|u_{i+1}> = integral A = F*area (counterclockwise)
e=1e-3; pts=[(-e/2,-e/2),(e/2,-e/2),(e/2,e/2),(-e/2,e/2),(-e/2,-e/2)]  # centred plaquette: O(e^2) error
us=[u0(*p) for p in pts]; prod=np.prod([np.vdot(us[i],us[i+1]) for i in range(4)])
F_wilson=-np.angle(prod)/e**2
Pa=lambda a,b: proj(Hab(a,b),1)
dPa=(Pa(h,0)-Pa(-h,0))/(2*h); dPb=(Pa(0,h)-Pa(0,-h))/(2*h); P1=Pa(0,0)
Om_A=(1j*np.trace(P1@(dPa@dPb-dPb@dPa))).real
rep("K6a Aug Omega = standard Berry curvature F_ab (k=1, Wilson loop)",abs(Om_A-F_wilson)/abs(F_wilson),1e-4,f"Omega={Om_A:.6f} F={F_wilson:.6f}")
print("     June sign gives", -Om_A, "i.e. -F: same tensor, opposite orientation of the 2-form")
# k=2 non-abelian: phase of det of overlaps
def U2(a,b): w,U=np.linalg.eigh(Hab(a,b)); return U[:,:2]
Us=[U2(*p) for p in pts]; prod2=np.prod([np.linalg.det(Us[i].conj().T@Us[i+1]) for i in range(4)])
F2=-np.angle(prod2)/e**2
P2=lambda a,b: proj(Hab(a,b),2)
dA=(P2(h,0)-P2(-h,0))/(2*h); dB=(P2(0,h)-P2(0,-h))/(2*h); PP=P2(0,0)
Om2=(1j*np.trace(PP@(dA@dB-dB@dA))).real
rep("K6b k=2: Omega = Tr of non-abelian curvature (det Wilson loop)",abs(Om2-F2)/abs(F2),1e-4,f"Omega={Om2:.6f} F={F2:.6f}")
# K7 parametric pullback matches <d_i u|(1-P)|d_j u> for k=1
du_a=(u0(h,0)*np.exp(-1j*np.angle(np.vdot(u0(0,0),u0(h,0))))-u0(-h,0)*np.exp(-1j*np.angle(np.vdot(u0(0,0),u0(-h,0)))))/(2*h)
du_b=(u0(0,h)*np.exp(-1j*np.angle(np.vdot(u0(0,0),u0(0,h))))-u0(0,-h)*np.exp(-1j*np.angle(np.vdot(u0(0,0),u0(0,-h)))))/(2*h)
Qstd=np.vdot(du_a,(np.eye(n)-P1)@du_b); Qint=np.trace(P1@dPa@dPb)
rep("K7 Q_Pi(d_aPi,d_bPi) = <d_a u|(1-P)|d_b u>",abs(Qstd-Qint),1e-6,f"{Qint:.6f} vs {Qstd:.6f}")
# K8 Remark 2.3: Tr(PVW)+Tr(PWV)=Tr(VW) holds for tangent, fails for generic
G1,G2=herm(n),herm(n)
rep("K8a identity holds on tangent vectors",abs(np.trace(P@V1@W1)+np.trace(P@W1@V1)-np.trace(V1@W1)),1e-12,"")
gen=abs(np.trace(P@G1@G2)+np.trace(P@G2@G1)-np.trace(G1@G2))
rep("K8b identity FAILS for generic operators (kill if residual < 1e-3)",1e-3/gen,1.0,f"generic residual {gen:.3f}")
# K9 Bloch sphere normalisation
sx=np.array([[0,1],[1,0]]);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1,-1])
nv=lambda th: np.array([np.sin(th),0,np.cos(th)])
Pm=lambda th: 0.5*(np.eye(2)-sum(c*s for c,s in zip(nv(th),(sx,sy,sz))))
th=0.6; Vb=(Pm(th+h)-Pm(th-h))/(2*h)
gtt=0.5*np.trace(Vb@Vb).real
rep("K9 g_thth = 1/4 (Fubini-Study)",abs(gtt-0.25),1e-8,f"{gtt}; Tr(VW) would give {2*gtt}")
# K10 Riesz contour
w=np.linalg.eigvalsh(H); 
def riesz(c,r,N=4000):
    z=c+r*np.exp(2j*np.pi*np.arange(N)/N); dz=1j*r*np.exp(2j*np.pi*np.arange(N)/N)*2*np.pi/N
    return sum(np.linalg.inv(zi*np.eye(n)-H)*d for zi,d in zip(z,dz))/(2j*np.pi)
c1=(w[0]+w[1])/2; gapmid=(w[1]+w[2])/2
Pr1=riesz(c1,gapmid-c1); Pr2=riesz(c1-0.3,(gapmid-c1)+0.3*0.99+0.0)
rep("K10a Riesz integral = eigenprojector",np.linalg.norm(Pr1-P),1e-9,"")
rep("K10b contour independence within the gap",np.linalg.norm(Pr2-P),1e-9,"")
Pbad=riesz(c1,(w[2]+w[3])/2-c1)   # contour now encloses 3 levels
rep("K10c opposite: contour past the gap gives a different projector (kill if |diff| < 0.5)",0.5/np.linalg.norm(Pbad-P),1.0,f"|diff|={np.linalg.norm(Pbad-P):.3f}")
# K11 basis independence inside Ran(P)
Uk=np.linalg.eigh(H)[1][:,:k]; Ur=Uk@sl.expm(1j*herm(k)); Prot=Ur@Ur.conj().T
rep("K11 Pi (hence Q) invariant under U(k) rotation of the occupied basis",np.linalg.norm(Prot-P),1e-12,"")
print("\nSUMMARY", sum(r[3] for r in R),"/",len(R),"pass")
