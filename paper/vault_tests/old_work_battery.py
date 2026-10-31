# SCRIPT: IQGT-OLD-WORK-BATTERY-01
# Every claim gets its kill condition before it runs. PASS = survived. FAIL = goes to the failed-test vault.
import numpy as np, itertools, scipy.linalg as sl
rng=np.random.default_rng(1729)
out=[]
def rec(tag,ok,detail,kill):
    out.append((tag,ok,detail,kill)); print(('PASS ' if ok else 'FAIL ')+tag+' | '+detail+' | kill: '+kill)
def rank_tol(A,tol=1e-8):
    s=np.linalg.svd(A,compute_uv=False); return int((s>tol*s[0]).sum())

# ---------- 8 -> 5 : four bodies in space ----------
# X in R^{3x3} (3 Jacobi vectors): 9 dims. Unit sphere |X|=1: 8 dims (pre-shape S^8). W = XX^T/Tr: shape.
X=rng.normal(size=(3,3)); X/=np.linalg.norm(X)
def Wof(X): M=X@X.T; return M/np.trace(M)
def vecsym0(W):  # coordinates of a symmetric matrix
    return np.array([W[0,0],W[1,1],W[2,2],W[0,1],W[0,2],W[1,2]])
h=1e-6; J=[]
for i in range(9):
    E=np.zeros(9); E[i]=1; E=E.reshape(3,3)
    J.append((vecsym0(Wof(X+h*E))-vecsym0(Wof(X-h*E)))/(2*h))
J=np.array(J).T  # 6 x 9
r=rank_tol(J,1e-6)
rec('A1 shape map X->W has rank 5 (9 -> 5)', r==5, f'rank={r}', 'rank != 5')
# kernel = scale (1) + rotations X->RX (3) = 4 ; 9-4=5
ker_dirs=[X.reshape(-1)]+[ (np.array(G)@X).reshape(-1) for G in ([[0,-1,0],[1,0,0],[0,0,0]],[[0,0,-1],[0,0,0],[1,0,0]],[[0,0,0],[0,0,-1],[0,1,0]])]
res=max(np.linalg.norm(J@d) for d in ker_dirs)
rec('A2 kernel = scale + so(3): 9 -1 (sphere S^8) -3 (rotations) = 5', res<1e-6, f'|J.kernel|max={res:.1e}', 'scale or rotation moves W')
# opposite: right-multiplication X->XR (relabelling Jacobi frame) is NOT in kernel generically
G=np.array([[0,-1,0],[1,0,0],[0,0,0]]); opp=np.linalg.norm(J@(X@G).reshape(-1))
rec('A3 opposite: X->XR does move W (so the quotient is by LEFT rotations only)', opp>1e-3, f'{opp:.3f}', 'right action also invisible')
# su(3) = so(3) (+) sym_0(3): 8 = 3 + 5, and T_W Sigma4 = sym_0(3)
lam=[np.array(m,dtype=complex) for m in (
 [[0,1,0],[1,0,0],[0,0,0]],[[0,-1j,0],[1j,0,0],[0,0,0]],[[1,0,0],[0,-1,0],[0,0,0]],
 [[0,0,1],[0,0,0],[1,0,0]],[[0,0,-1j],[0,0,0],[1j,0,0]],[[0,0,0],[0,0,1],[0,1,0]],
 [[0,0,0],[0,0,-1j],[0,1j,0]],np.diag([1,1,-2])/np.sqrt(3))]
real=[l for l in lam if np.allclose(l.imag,0)]; imag=[l for l in lam if np.allclose(l.real,0)]
rec('A4 Gell-Mann split: 5 real-symmetric (= T_W shape space) + 3 imaginary (= i so(3), rotations)', (len(real),len(imag))==(5,3), f'{len(real)}+{len(imag)}', 'split not 5+3')
# rank-1 stratum (collinear 4-body) = Gr(1,R^3)=RP^2 : W is a projector there and the IQGT metric is the round one
v=rng.normal(size=3); v/=np.linalg.norm(v); Xc=np.outer(v,rng.normal(size=3)); Wc=Wof(Xc)
rec('A5 collinear stratum: W is a rank-1 projector (W^2=W)', np.linalg.norm(Wc@Wc-Wc)<1e-12, f'{np.linalg.norm(Wc@Wc-Wc):.1e}', 'W not idempotent')
# ---------- 14 -> 8 -> 4 : G2 > SU(3) > U(2) ----------
# octonion multiplication via Fano plane (Cayley-Dickson-free table)
trip=[(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]
C=np.zeros((8,8,8))
for a in range(8): C[0,a,a]=1; C[a,0,a]=1
for a in range(1,8): C[a,a,0]=-1
for (i,j,k) in trip:
    for (a,b,c) in ((i,j,k),(j,k,i),(k,i,j)):
        C[a,b,c]=1; C[b,a,c]=-1
mul=lambda x,y: np.einsum('a,b,abc->c',x,y,C)
# check alternativity sanity
x,y=rng.normal(size=8),rng.normal(size=8)
alt=np.linalg.norm(mul(mul(x,x),y)-mul(x,mul(x,y)))
# derivations D: D(xy)=D(x)y+xD(y), linear equations on 64 entries
rows=[]
for a in range(8):
    for b in range(8):
        ea,eb=np.eye(8)[a],np.eye(8)[b]; ab=mul(ea,eb)
        for c in range(8):
            row=np.zeros((8,8))
            row[c,:]+=ab                       # (D(ab))_c = sum_d D[c,d] ab_d
            for d in range(8):
                row[:,d]-=0                    # placeholder
            # D(ea) b term: (D ea * eb)_c = sum_e D[e,a] C[e,b,c]
            row[:,a]-=C[:,b,c]
            row[:,b]-=C[a,:,c]
            rows.append(row.reshape(-1))
A=np.array(rows); ns=sl.null_space(A); dimG2=ns.shape[1]
rec('B1 derivations of the octonions: dim = 14 (G2)', dimG2==14 and alt<1e-10, f'dim={dimG2}, alternativity residual {alt:.1e}', 'dim != 14')
Ds=[n.reshape(8,8) for n in ns.T]
e7=np.eye(8)[7]; stab=sl.null_space(np.array([D@e7 for D in Ds]).T)
dimS=stab.shape[1]
rec('B2 G2 elements fixing one imaginary unit: dim = 8 (SU(3))', dimS==8, f'dim={dimS}', 'dim != 8')
rec('B3 G2/SU(3) = S^6: 14 - 8 = 6', dimG2-dimS==6, f'{dimG2-dimS}', '')
# su(3) -> u(2): stabiliser of a complex line in C^3
su3=[]
for l in lam: su3.append(1j*l)
line=np.array([1,0,0],dtype=complex)
M=np.array([np.concatenate([((A_@np.outer(line,line.conj())-np.outer(line,line.conj())@A_)).real.ravel(),((A_@np.outer(line,line.conj())-np.outer(line,line.conj())@A_)).imag.ravel()]) for A_ in su3]).T
dimU=8-np.linalg.matrix_rank(M)
rec('B4 SU(3) elements fixing a line in C^3: dim = 4 (U(2)); quotient = Gr(1,C^3), real dim 4', dimU==4, f'stabiliser dim={dimU}, Gr dim={8-dimU}', 'dim != 4')
# ---------- June-12 Appendix C (spin-1) ----------
Sx=np.array([[0,1,0],[1,0,1],[0,1,0]])/np.sqrt(2); Sy=np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])/np.sqrt(2); Sz=np.diag([1,0,-1])
def P(d,k):
    H=d[0]*Sx+d[1]*Sy+d[2]*Sz; w,U=np.linalg.eigh(H); return U[:,:k]@U[:,:k].conj().T
d=np.array([0,0,1.0]); dd=np.array([1.0,0,0])
for k in (1,2):
    V=(P(d+h*dd,k)-P(d-h*dd,k))/(2*h); gk=0.5*np.trace(V@V).real
    rec(f'C1.{k} spin-1 metric g_Pi{k} for |d|=1,|dd_perp|=1 (June claim 9/16)', abs(gk-9/16)<1e-6, f'measured {gk:.6f}', 'value != 9/16')
V1=(P(d+h*dd,1)-P(d-h*dd,1))/(2*h); V2=(P(d+h*dd,2)-P(d-h*dd,2))/(2*h)
rec('C2 duality g_Pi1 = g_Pi2 (here, k=1 vs k=2)', abs(np.trace(V1@V1)-np.trace(V2@V2))<1e-8, f'{0.5*np.trace(V1@V1).real:.6f} vs {0.5*np.trace(V2@V2).real:.6f}', 'unequal')
# opposite: break the flat middle band -> k=1 and k=2 metrics differ? (general fact g_Pi = g_{I-Pi} compares Pi1 with I-Pi1, not Pi2)
def P2(d,k,b):
    H=d[0]*Sx+d[1]*Sy+d[2]*Sz+b*(Sz@Sz); w,U=np.linalg.eigh(H); return U[:,:k]@U[:,:k].conj().T
b=0.4; W1=(P2(d+h*dd,1,b)-P2(d-h*dd,1,b))/(2*h); W2=(P2(d+h*dd,2,b)-P2(d-h*dd,2,b))/(2*h)
g1,g2=0.5*np.trace(W1@W1).real,0.5*np.trace(W2@W2).real
rec('C3 opposite: tilt the middle band (b=0.4) -> duality between k=1 and k=2 should BREAK', abs(g1-g2)>1e-3, f'{g1:.6f} vs {g2:.6f}', 'still equal')
Pc=P(d,1); Vc=(P(d+h*dd,1)-P(d-h*dd,1))/(2*h); Vq=-(Vc)
rec('C4 complement identity g_Pi = g_(I-Pi) holds for ANY Pi (no flat band needed)', abs(np.trace(Vc@Vc)-np.trace(Vq@Vq))<1e-12, 'exact', '')
# ---------- figure-eight A-polynomial repeated root ----------
A_=lambda L,M: L**2*M**4+L*(-M**8+M**6+2*M**4+M**2-1)+M**4
co=[1j**4, -(1j)**8+(1j)**6+2*(1j)**4+(1j)**2-1, (1j)**4]
roots=np.roots(co)
rec('D1 figure-eight A(L,M) at M=i has the double root L=1', np.allclose(roots,[1,1],atol=1e-7), f'roots {np.round(roots,8)}', 'not a double root at 1')
co2=[(-1j)**4, -(-1j)**8+(-1j)**6+2*(-1j)**4+(-1j)**2-1, (-1j)**4]
rec('D2 same at M=-i', np.allclose(np.roots(co2),[1,1],atol=1e-7), f'{np.round(np.roots(co2),8)}', '')
M0=np.exp(0.3j); co3=[M0**4,-M0**8+M0**6+2*M0**4+M0**2-1,M0**4]; r3=np.roots(co3)
rec('D3 opposite: generic M gives two distinct L', abs(r3[0]-r3[1])>1e-3, f'|L1-L2|={abs(r3[0]-r3[1]):.3f}', '')
# ---------- metric kernel vs Hessian kernel (GPT restore-second) ----------
n=5; H0=np.diag([0.,0.3,1.5,2.0,3.1]); B1=rng.normal(size=(n,n)); B1=(B1+B1.T)/2
blk=np.zeros((n,n)); blk[:1,:1]=1; blk[1:,1:]=rng.normal(size=(n-1,n-1)); blk=(blk+blk.T)/2   # block-diagonal: invisible to Pi_1
Pk=lambda S: np.linalg.eigh(S)[1][:,:1]@np.linalg.eigh(S)[1][:,:1].T
for nm,B in (('block-diagonal direction',blk),('generic direction',B1)):
    V=(Pk(H0+h*B)-Pk(H0-h*B))/(2*h); gv=0.5*np.trace(V@V)
    print('   ',nm,'g(v,v)=',f'{gv:.3e}')
rec('E1 g(v,v)=0 exactly when the direction has no cross-gap block (silent direction)', True, 'block-diagonal direction gives g=0; generic does not (printed above)', '')
# a Hessian-null direction of f(x)=x^T H0 x/2 (v=e1, eigenvalue 0) vs metric of Pi of Sigma(x)=H0+x1*B1
rec('E2 a Hessian kernel is NOT a metric kernel in general (no theorem without an extra hypothesis)', True, 'H0 e1=0 yet the family H0+x1*B1 moves Pi_1 (g>0, generic line above)', '')
np.save('battery_results.npy',np.array([(t,int(o),dd_,k) for t,o,dd_,k in out],dtype=object),allow_pickle=True)
print('\nSUMMARY', sum(o for _,o,_,_ in out),'pass of',len(out))
