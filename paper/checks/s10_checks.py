# SCRIPT: S10-CHECKS-V1
# Section 10 (protocol, gap law, heat engine, BEC). Kill conditions written BEFORE the run:
# P1 his reporting inequality (17) RHS equals the sharp RHS 2 sqrt(g_VV g_WW - g_VW^2): kill if rel diff > 1e-12
# P2 (17) holds on 2000 random tangent pairs: kill on any violation
# P3 theta_G unchanged by a frame change U: kill if > 1e-10 deg.  Control: raw readout angle under axis
#    rescaling diag(1,2) must move > 1 deg, else the control cannot discriminate (kill the control)
# G1 gap identity (Q'A'Q')(Q'P) - (Q'P)(PAP) = Q'(A'-A)P: kill if residual > 1e-12
# G2 Davis-Kahan sin-theta: ||Q'P|| <= ||Q'EP|| / sep on 2000 random: kill on any violation
# G3 exact tie: vector angle median > 10 deg while pair projector <= 2x its bound: else kill
# G4 voting scatter at 9949 samples, ratios 4/1.2/1.02, simulated vs first-order sqrt(l1 l2)/((l1-l2)sqrt n):
#    kill if any simulated/predicted outside [0.6, 1.6] (ratio 1.02 is past first order; reported, not graded)
# H1 heat-engine core axis (Untitled21 cell 9 eigenvalues, n=9950): first-order axis scatter, reported
# H2 null band for core-vs-edge axis angle (his cell 10 null, 300 seeds): kill "tail axis is not noise"
#    if the measured 75.43 deg lies inside the null 99.9th percentile
import numpy as np
rng = np.random.default_rng(10)
res = {}
def herm(n):
    a = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n)); return (a+a.conj().T)/2
def proj(A, idx):
    w, U = np.linalg.eigh(A); Us = U[:, idx]; return Us@Us.conj().T, w
def tangent(P):
    n = P.shape[0]; X = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    I = np.eye(n); B = P@X@(I-P); return B + B.conj().T
def g(V,W): return 0.5*np.trace(V@W).real
def Om(P,V,W): return (1j*np.trace(P@(V@W-W@V))).real
# P1/P2
worst1 = 0; viol = 0
for t in range(2000):
    n = rng.integers(3,8); k = rng.integers(1,n)
    P,_ = proj(herm(n), list(range(k))); V = tangent(P); W = tangent(P)
    C1, C2, gvw = g(V,V), g(W,W), g(V,W)
    cos = gvw/np.sqrt(C1*C2); sin = np.sqrt(max(0,1-cos**2))
    r17 = 2*np.sqrt(C1*C2)*sin; rs = 2*np.sqrt(C1*C2-gvw**2)
    worst1 = max(worst1, abs(r17-rs)/rs)
    if abs(Om(P,V,W)) > r17*(1+1e-12): viol += 1
res['P1'] = (worst1, worst1 < 1e-12); res['P2'] = (viol, viol == 0)
# P3
n=6; P,_=proj(herm(n),[0,1]); V=tangent(P); W=tangent(P)
U,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
th = lambda P,V,W: np.degrees(np.arccos(g(V,W)/np.sqrt(g(V,V)*g(W,W))))
d3 = abs(th(P,V,W) - th(U@P@U.conj().T, U@V@U.conj().T, U@W@U.conj().T))
a = np.array([np.cos(0.4), np.sin(0.4)]); b = np.array([np.cos(1.3), np.sin(1.3)]); S = np.diag([1,2])
raw = lambda x,y: np.degrees(np.arccos(abs(x@y)/np.linalg.norm(x)/np.linalg.norm(y)))
dctl = abs(raw(a,b) - raw(S@a,S@b))
res['P3'] = ((d3, dctl), d3 < 1e-10 and dctl > 1)
# G1/G2
worstG1 = 0; violG2 = 0; worstG2 = 0
for t in range(2000):
    n = rng.integers(3,9); k = rng.integers(1,n)
    A = herm(n); E = herm(n)*10**rng.uniform(-3,0); A2 = A+E
    P, wA = proj(A, list(range(k))); P2, wA2 = proj(A2, list(range(k)))
    I = np.eye(n); Q2 = I-P2
    lhs = (Q2@A2@Q2)@(Q2@P) - (Q2@P)@(P@A@P); rhs = Q2@(A2-A)@P
    worstG1 = max(worstG1, np.abs(lhs-rhs).max())
    sep = wA2[k] - wA[k-1]   # spectrum of A on P sits below spectrum of A' on Q'
    if sep <= 0: continue
    s = np.linalg.norm(Q2@P, 2); b = np.linalg.norm(Q2@E@P, 2)/sep
    worstG2 = max(worstG2, s/b)
    if s > b*(1+1e-10): violG2 += 1
res['G1'] = (worstG1, worstG1 < 1e-12); res['G2'] = ((violG2, worstG2), violG2 == 0)
# G3 exact tie
vec, pj, bd = [], [], []
for t in range(500):
    n=6; U,_=np.linalg.qr(rng.normal(size=(n,n))); w=np.array([0,0,3,4,5,6.])
    A=U@np.diag(w)@U.T; E=herm(n).real*1e-3; A2=A+E
    wv,Uv=np.linalg.eigh(A2); v=Uv[:,0]
    u0=U[:,0]; vec.append(np.degrees(np.arccos(min(1,abs(u0@v)))))
    P=U[:,:2]@U[:,:2].T; P2=Uv[:,:2]@Uv[:,:2].T; I=np.eye(n)
    pj.append(np.linalg.norm((I-P2)@P,2)); bd.append(np.linalg.norm((I-P2)@E@P,2)/(wv[2]-0))
res['G3'] = ((np.median(vec), max(np.array(pj)/np.array(bd))), np.median(vec) > 10 and max(np.array(pj)/np.array(bd)) <= 2)
# G4 voting scatter
out = []
for r in (4, 1.2, 1.02):
    l1, l2, N = r, 1.0, 9949; angs = []
    for t in range(400):
        x = rng.normal(size=(N,2))*np.sqrt([l1,l2]); C = np.cov(x.T); w,U = np.linalg.eigh(C)
        angs.append(np.arccos(min(1,abs(U[0,-1]))))
    sim = np.degrees(np.sqrt(np.mean(np.square(angs)))); pred = np.degrees(np.sqrt(l1*l2)/((l1-l2)*np.sqrt(N)))
    out.append((r, round(sim,2), round(pred,2)))
okG4 = all(0.6 <= s/p <= 1.6 for r,s,p in out if r != 1.02)
res['G4'] = (out, okG4)
# H1 heat-engine core axis
l_lo, l_hi, N = 0.00051245, 0.00055124, 9950
res['H1'] = (np.degrees(np.sqrt(l_lo*l_hi)/((l_hi-l_lo)*np.sqrt(N))), True)
# H2 null band, his cell 10 construction
evecs = np.array([[0.22694602,-0.97390734],[-0.97390734,-0.22694602]])
cov_core = evecs@np.diag([l_lo,l_hi])@evecs.T
nulls = []
for s in range(300):
    x = rng.multivariate_normal([0,0], cov_core, size=19899)
    d = np.linalg.norm(x-x.mean(0),axis=1)
    core = x[d <= np.percentile(d,50)]; edge = x[d >= np.percentile(d,90)]
    vc = np.linalg.eigh(np.cov(core.T))[1][:,-1]; ve = np.linalg.eigh(np.cov(edge.T))[1][:,-1]
    nulls.append(np.degrees(np.arccos(min(1,abs(vc@ve)))))
q999 = np.percentile(nulls, 99.9)
res['H2'] = ((round(np.median(nulls),2), round(q999,2), round(max(nulls),2)), 75.43 > q999)
for k,(v,ok) in res.items(): print(k, 'PASS' if ok else 'FAIL', v)
print(sum(ok for v,ok in res.values()), '/', len(res))
