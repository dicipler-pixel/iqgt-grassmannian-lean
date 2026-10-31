# SCRIPT: IQGT-S12-CLOSING-CHECKS
# Kill conditions written before the run:
# W1 from the saved BEC numbers (odd/even + 2 random splits): every cut with relative gap >= 0.10 has
#    measured/predicted in [0.93, 1.12] (claimed 18 of 18). Kill if any falls outside or the count is not 18.
# W2 at the cut that closes each pair (k=3, 8, 11) the subspace angle is below the single-vector angle in all
#    three splits; at the cut inside each pair (k=2, 7, 10) the two are within 1 deg (a rank-k projector whose
#    last vector is free moves with it). Kill otherwise.
# W3 corridor: three-level cluster inside an outer spectrum; for random perturbations E of fixed size, the
#    cluster projector moves at most ||E||/Delta_out (Davis-Kahan) at every inner gap, while the single vector
#    exceeds that bound once the inner gap falls below Delta_out. Kill if the tail bound is ever violated or the
#    tip never exceeds it.
# W4 speed ceiling: along a driven path the metric speed sqrt(g) never exceeds sqrt(k)||dSigma/dt||/Delta (Thm 5.2);
#    a two-level off-diagonal drive reaches it (ratio > 0.99). Kill otherwise.
import json, numpy as np
d=json.load(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','data','bec_numbers.json'))); res={}
cuts=[x for r in d['T1'] if r['label']!='July vs August' for x in r['rows']]
gapped=[x for x in cuts if x['rel_gap']>=0.10]
rs=[x['ratio'] for x in gapped]
res['W1']=((len(gapped),round(min(rs),3),round(max(rs),3)),len(gapped)==18 and min(rs)>=0.93 and max(rs)<=1.12)
ok=True; rows=[]
for r in d['T1'][:3]:
    R={x['k']:x for x in r['rows']}
    for c in (3,8,11):
        ok&= R[c]['projector_deg']<R[c]['vector_deg']; rows.append((r['label'][:6],c,round(R[c]['vector_deg'],1),round(R[c]['projector_deg'],1)))
    for c in (2,7,10):
        ok&= abs(R[c]['projector_deg']-R[c]['vector_deg'])<1.0
res['W2']=(rows,ok)
rng=np.random.default_rng(4); viol=0; tip_exceeds=False; tab=[]
for din in (1.0,0.3,0.1,0.03,0.01):
    w=np.array([0,din,2*din,3.0,3.5,4.0,5.0]); n=7; Dout=3.0-2*din
    tips=[];tails=[];bnd=[]
    for t in range(300):
        U,_=np.linalg.qr(rng.normal(size=(n,n))); A=U@np.diag(w)@U.T
        E=rng.normal(size=(n,n)); E=(E+E.T)/2; E*=0.005/np.linalg.norm(E,2)
        w2,U2=np.linalg.eigh(A+E)
        P=U[:,:3]@U[:,:3].T; P2=U2[:,:3]@U2[:,:3].T
        tail=np.linalg.norm((np.eye(n)-P2)@P,2); b=np.linalg.norm((np.eye(n)-P2)@E@P,2)/(w2[3]-w[2])
        if tail>b*(1+1e-9): viol+=1
        tip=np.sqrt(max(0,1-(U[:,1]@U2[:,1])**2))
        tips.append(np.degrees(np.arcsin(min(1,tip)))); tails.append(np.degrees(np.arcsin(tail))); bnd.append(np.degrees(np.arcsin(min(1,0.005/Dout))))
    tab.append((din,round(np.median(tips),3),round(np.median(tails),4),round(bnd[0],4)))
    if np.median(tips)>bnd[0]: tip_exceeds=True
res['W3']=((viol,tab),viol==0 and tip_exceeds)
worst=0
for t in range(200):
    n=6; k=2; H0=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); H0=(H0+H0.conj().T)/2; H1=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); H1=(H1+H1.conj().T)/2
    s=rng.uniform(-1,1); f=lambda s: (lambda w,U:(U[:,:k]@U[:,:k].conj().T,w))(*np.linalg.eigh(H0+s*H1))
    (P,w)=f(s); V=(f(s+1e-6)[0]-f(s-1e-6)[0])/2e-6; g=0.5*np.trace(V@V).real
    ceil=np.sqrt(k)*np.linalg.norm(H1,2)/(w[k]-w[k-1]); worst=max(worst,np.sqrt(g)/ceil)
sx=np.array([[0,1],[1,0]]); sz=np.diag([1.,-1.]); r2=[]
for D in (0.5,1,2):
    f=lambda s:(lambda w,U:U[:,:1]@U[:,:1].T)(*np.linalg.eigh(D/2*sz+s*sx)); V=(f(1e-7)-f(-1e-7))/2e-7
    r2.append(np.sqrt(0.5*np.trace(V@V))/(1*np.linalg.norm(sx,2)/D))
res['W4']=((round(worst,3),[round(x,4) for x in r2]),worst<=1 and min(r2)>0.99)
for k,(v,o) in res.items(): print(k,'PASS' if o else 'FAIL',v)
