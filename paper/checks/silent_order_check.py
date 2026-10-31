# SCRIPT: IQGT-SILENT-ORDER
# Claim (extends Proposition 8.6): along a direction where the projector's derivatives vanish through order m-1,
# (i) the m-th derivative is block off-diagonal (a tangent vector);
# (ii) d(eps) = |eps|^m/m! * sqrt(g(P^(m),P^(m))) + O(eps^(m+1));
# (iii) m even: Pi(+eps) and Pi(-eps) agree to order eps^m (two sheets meet: a genuine fold);
#       m odd: they separate at order eps^m (the curve passes through; t = eps^m removes the silence).
# Kill: a diagonal block at order m, the ratio in (ii) not -> 1, or the even/odd separation test failing either way.
import numpy as np, math
sx=np.array([[0,1],[1,0]]);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1])
def P(th,ph): w,U=np.linalg.eigh(np.cos(ph)*np.sin(th)*sx+np.sin(ph)*np.sin(th)*sy+np.cos(th)*sz); return np.outer(U[:,0],U[:,0].conj())
th0=1.4; ok=True
for m in (2,3,4,5):
    f=lambda e: e**m/math.factorial(m)          # phi(eps) = eps^m/m!: silent through order m-1
    Pe=lambda e: P(th0,f(e)); P0=Pe(0.0)
    # m-th derivative of Pi along eps at 0 equals dPi/dphi at phi=0 (chain rule, since phi^(m)(0)=1)
    h=1e-6; Pm=(P(th0,h)-P(th0,-h))/(2*h)
    diag=np.linalg.norm(P0@Pm@P0)+np.linalg.norm((np.eye(2)-P0)@Pm@(np.eye(2)-P0))
    gm=0.5*np.trace(Pm@Pm).real
    rat=[ (np.linalg.norm(Pe(e)-P0)/np.sqrt(2))/(e**m/math.factorial(m)*np.sqrt(gm)) for e in (0.2,0.05) ]
    sep=[ np.linalg.norm(Pe(e)-Pe(-e))/np.linalg.norm(Pe(e)-P0) for e in (0.2,0.05) ]
    even=(m%2==0)
    sep_ok = (sep[1]<1e-12) if even else (abs(sep[1]-2)<1e-3)
    good= diag<1e-8 and abs(rat[1]-1)<1e-3 and sep_ok; ok&=good
    print(f'm={m}: diag blocks {diag:.1e}; d/(eps^m/m! sqrt g) {rat[0]:.6f}, {rat[1]:.6f}; |Pi(+e)-Pi(-e)|/|Pi(e)-Pi(0)| {sep[1]:.3e} ->',
          'two sheets meet (fold)' if sep[1]<1e-6 else 'passes through (removable)', 'OK' if good else 'FAIL')
print('RESULT','PASS' if ok else 'FAIL')
