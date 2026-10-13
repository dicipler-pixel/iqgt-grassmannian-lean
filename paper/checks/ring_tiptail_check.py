# SCRIPT: IQGT-RING-TIP-TAIL
# A worked example (reproduction level): light going both ways round a ring.
# Model: two longitudinal orders, each a forward/backward pair; Sagnac splitting Om between the directions,
# backscatter coupling kap within each pair, free spectral range F between orders. Perturbation E (random, norm eps).
# Claims and kill conditions, fixed before running:
#  (a) tail: the pair subspace of the lower order turns by at most ||E||/F_eff (Davis-Kahan), at every Om. Kill: any violation.
#  (b) tip: the single lowest mode turns ~ ||E||/delta, delta = sqrt(Om^2 + 4 kap^2); log-log slope vs delta in [-1.15,-0.85]. Kill otherwise.
#  (c) mixing angle (1/2) arctan(2 kap/Om) = 22.5 deg exactly at Om = 2 kap: standing waves below, travelling waves above. Kill if not.
#  (d) active laser (Adler): psi' = Om - OmL sin psi locks iff |Om| <= OmL. Kill: 0.995 OmL fails to lock or 1.005 OmL locks.
import numpy as np
rng=np.random.default_rng(5)
def H(Om,kap,F=10.0):
    h=np.zeros((4,4))
    for o,base in ((0,0.0),(2,F)):
        h[o,o]=base+Om/2; h[o+1,o+1]=base-Om/2; h[o,o+1]=h[o+1,o]=kap
    return h
def lowproj(h,k): w,U=np.linalg.eigh(h); return U[:,:k]
def angle(A,B):
    s=np.linalg.svd(A.T@B,compute_uv=False); return np.degrees(np.arccos(np.clip(s.min(),-1,1)))
kap=0.05; eps=0.002; F=10.0; ok=True; rows=[]
for Om in np.geomspace(0.01,3.0,12):
    h=H(Om,kap,F); tip=[];tail=[]
    for _ in range(60):
        E=rng.normal(size=(4,4)); E=(E+E.T)/2; E*=eps/np.linalg.norm(E,2)
        tip.append(angle(lowproj(h,1),lowproj(h+E,1))); tail.append(angle(lowproj(h,2),lowproj(h+E,2)))
    delta=np.sqrt(Om**2+4*kap**2); w=np.linalg.eigvalsh(h); Fe=w[2]-w[1]
    bound=np.degrees(np.arcsin(min(1,eps/Fe)))
    ok&= max(tail)<=bound*1.0001
    rows.append((Om,delta,np.median(tip),np.median(tail),bound))
    print(f'Om={Om:7.4f}  delta={delta:.4f}  tip median {np.median(tip):8.4f} deg  tail median {np.median(tail):.5f} deg (bound {bound:.5f})')
r=np.array(rows); sl=np.polyfit(np.log(r[:,1]),np.log(r[:,2]),1)[0]; print('tip slope vs delta',round(sl,3)); ok&= -1.15<=sl<=-0.85
th=0.5*np.degrees(np.arctan2(2*kap,2*kap)); print('mixing angle at Om=2kap',th); ok&= abs(th-22.5)<1e-12
def locks(Om,OmL=1.0,T=400,dt=0.01):
    p=0.3
    for _ in range(int(T/dt)): p+=dt*(Om-OmL*np.sin(p))
    return abs(Om-OmL*np.sin(p))<1e-6
l1,l2=locks(0.995),locks(1.005); print('Adler: 0.995 locks',l1,' 1.005 locks',l2); ok&= l1 and not l2
np.save('ring_rows.npy',r)
print('RESULT','PASS' if ok else 'FAIL')
