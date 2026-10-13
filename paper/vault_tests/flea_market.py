# SCRIPT: FLEA-MARKET-01  (play with the old chat's items; each check states what would sink it)
import numpy as np, cmath
# A. Kuramoto locking vs two-level mixing. Kill: the locking threshold |dw|=2K does not land on a fixed mixing angle.
def locks(dw,K,T=400,dt=0.01):
    D=0.3
    for _ in range(int(T/dt)): D+=dt*(dw-2*K*np.sin(D))
    v=dw-2*K*np.sin(D); return abs(v)<1e-6
K=1.0
for dw in [1.0,1.9,1.99,2.01,2.5]:
    th=0.5*np.degrees(np.arctan2(2*K,dw))
    print(f'A  dw/2K={dw/2:.3f}  locks={locks(dw,K)}  two-level mixing angle={th:.2f} deg')
# B. The old knot expression JP = t^-3/2 + t^3/2. Kill for "it is a Jones polynomial": V(1) must be (-2)^(c-1) and V(e^{2pi i/3}) must be 1.
f=lambda t,s: s*(t**-1.5+t**1.5)
for s in (1,-1):
    print('B  sign',s,' V(1)=',f(1.0,s),'  V(w)=',np.round(f(cmath.exp(2j*np.pi/3),s),6))
# C. 3-6-9: doubling mod 9 as a permutation of the nine residues; its graph-subspace angles (paper Sec 12.4).
perm=[(2*r)%9 for r in range(9)]; n=9
G=np.zeros((2*n,n)); I=np.zeros((2*n,n))
for k in range(n): G[k,k]=1; G[n+perm[k],k]=1; I[k,k]=1; I[n+k,k]=1
s=np.linalg.svd((I/np.sqrt(2)).T@(G/np.sqrt(2)),compute_uv=False); th=np.arccos(np.clip(s,-1,1)); th=np.sort(th[th>1e-9])[::-1]
moved=sum(perm[k]!=k for k in range(n))
print('C  perm',perm,' cycles: (1 2 4 8 7 5)(3 6)(0)  moved',moved)
print('C  angles deg',np.round(np.degrees(th),2),' sum sin^2 =',round((np.sin(th)**2).sum(),4),' = moved/2 ->',moved/2,'  length sqrt(sum th^2) =',round(np.sqrt((th**2).sum()),4))
# compare with the same count rearranged: four swaps moving 8 points
G2=np.zeros((2*n,n)); p2=[0,2,1,4,3,6,5,8,7]
for k in range(n): G2[k,k]=1; G2[n+p2[k],k]=1
s2=np.linalg.svd((I/np.sqrt(2)).T@(G2/np.sqrt(2)),compute_uv=False); t2=np.arccos(np.clip(s2,-1,1))
print('C  four swaps, also 8 moved: sum sin^2 =',round((np.sin(t2)**2).sum(),4),' length =',round(np.sqrt((t2**2).sum()),4))
# D. Memory without a changed present state: Berry holonomy on a loop that returns the projector to where it started.
def lower(theta,phi):
    return np.array([np.sin(theta/2)*np.exp(-1j*phi)*-1,np.cos(theta/2)])  # any smooth gauge
th0=np.radians(60); N=4000; ph=np.linspace(0,2*np.pi,N+1)
vs=[np.array([np.cos(th0/2),np.sin(th0/2)*np.exp(1j*p)]) for p in ph]
prod=np.prod([np.vdot(vs[i],vs[i+1]) for i in range(N)])
print('D  loop at 60 deg: projector returns exactly; Berry phase =',round(-np.angle(prod)%(2*np.pi),4),' half solid angle =',round(np.pi*(1-np.cos(th0)),4))
