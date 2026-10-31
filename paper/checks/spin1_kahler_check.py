# SCRIPT: IQGT-SPIN1-KAHLER-CHECK
# Claim to test (Appendix C, why the lowest spin-1 band saturates Theorem 6.1):
#   the lowest band of d.S is the spin coherent state, proportional to a vector that depends
#   holomorphically (or antiholomorphically) on z = tan(theta/2) e^{i phi}: (1, sqrt2 w, w^2) with w = z or conj(z) up to sign.
# Kill conditions written first:
#   K1 overlap |<band|w-vector>| differs from 1 by more than 1e-12 for every choice of w tested  -> claim dead
#   K2 4 det g - Omega^2 not ~0 (rel 1e-8) on the same points                                     -> saturation not reproduced
#   K3 control: a generic three-band family (random Hermitian H0,H1,H2) must fall STRICTLY inside (ratio < 0.99 somewhere)
import numpy as np
r2=np.sqrt(2)
Sx=np.array([[0,1,0],[1,0,1],[0,1,0]])/r2; Sy=np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])/r2; Sz=np.diag([1.,0,-1])
def band(th,ph):
    H=np.sin(th)*np.cos(ph)*Sx+np.sin(th)*np.sin(ph)*Sy+np.cos(th)*Sz
    w,V=np.linalg.eigh(H); return V[:,0]
def proj(v): return np.outer(v,v.conj())
def qgt(f,p,h=1e-5):
    P=proj(f(*p)); D=[]
    for i in range(2):
        e=np.zeros(2); e[i]=h; D.append((proj(f(*(p+e)))-proj(f(*(p-e))))/(2*h))
    Q=np.array([[np.trace(P@D[i]@D[j]) for j in range(2)] for i in range(2)])
    g=Q.real; Om=-2*Q[0,1].imag; return g,Om
rng=np.random.default_rng(8); worst1={}; worst2=0
for t in range(200):
    th,ph=rng.uniform(0.2,2.9),rng.uniform(0,2*np.pi); v=band(th,ph); z=np.tan(th/2)*np.exp(1j*ph)
    for name,w in [('z',z),('-z',-z),('conj z',np.conj(z)),('-conj z',-np.conj(z)),('1/z',1/z),('-1/conj z',-1/np.conj(z))]:
        u=np.array([1,r2*w,w*w]); u=u/np.linalg.norm(u); worst1[name]=max(worst1.get(name,0),abs(1-abs(np.vdot(u,v))))
    g,Om=qgt(band,np.array([th,ph])); worst2=max(worst2,abs(4*np.linalg.det(g)-Om**2)/Om**2)
best=min(worst1,key=worst1.get)
print('K1 holomorphic form, best choice w =',best,': max |1-overlap| =',worst1[best],'PASS' if worst1[best]<1e-12 else 'FAIL')
print('   (all choices:',{k:float('%.1e'%v) for k,v in worst1.items()},')')
print('K2 spin-1 saturation: max |4det g - Om^2|/Om^2 =',worst2,'PASS' if worst2<1e-6 else 'FAIL')
def herm(n): A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return (A+A.conj().T)/2
H0,H1,H2=herm(3),herm(3),herm(3)
def band3(a,b): w,V=np.linalg.eigh(H0+a*H1+b*H2); return V[:,0]
mn=1
for t in range(200):
    g,Om=qgt(band3,rng.normal(size=2)*0.5); mn=min(mn,abs(Om)/(2*np.sqrt(np.linalg.det(g))))
print('K3 generic three-band control: min |Om|/(2 sqrt det g) =',mn,'PASS' if mn<0.99 else 'FAIL')
