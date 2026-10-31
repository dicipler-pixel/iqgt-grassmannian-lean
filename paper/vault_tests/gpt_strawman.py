# SCRIPT: IQGT-GPT-RECOVERY-STRAWMAN-01
# Kill conditions written first:
#  F1 F=(I-P)OP+PO(I-P) equals OP+PO-2POP and [[O,P],P]   kill: residual > 1e-12
#  F2 F=0 iff [O,P]=0 for ARBITRARY (non-Hermitian) O    kill: a case with F=0 but [O,P]!=0, or reverse
#  F3 F is a Grassmannian tangent vector for Hermitian O; F is the g-gradient: g(F,V)=1/2 Tr(O V) for all tangent V   kill: >1e-12
#  F4 flow dP/dt = F(O,P) (Brockett double bracket) ends at the projector onto the TOP-k eigenvectors of O  kill: final distance > 1e-6
#     opposite: dP/dt = -F ends at the BOTTOM-k projector
#  E1 the 240 E8 roots are a spherical 7-design: sum_r (r.u)^m is the same for every unit u, m=2,4,6   kill: spread > 1e-9
#     opposite: m=8 must NOT be constant (spread clearly > 0)
import numpy as np, itertools
rng=np.random.default_rng(8)
def herm(n): A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return (A+A.conj().T)/2
n,k=6,2; H=herm(n); w,U=np.linalg.eigh(H); P=U[:,:k]@U[:,:k].conj().T; I=np.eye(n)
O=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
F=(I-P)@O@P+P@O@(I-P)
print('F1', np.linalg.norm(F-(O@P+P@O-2*P@O@P)), np.linalg.norm(F-((O@P-P@O)@P-P@(O@P-P@O))))
# F2: build O commuting with P (block diagonal) -> F=0 ; generic -> both nonzero
Ob=P@O@P+(I-P)@O@(I-P); Fb=(I-P)@Ob@P+P@Ob@(I-P)
print('F2 commuting: |F|=%.1e |[O,P]|=%.1e ; generic: |F|=%.2f |[O,P]|=%.2f'%(np.linalg.norm(Fb),np.linalg.norm(Ob@P-P@Ob),np.linalg.norm(F),np.linalg.norm(O@P-P@O)))
Oh=herm(n); Fh=(I-P)@Oh@P+P@Oh@(I-P)
print('F3 tangent:',np.linalg.norm(Fh-P@Fh-Fh@P), ' gradient:',max(abs(0.5*np.trace(Fh@V)-0.5*np.trace(Oh@V)) for V in [(lambda X:X+X.conj().T)(P@herm(n)@(I-P)) for _ in range(50)]))
def flow(sign):
    Pt=P.copy(); dt=1e-3
    for _ in range(40000):
        Pt=Pt+sign*dt*((I-Pt)@Oh@Pt+Pt@Oh@(I-Pt))
        ww,UU=np.linalg.eigh(Pt); Pt=UU[:,-k:]@UU[:,-k:].conj().T   # retract to Gr(k)
    return Pt
wo,Uo=np.linalg.eigh(Oh)
print('F4 +F -> top-k:',np.linalg.norm(flow(+1)-Uo[:,-k:]@Uo[:,-k:].conj().T),'   -F -> bottom-k:',np.linalg.norm(flow(-1)-Uo[:,:k]@Uo[:,:k].conj().T))
# E8 roots
R=[]
for i,j in itertools.combinations(range(8),2):
    for si in (1,-1):
        for sj in (1,-1):
            v=np.zeros(8); v[i]=si; v[j]=sj; R.append(v)
for s in itertools.product((0.5,-0.5),repeat=8):
    if sum(x<0 for x in s)%2==0: R.append(np.array(s))
R=np.array(R); print('E8 roots',len(R),'norms',set(np.round((R**2).sum(1),12)))
us=rng.normal(size=(200,8)); us/=np.linalg.norm(us,axis=1)[:,None]
for m in (2,4,6,8):
    vals=((us@R.T)**m).sum(1); print(f'E1 m={m}: mean {vals.mean():.6f} spread {vals.max()-vals.min():.3e}')
