# SCRIPT: PROP35-HS-REPAIR
# Kill: any case with ||V||_HS > sqrt(2k)||dS||/Delta, or with ||X||_HS > ||Pi dS (I-Pi)||_HS/Delta, using V from finite differences.
import numpy as np
rng=np.random.default_rng(1); worst1=worst2=0
for t in range(400):
    n=rng.integers(4,10); k=rng.integers(1,n)
    w=np.sort(rng.normal(size=n)*3); Q=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))[0]
    I=rng.choice(n,k,replace=False)  # non-adjacent allowed
    H=Q@np.diag(w)@Q.conj().T; A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); dS=(A+A.conj().T)/2
    P=lambda M: (lambda e: (e[1][:,I]@e[1][:,I].conj().T))(np.linalg.eigh(M))
    h=1e-6; V=(P(H+h*dS)-P(H-h*dS))/(2*h)
    Pi=P(H); Ic=[j for j in range(n) if j not in I]; Delta=min(abs(w[a]-w[b]) for a in I for b in Ic)
    X=Pi@V@(np.eye(n)-Pi)
    worst1=max(worst1,np.linalg.norm(V)/(np.sqrt(2*k)*np.linalg.norm(dS,2)/Delta))
    worst2=max(worst2,np.linalg.norm(X)/(np.linalg.norm(Pi@dS@(np.eye(n)-Pi))/Delta))
print('ratio to stated bound', worst1, 'PASS' if worst1<=1+1e-6 else 'FAIL')
print('ratio to HS step', worst2, 'PASS' if worst2<=1+1e-6 else 'FAIL')
