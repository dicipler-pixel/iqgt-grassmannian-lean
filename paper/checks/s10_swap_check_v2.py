# SCRIPT: S10-SWAP-CHECK-V2
# Kill (before run): mixing angle of H=[[D/2,g],[g,-D/2]] must satisfy tan 2theta = 2g/D to 1e-12;
# a crossing narrower than the sampling step must show one step > 80 deg in 29 samples; a wide one none > 15 deg.
# Edge/core e-h balance shift (Untitled21 cell 33, ~995 edge shots each): binomial z for e share 0.623 -> 0.387.
import numpy as np
err=0
for t in range(200):
    D,g=np.random.uniform(-3,3),np.random.uniform(0.01,2)
    w,U=np.linalg.eigh([[D/2,g],[g,-D/2]]); v=U[:,1]
    v=v*np.sign(v[0]) if v[0]!=0 else v; th=np.arctan2(v[1],v[0])
    err=max(err,abs(np.sin(2*th)*D-2*g*np.cos(2*th)))
print('tan2theta err',err, 'PASS' if err<1e-9 else 'FAIL')
def maxstep(g):
    x=np.linspace(-1,1,29)+0.5*(2/28)
    vs=[np.linalg.eigh([[a,g],[g,-a]])[1][:,1] for a in x]
    return max(np.degrees(np.arccos(min(1,abs(vs[i]@vs[i+1])))) for i in range(28))
n,w=maxstep(0.002),maxstep(0.4)
print('narrow max step %.1f  wide %.1f'%(n,w),'PASS' if n>80 and w<15 else 'FAIL')
p1,p2,N=0.623,0.387,995; se=np.sqrt(p1*(1-p1)/N+p2*(1-p2)/N); print('edge e-share shift z = %.1f'%((p1-p2)/se))
