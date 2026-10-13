# SCRIPT: IQGT-ETA-DIPOLE-CHECK
# Claim (Section 7.5): regularising (7.2) by 1/D^2 -> 1/(D^2+eta^2) turns the gap-derivative weight 2 dD/D^3 into
# 2 D dD/(D^2+eta^2)^2, which (i) vanishes at D = 0, (ii) is d/dk of -1/(D^2+eta^2) (so it is the derivative of the regularised
# pair weight, as (7.2) is of 1/D^2), (iii) peaks at D = eta/sqrt3 with value (3 sqrt3/8)/eta^3 per unit |dD| and coupling,
# (iv) returns 2/D^3 as eta -> 0 on an open gap. Kill: any of (i)-(iv) fails numerically.
import numpy as np
ok=True
for eta in (0.05,0.2,1.0):
    D=np.linspace(0,5,200001); w=2*D/(D**2+eta**2)**2
    peak=w.max(); Dp=D[w.argmax()]; pred=3*np.sqrt(3)/8/eta**3
    h=1e-6; der=-(1/((D[1:50]+h)**2+eta**2)-1/((D[1:50]-h)**2+eta**2))/(2*h)
    ok&= w[0]==0 and abs(peak/pred-1)<1e-6 and abs(Dp-eta/np.sqrt(3))<1e-3 and np.allclose(der,w[1:50],rtol=1e-6)
    print(f'eta={eta}: w(0)={w[0]}, peak {peak:.6g} at D={Dp:.5f} (pred {pred:.6g} at {eta/np.sqrt(3):.5f})')
D0=0.7; lim=[2*D0/(D0**2+e**2)**2 for e in (1e-2,1e-4)]; print('eta->0 at D=0.7:',lim,'vs 2/D^3 =',2/D0**3); ok&=abs(lim[-1]/(2/D0**3)-1)<1e-6
print('RESULT','PASS' if ok else 'FAIL')
