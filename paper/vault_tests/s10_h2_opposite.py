# SCRIPT: S10-H2-OPPOSITE-V1
# Opposite of H2: test the TAIL'S OWN SHAPE, not its angle to the ill-defined core axis.
# Kill (written before run): "the g-cloud tail is anisotropic beyond noise" dies if the measured
#   edge eigenvalue ratio 2.2163 (Untitled21 cell 9) lies inside the null 99.9th percentile.
# Second: the tail axis against the TRUE (not estimated) noise axis; kill if 75.43-ish measured edge
#   axis angle to the true core axis (computed from printed eigenvectors) is inside the null 99.9th pct.
import numpy as np
rng=np.random.default_rng(11)
l_lo,l_hi=0.00051245,0.00055124
ev=np.array([[0.22694602,-0.97390734],[-0.97390734,-0.22694602]])
C=ev@np.diag([l_lo,l_hi])@ev.T; vtrue=ev[:,1]
ve_meas=np.array([-0.46464346,0.88549786])
ang_meas=np.degrees(np.arccos(abs(vtrue@ve_meas)))
ratios=[];angs=[]
for s in range(300):
    x=rng.multivariate_normal([0,0],C,size=19899); d=np.linalg.norm(x-x.mean(0),axis=1)
    e=x[d>=np.percentile(d,90)]; w,U=np.linalg.eigh(np.cov(e.T))
    ratios.append(w[1]/w[0]); angs.append(np.degrees(np.arccos(min(1,abs(vtrue@U[:,1])))))
q=np.percentile(ratios,99.9); qa=np.percentile(angs,99.9)
print('edge ratio measured 2.2163 | null median %.3f  99.9pct %.3f  max %.3f ->'%(np.median(ratios),q,max(ratios)),'LIVES' if 2.2163>q else 'KILLED')
print('edge axis vs true noise axis measured %.2f deg | null median %.2f 99.9pct %.2f ->'%(ang_meas,np.median(angs),qa),'LIVES' if ang_meas>qa else 'KILLED')
