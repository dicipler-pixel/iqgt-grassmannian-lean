# SCRIPT: POLARITON-SATURATION-CHECK
# Question: can the Gianfrate et al. (Nature 578, 381) extraction test the bound 4 det g >= Bz^2 at all?
# Kill (written first): if 4 det g - Bz^2 computed with their Eq. (4) from an ARBITRARY unit pseudospin field is not
# zero to rounding, their maps carry an independent test of the bound. If it is zero, the bound is forced by the method.
import numpy as np
rng=np.random.default_rng(3)
def qgt_from_angles(th,ph,dx):
    tx,ty=np.gradient(th,dx,axis=1),np.gradient(th,dx,axis=0); px,py=np.gradient(ph,dx,axis=1),np.gradient(ph,dx,axis=0)
    gxx=(tx*tx+np.sin(th)**2*px*px)/4; gyy=(ty*ty+np.sin(th)**2*py*py)/4; gxy=(tx*ty+np.sin(th)**2*px*py)/4
    B=0.5*np.sin(th)*(tx*py-ty*px); return gxx,gyy,gxy,B
# (1) arbitrary smooth random unit-vector field (no Hamiltonian at all)
x=np.linspace(-5,5,301); X,Y=np.meshgrid(x,x); dx=x[1]-x[0]
th=1.2+0.6*np.sin(0.7*X+0.3*Y)+0.4*np.cos(0.5*X*Y/3); ph=0.9*X-0.4*Y+0.5*np.sin(Y)
gxx,gyy,gxy,B=qgt_from_angles(th,ph,dx); det=gxx*gyy-gxy**2
r1=np.abs(4*det-B**2).max()/np.abs(B**2).max()
print('arbitrary unit field: max |4det g - Bz^2| / max Bz^2 =',r1, 'FORCED' if r1<1e-12 else 'NOT FORCED')
# (2) their model, Eq. (3), parameters from the paper: alpha=8 ueV, beta=5.92 ueV um^2, Dz=50.2 ueV, phi0 set to 0
a,b,Dz,p0=8.0,5.92,50.2,0.0
k=np.hypot(X,Y); f=np.arctan2(Y,X)
Om=np.stack([a*np.cos(p0)+b*k**2*np.cos(2*f), a*np.sin(p0)-b*k**2*np.sin(2*f), Dz*np.ones_like(k)])
n=Om/np.linalg.norm(Om,axis=0); thm=np.arccos(n[2]); phm=np.arctan2(n[1],n[0])
gxx,gyy,gxy,B=qgt_from_angles(thm,np.unwrap(np.unwrap(phm,axis=1),axis=0),dx); det=gxx*gyy-gxy**2
r2=np.abs(4*det-B**2).max()/np.abs(B**2).max()
print('their model:          max |4det g - Bz^2| / max Bz^2 =',r2)
