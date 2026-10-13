# SCRIPT: IQGT-REBUILD-S8-CHECKS-01
# Section 8 (transport across degenerate strata). Two-band family d(x,y)=(sin y cos kx, sin y sin kx, cos y), k=3, |d|=1 (gap 2 everywhere).
import numpy as np
from scipy.integrate import solve_ivp
R=[]
def rep(tag,val,thr,note=''):
    ok=val<thr; R.append(ok); print(('PASS ' if ok else 'FAIL ')+f'{tag}: {val:.3e} (kill if >= {thr:g}) {note}')
K=3
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.diag([1.,-1.]).astype(complex)
def Pm(x,y):
    d=np.array([np.sin(y)*np.cos(K*x),np.sin(y)*np.sin(K*x),np.cos(y)])
    return 0.5*(np.eye(2)-(d[0]*sx+d[1]*sy+d[2]*sz))
def gnum(x,y,h=1e-5):
    Vx=(Pm(x+h,y)-Pm(x-h,y))/(2*h); Vy=(Pm(x,y+h)-Pm(x,y-h))/(2*h)
    return np.array([[0.5*np.trace(Vx@Vx).real,0.5*np.trace(Vx@Vy).real],[0.5*np.trace(Vy@Vx).real,0.5*np.trace(Vy@Vy).real]])
pts=[(0.3,0.4),(1.1,0.9),(2.0,1.3)]
err=max(np.abs(gnum(x,y)-np.diag([K*K*np.sin(y)**2/4,0.25])).max() for x,y in pts)
rep('R0 numerical metric from Pi equals diag(k^2 sin^2 y /4, 1/4)',err,1e-8)
# principal angle reading (Prop 8.3): sqrt(g_par) dx equals the principal angle between Pi(x,y) and Pi(x+dx,y)
x,y,dx=0.4,0.7,1e-4
w1=np.linalg.eigh(Pm(x,y))[1][:,-1]; w2=np.linalg.eigh(Pm(x+dx,y))[1][:,-1]
ang=np.arccos(min(1,abs(np.vdot(w1,w2))))
rep('R1 principal angle between neighbouring projectors = sqrt(g_par) dx (eigenvector overlaps, no trace identity)',abs(ang-np.sqrt(K*K*np.sin(y)**2/4)*dx)/ang,1e-6)
# geodesics of g = gpar(y) dx^2 + gperp dy^2, no conservation law imposed
gpar=lambda y: K*K*np.sin(y)**2/4; dgpar=lambda y: K*K*np.sin(2*y)/4; gperp=0.25
def rhs(t,s):
    x,y,vx,vy=s
    ax=-(dgpar(y)/gpar(y))*vx*vy
    ay=(dgpar(y)/(2*gperp))*vx*vx
    return [vx,vy,ax,ay]
def launch(y0,theta):
    # unit speed; theta from the y (normal) direction, heading toward y=0
    vx=np.sin(theta)/np.sqrt(gpar(y0)); vy=-np.cos(theta)/np.sqrt(gperp)
    return [0.0,y0,vx,vy]
s0=launch(0.15,np.deg2rad(60))
sol=solve_ivp(rhs,(0,1.3),s0,rtol=1e-12,atol=1e-14,dense_output=True)
x_,y_,vx_,vy_=sol.y
spd=np.sqrt(gpar(y_)*vx_**2+gperp*vy_**2); sin_th=np.sqrt(gpar(y_))*np.abs(vx_)/spd
lam=np.minimum(gpar(y_),gperp)
cands={'sqrt(gpar) sin':np.sqrt(gpar(y_))*sin_th,'gpar sin^2':gpar(y_)*sin_th**2,'lmin sin^2':lam*sin_th**2,'sin^2/lmin':sin_th**2/lam}
drift={k_:(v.max()-v.min())/np.mean(v) for k_,v in cands.items()}
print('     drifts:',{k_:f'{v:.2e}' for k_,v in drift.items()})
rep('R2 Clairaut: g_par sin^2(theta) conserved along a free geodesic',drift['gpar sin^2'],1e-8)
rep('R3 opposite: index = lambda_min is NOT conserved here (kill if conserved)',1e-3/drift['lmin sin^2'],1.0)
rep('R4 opposite: index = 1/lambda_min is NOT conserved (kill if conserved)',1e-3/drift['sin^2/lmin'],1.0)
print('     y crosses the lambda_min branch switch sin y = 1/k at y =',np.arcsin(1/K).round(5),'; trajectory y range',y_.min().round(3),y_.max().round(3))
# total internal reflection: launch inward from y=0.60 at 35 deg; C = gpar sin^2; turning where gpar = C
s0=launch(0.60,np.deg2rad(35)); C=gpar(0.60)*np.sin(np.deg2rad(35))**2
ypred=np.arcsin(np.sqrt(4*C)/K)
ev=lambda t,s: s[3]; ev.terminal=True; ev.direction=1
sol=solve_ivp(rhs,(0,5),s0,rtol=1e-13,atol=1e-15,events=ev,dense_output=True)
ymin=sol.y_events[0][0][1] if len(sol.y_events[0]) else sol.y[1].min()
st=sol.y_events[0][0]; sin_t=np.sqrt(gpar(st[1]))*abs(st[2])/np.sqrt(gpar(st[1])*st[2]**2+gperp*st[3]**2)
rep('R5 total internal reflection: integrated turning point = predicted root of g_par = C',abs(ymin-ypred),1e-8,f'C={C:.8f} predicted {ypred:.8f} observed {ymin:.8f}, sin(theta) at turn {sin_t:.8f}')
reach=[]
for th in (20,35,50,65,89):
    s0=launch(0.6,np.deg2rad(th)); sol=solve_ivp(rhs,(0,5),s0,rtol=1e-12,atol=1e-14,events=ev); reach.append(sol.y_events[0][0][1] if len(sol.y_events[0]) else 0)
rep('R6 no oblique geodesic reaches the stratum y=0 (min turning y over 20..89 deg)',1e-6/min(reach),1.0,f'turning y {np.round(reach,4)}')
s0=launch(0.6,0.0); sol=solve_ivp(rhs,(0,1.2),s0,rtol=1e-12,atol=1e-14)
rep('R6b opposite: the exactly normal geodesic (C=0) does reach y=0',max(0,sol.y[1].min()),1e-9,f'min y {sol.y[1].min():.2e}')
# two kinds of stratum: polar (kernel along the locus) vs fold (kernel transverse)
# polar: the family above at y=0: Pi constant along the locus
Pd=np.abs(Pm(0.0,0.0)-Pm(1.0,0.0)).max()
rep('R7 polar stratum: Pi does not move along the locus (kernel = along-interface direction)',Pd,1e-14)
# fold: theta = y + 1.1, phi = x^2/2 on the Bloch sphere; locus x=0 has kernel along d/dx (transverse), Pi moves along it
def Pf(x,y):
    th=y+1.1; ph=x*x/2; d=np.array([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
    return 0.5*(np.eye(2)-(d[0]*sx+d[1]*sy+d[2]*sz))
gf=lambda x,y: gnum.__wrapped__(x,y) if hasattr(gnum,'__wrapped__') else None
h=1e-5; Vx=(Pf(h,0.3)-Pf(-h,0.3))/(2*h); Vy=(Pf(0,0.3+h)-Pf(0,0.3-h))/(2*h)
rep('R8 fold stratum: metric degenerates at x=0 (d_x Pi = 0) ...',np.abs(Vx).max(),1e-9)
rep('R8b ... while Pi moves along the locus (|d_y Pi| > 0)',1e-3/np.abs(Vy).max(),1.0,f'|d_y Pi| = {np.abs(Vy).max():.3f}')
print('SUMMARY',sum(R),'/',len(R))
