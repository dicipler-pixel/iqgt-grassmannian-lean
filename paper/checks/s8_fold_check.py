# SCRIPT: IQGT-REBUILD-S8-FOLD-01
# Kill: some oblique geodesic of the fold metric g = 1/4 [ s(y)^2 x^2 dx^2 + dy^2 ], s = sin(y+1.1), fails to reach the fold x = 0.
import numpy as np
from scipy.integrate import solve_ivp
s=lambda y: np.sin(y+1.1); ds=lambda y: np.cos(y+1.1)
# metric components: gxx = s^2 x^2 /4, gyy = 1/4
def rhs(t,u):
    x,y,vx,vy=u
    gxx=s(y)**2*x**2/4; dgx=s(y)**2*x/2; dgy=s(y)*ds(y)*x**2/2; gyy=0.25
    ax=-(dgx*vx*vx/2 + dgy*vx*vy)/gxx      # Christoffel for diagonal metric
    ay=(dgy*vx*vx/2)/gyy
    return [vx,vy,ax,ay]
hit=lambda t,u: u[0]-1e-6; hit.terminal=True
res=[]
for th in (10,30,50,70,85):
    x0,y0=0.5,0.3; gxx=s(y0)**2*x0**2/4
    t=np.deg2rad(th); vx=-np.cos(t)/np.sqrt(gxx); vy=np.sin(t)/np.sqrt(0.25)   # angle from the normal (x direction)
    sol=solve_ivp(rhs,(0,20),[x0,y0,vx,vy],events=hit,rtol=1e-10,atol=1e-12)
    res.append((th,len(sol.t_events[0])>0, sol.y[0].min()))
print(res)
print('PASS' if all(r[1] for r in res) else 'FAIL','all oblique geodesics reach the fold')
