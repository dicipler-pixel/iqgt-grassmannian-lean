# SCRIPT: IQGT-FIG-APPENDIX
# (a) fig_narrow: geodesics cross where the gap narrows (index rises) vs turn at a polar stratum (index falls).
# (b) fig_saturation: |Omega| against the bound for two-band, spin-1 lowest band, generic three-band.
# (c) fig_krylov: regularised Krylov metric error against eta.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
BG='#0c1220'; FG='#b9dcff'; HEAD='#4cc3ff'; MUTED='#7088a8'; LINE='#1a2740'; WARM='#e0b85c'; VIO='#a99bff'; ROSE='#ff7a85'
plt.rcParams.update({'svg.fonttype':'none','font.family':'serif','font.serif':['Newsreader','DejaVu Serif'],'font.size':11,
 'axes.facecolor':BG,'figure.facecolor':BG,'savefig.facecolor':BG,'text.color':FG,'axes.labelcolor':FG,'axes.edgecolor':MUTED,
 'xtick.color':MUTED,'ytick.color':MUTED,'axes.grid':False,'axes.spines.top':False,'axes.spines.right':False})
def geo(gpar,dgpar,gyy,dgyy,y0,th,T):
    rhs=lambda t,s:[s[2],s[3],-(dgpar(s[1])/gpar(s[1]))*s[2]*s[3],(dgpar(s[1])*s[2]**2-dgyy(s[1])*s[3]**2)/(2*gyy(s[1]))]
    s0=[0,y0,np.sin(th)/np.sqrt(gpar(y0)),-np.cos(th)/np.sqrt(gyy(y0))]
    sol=solve_ivp(rhs,(0,T),s0,rtol=1e-11,atol=1e-13,dense_output=True,events=None); t=np.linspace(0,T,800); return sol.sol(t)
fig,axs=plt.subplots(1,2,figsize=(10.6,4.0))
K=3.0; dm=2.5
gp=lambda y:K*K/(4*(1+(dm*y)**2)); dgp=lambda y:-K*K*2*dm*dm*y/(4*(1+(dm*y)**2)**2)
gy=lambda y:dm**2/(4*(1+(dm*y)**2)**2); dgy=lambda y:-dm**2*4*dm*dm*y/(4*(1+(dm*y)**2)**3)
ax=axs[0]; ax.axhspan(-0.06,0.06,color=HEAD,alpha=0.10)
for th in (35,50,65):
    x,y,_,_=geo(gp,dgp,gy,dgy,0.8,np.deg2rad(th),4.0)
    ax.plot(x,y,color=HEAD,lw=1.8)
    yt=np.sqrt((1+(dm*0.8)**2)/np.sin(np.deg2rad(th))**2-1)/dm
    ax.axhline(-yt,color=MUTED,lw=0.8,ls=':'); ax.axhline(yt,color=MUTED,lw=0.8,ls=':')
ax.text(1.32,-0.32,'gap narrowest: index largest',color=HEAD,fontsize=9.5)
ax.set_xlabel('x (along the interface)'); ax.set_ylabel('y'); ax.set_title('gap narrows: trajectories cross and are guided',color=FG,fontsize=10.5)
Kp=3.0; gp2=lambda y:Kp*Kp*np.sin(y)**2/4; dgp2=lambda y:Kp*Kp*np.sin(2*y)/4; gy2=lambda y:0.25; dgy2=lambda y:0.0
ax=axs[1]; ax.axhline(0,color=ROSE,lw=2)
for th in (20,35,50,65):
    x,y,_,_=geo(gp2,dgp2,gy2,dgy2,0.6,np.deg2rad(th),1.0); m=y>0.0
    ax.plot(x[m],y[m],color=WARM,lw=1.8)
ax.text(0.02,0.03,'polar stratum: subspace stops rotating',color=ROSE,fontsize=9.5)
ax.set_xlabel('x (along the interface)'); ax.set_ylim(-0.05,0.65); ax.set_title('subspace freezes: trajectories turn back',color=FG,fontsize=10.5)
fig.tight_layout(); fig.savefig('fig_narrow.svg',bbox_inches='tight'); fig.savefig('fig_narrow.png',dpi=110,bbox_inches='tight'); plt.close(fig)
# (b)
rng=np.random.default_rng(5)
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.diag([1,-1]).astype(complex)
Sx=np.array([[0,1,0],[1,0,1],[0,1,0]])/np.sqrt(2); Sy=np.array([[0,-1j,0],[1j,0,-1j],[0,1j,0]])/np.sqrt(2); Sz=np.diag([1,0,-1]).astype(complex)
def herm(n):
    a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return (a+a.conj().T)/2
H0,H1,H2=herm(3),herm(3),herm(3)
def low(M): w,U=np.linalg.eigh(M); return U[:,:1]@U[:,:1].conj().T
fams={'two-band':lambda p:low(p[0]*sx+p[1]*sy+p[2]*sz),'spin-1, lowest band':lambda p:low(p[0]*Sx+p[1]*Sy+p[2]*Sz),'generic three-band':lambda p:low(H0+p[0]*H1+p[1]*H2+p[2]*0*H1)}
cols={'two-band':HEAD,'spin-1, lowest band':VIO,'generic three-band':WARM}
fig,ax=plt.subplots(figsize=(5.6,4.6))
for name,f in fams.items():
    B=[];O=[]
    for t in range(150):
        p=rng.normal(size=3); v1=rng.normal(size=3); v2=rng.normal(size=3)
        if name=='generic three-band': p=p*0.6; v1[2]=0; v2[2]=0
        h=1e-6; P=f(p); V=(f(p+h*v1)-f(p-h*v1))/(2*h); W=(f(p+h*v2)-f(p-h*v2))/(2*h)
        g11=0.5*np.trace(V@V).real; g22=0.5*np.trace(W@W).real; g12=0.5*np.trace(V@W).real
        B.append(2*np.sqrt(max(g11*g22-g12**2,0))); O.append(abs((1j*np.trace(P@(V@W-W@V))).real))
    ax.loglog(B,O,'o',ms=4,color=cols[name],alpha=0.8,label=name)
xx=np.logspace(-4,2,10); ax.loglog(xx,xx,color=MUTED,lw=1,ls=':'); ax.text(3e-3,6e-3,'|Ω| = bound',color=MUTED,fontsize=9,rotation=38)
ax.set_xlabel('bound  2√(g g − g²)'); ax.set_ylabel('|Ω|'); ax.set_ylim(1e-5,1e2); ax.set_xlim(1e-4,1e2); ax.legend(frameon=False,fontsize=9.5,loc='lower right')
fig.tight_layout(); fig.savefig('fig_saturation.svg',bbox_inches='tight'); fig.savefig('fig_saturation.png',dpi=110,bbox_inches='tight'); plt.close(fig)
# (c) from app_checks_v2 numbers
eta=np.array([1e-1,1e-2,1e-3,1e-6]); eg=np.array([7.6e-2,8.3e-4,8.3e-6,8.3e-12]); eo=np.array([1.2e-1,1.4e-3,1.4e-5,1.4e-11])
fig,ax=plt.subplots(figsize=(5.6,4.0))
ax.loglog(eta,eg,'o-',color=HEAD,lw=2,label='metric g'); ax.loglog(eta,eo,'s-',color=VIO,lw=2,label='curvature Ω')
ax.loglog(eta,8.3*eta**2,color=MUTED,ls=':',lw=1); ax.text(2e-4,5e-6,'∝ η²',color=MUTED,fontsize=10)
ax.set_xlabel('regulator η'); ax.set_ylabel('relative error vs exact pair sum'); ax.legend(frameon=False,fontsize=9.5,loc='upper left')
fig.tight_layout(); fig.savefig('fig_krylov.svg',bbox_inches='tight'); fig.savefig('fig_krylov.png',dpi=110,bbox_inches='tight')
