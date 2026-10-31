# SCRIPT: IQGT-FIG-S7S8
# (a) weights ladder: one pair sum, read with different powers of the gap. (b) three candidate indices along one geodesic.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
BG='#0c1220'; FG='#b9dcff'; HEAD='#4cc3ff'; MUTED='#7088a8'; LINE='#1a2740'; WARM='#e0b85c'; VIO='#a99bff'; ROSE='#ff7a85'
plt.rcParams.update({'svg.fonttype':'none','font.family':'serif','font.serif':['Newsreader','DejaVu Serif'],'font.size':11,
 'axes.facecolor':BG,'figure.facecolor':BG,'savefig.facecolor':BG,'text.color':FG,'axes.labelcolor':FG,'axes.edgecolor':MUTED,
 'xtick.color':MUTED,'ytick.color':MUTED,'axes.grid':False,'axes.spines.top':False,'axes.spines.right':False})
# (a) ladder
fig,ax=plt.subplots(figsize=(7.4,3.9)); ax.set_xlim(-1.6,1.6); ax.set_ylim(-0.6,2.4); ax.axis('off')
for p in (-1,0,1): ax.axvline(p,color=LINE,lw=1); ax.text(p,2.25,{-1:'Δ⁻¹',0:'Δ⁰',1:'Δ¹'}[p],ha='center',color=MUTED,fontsize=12)
ax.text(-1.55,1.45,'metric\nfamily  gₙ',color=HEAD,fontsize=11,va='center'); ax.text(-1.55,0.25,'curvature\nfamily  Ωₙ',color=VIO,fontsize=11,va='center')
ax.plot([-1,0,1],[1.45]*3,color=HEAD,lw=1.2,alpha=0.6); ax.plot([0,1],[0.25]*2,color=VIO,lw=1.2,alpha=0.6)
for x_,lab in ((-1,'Berry connection\npolarizability\n(nonlinear Hall input)'),(0,'quantum metric'),(1,'static response\n−E₀″ and\nrotation friction')):
    ax.scatter([x_],[1.45],s=140,color=HEAD,zorder=3); ax.text(x_,1.72,lab,ha='center',va='bottom',fontsize=9.2,color=FG)
for x_,lab in ((0,'Berry curvature'),(1,'orbital-moment\nkernel')):
    ax.scatter([x_],[0.25],s=140,color=VIO,zorder=3); ax.text(x_,-0.05,lab,ha='center',va='top',fontsize=9.2,color=FG)
ax.text(0,0.85,'same cross-gap pairs  ⟨0|∂Σ|n⟩,  different power of the gap',ha='center',color=WARM,fontsize=10)
fig.savefig('fig_weights.svg',bbox_inches='tight'); fig.savefig('fig_weights.png',dpi=110,bbox_inches='tight'); plt.close(fig)
# (b) candidates along one geodesic
K=3; gpar=lambda y:K*K*np.sin(y)**2/4; dg=lambda y:K*K*np.sin(2*y)/4; gperp=0.25
rhs=lambda t,s:[s[2],s[3],-(dg(s[1])/gpar(s[1]))*s[2]*s[3],(dg(s[1])/(2*gperp))*s[2]**2]
y0=0.15; th=np.deg2rad(60); s0=[0,y0,np.sin(th)/np.sqrt(gpar(y0)),-np.cos(th)/np.sqrt(gperp)]
sol=solve_ivp(rhs,(0,1.3),s0,rtol=1e-12,atol=1e-14,dense_output=True); t=np.linspace(0,1.3,800); x,y,vx,vy=sol.sol(t)
sp=np.sqrt(gpar(y)*vx**2+gperp*vy**2); sn=np.sqrt(gpar(y))*abs(vx)/sp; lm=np.minimum(gpar(y),gperp)
fig,ax=plt.subplots(figsize=(7.2,4.0))
for v,c,l in ((gpar(y)*sn**2,HEAD,'g∥ sin²θ  (index √g∥)'),(lm*sn**2,WARM,'λ_min sin²θ'),(sn**2/lm,ROSE,'sin²θ / λ_min')):
    ax.semilogy(t,v/v[0],color=c,lw=2.2,label=l)
sw=t[np.argmin(abs(np.sin(y)-1/K))]; ax.axvline(sw,color=MUTED,ls=':',lw=1); ax.text(sw+0.02,2.2,'λ_min switches branch',color=MUTED,fontsize=9)
ax.set_xlabel('affine parameter along one geodesic'); ax.set_ylabel('value / initial value'); ax.legend(frameon=False,fontsize=9.5,loc='lower center',bbox_to_anchor=(0.5,1.02),ncol=3); ax.set_ylim(3e-3,4)
ax.grid(True,color=LINE,lw=0.6)
fig.savefig('fig_candidates.svg',bbox_inches='tight'); fig.savefig('fig_candidates.png',dpi=110,bbox_inches='tight')
