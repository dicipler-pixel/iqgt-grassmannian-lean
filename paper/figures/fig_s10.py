# SCRIPT: IQGT-FIG-S10
# (a) fig_direction: a direction exists only where a gap exists (heat-engine readout, printed covariances + null bands).
# (b) fig_gaplaw: vector vs subspace. Synthetic voting scatter vs first-order law; BEC closing cuts (ledger v2 ranges).
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
BG='#0c1220'; FG='#b9dcff'; HEAD='#4cc3ff'; MUTED='#7088a8'; LINE='#1a2740'; WARM='#e0b85c'; VIO='#a99bff'; ROSE='#ff7a85'
plt.rcParams.update({'svg.fonttype':'none','font.family':'serif','font.serif':['Newsreader','DejaVu Serif'],'font.size':11,
 'axes.facecolor':BG,'figure.facecolor':BG,'savefig.facecolor':BG,'text.color':FG,'axes.labelcolor':FG,'axes.edgecolor':MUTED,
 'xtick.color':MUTED,'ytick.color':MUTED,'axes.grid':False,'axes.spines.top':False,'axes.spines.right':False})
rng=np.random.default_rng(12)
l_lo,l_hi=0.00051245,0.00055124
ev=np.array([[0.22694602,-0.97390734],[-0.97390734,-0.22694602]]); C=ev@np.diag([l_lo,l_hi])@ev.T
angs=[];rats=[]
for s in range(300):
    x=rng.multivariate_normal([0,0],C,size=19899); d=np.linalg.norm(x-x.mean(0),axis=1)
    co=x[d<=np.percentile(d,50)]; ed=x[d>=np.percentile(d,90)]
    vc=np.linalg.eigh(np.cov(co.T))[1][:,-1]; w,U=np.linalg.eigh(np.cov(ed.T))
    angs.append(np.degrees(np.arccos(min(1,abs(vc@U[:,-1]))))); rats.append(w[1]/w[0])
fig,axs=plt.subplots(1,3,figsize=(11.2,3.7),gridspec_kw={'width_ratios':[1.05,1,1]})
ax=axs[0]; ax.set_aspect('equal'); ax.axis('off')
def ell(ax,lo,hi,v,col,lab,scale):
    a=np.degrees(np.arctan2(v[1],v[0])); e=Ellipse((0,0),2*scale*np.sqrt(hi),2*scale*np.sqrt(lo),angle=a,fill=False,ec=col,lw=2.2); ax.add_patch(e)
ell(ax,l_lo,l_hi,ev[:,1],HEAD,'core',2.0); ell(ax,0.00481701,0.01067572,np.array([-0.46464346,0.88549786]),WARM,'tail',1.0)
ax.set_xlim(-0.12,0.12); ax.set_ylim(-0.12,0.12)
ax.text(0,-0.115,'ground-state readout cloud\ncore (inner half)  ratio 1.08\ntail (outer tenth)  ratio 2.22',ha='center',va='top',fontsize=9.5,color=FG)
ax.text(-0.012,0.006,'core',color=HEAD,fontsize=10); ax.text(0.045,0.085,'tail',color=WARM,fontsize=10)
ax=axs[1]; ax.hist(angs,bins=30,color=MUTED,alpha=0.85); ax.axvline(75.43,color=ROSE,lw=2.2)
ax.set_xlabel('tail axis vs core axis (deg)'); ax.set_yticks([]); ax.text(74,ax.get_ylim()[1]*0.92,'measured\n75.4°',color=ROSE,ha='right',fontsize=9.5)
ax.set_title('angle: inside the noise',color=FG,fontsize=10.5)
ax=axs[2]; ax.hist(rats,bins=25,color=MUTED,alpha=0.85); ax.axvline(2.2163,color=HEAD,lw=2.2)
ax.set_xlabel('tail eigenvalue ratio'); ax.set_yticks([]); ax.set_xlim(1.0,2.35); ax.text(2.18,ax.get_ylim()[1]*0.92,'measured\n2.22',color=HEAD,ha='right',fontsize=9.5)
ax.set_title('shape: far outside it',color=FG,fontsize=10.5)
fig.tight_layout(); fig.savefig('fig_direction.svg',bbox_inches='tight'); fig.savefig('fig_direction.png',dpi=110,bbox_inches='tight'); plt.close(fig)
# (b)
rs=np.array([1.01,1.02,1.04,1.08,1.15,1.3,1.6,2.2,4]); N=9949; sim=[]
for r in rs:
    a=[]
    for t in range(300):
        x=rng.normal(size=(N,2))*np.sqrt([r,1.0]); w,U=np.linalg.eigh(np.cov(x.T)); a.append(np.arccos(min(1,abs(U[0,-1]))))
    sim.append(np.degrees(np.sqrt(np.mean(np.square(a)))))
rr=np.linspace(1.008,4,400); pred=np.degrees(np.sqrt(rr)/((rr-1)*np.sqrt(N)))
fig,axs=plt.subplots(1,2,figsize=(10.6,3.9))
ax=axs[0]; ax.loglog(rr-1,np.minimum(pred,60),color=HEAD,lw=2,label='first-order gap law'); ax.loglog(rs-1,sim,'o',color=WARM,ms=7,label='simulated, 9,949 samples')
ax.axhline(np.degrees(np.pi/(2*np.sqrt(3))),color=MUTED,ls=':',lw=1); ax.text(0.25,60,'no gap: direction uniform',color=MUTED,fontsize=9)
ax.set_xlabel('relative gap  (λ₁ − λ₂)/λ₂'); ax.set_ylabel('leading-vector scatter (deg)'); ax.legend(frameon=False,fontsize=9.5,loc='lower left'); ax.set_ylim(0.2,80)
ax=axs[1]
rows=[('k = 3',(26,64),(5.2,8.7)),('k = 8',(12,39),(7.7,11.4)),('k = 11',(69,81),(11.8,14.2))]
for i,(lab,v,p) in enumerate(rows):
    y=2-i; ax.plot(v,[y+0.13]*2,color=ROSE,lw=9,solid_capstyle='butt'); ax.plot(p,[y-0.13]*2,color=HEAD,lw=9,solid_capstyle='butt')
    ax.text(-3,y,lab,ha='right',va='center',color=FG,fontsize=10.5)
ax.set_xlim(0,90); ax.set_ylim(-0.6,2.7); ax.set_yticks([]); ax.spines['left'].set_visible(False)
ax.set_xlabel('angle between independent halves (deg)')
ax.text(46,2.48,'single eigen-image',color=ROSE,fontsize=10); ax.text(14,2.48,'subspace',color=HEAD,fontsize=10)
fig.tight_layout(); fig.savefig('fig_gaplaw.svg',bbox_inches='tight'); fig.savefig('fig_gaplaw.png',dpi=110,bbox_inches='tight')
print('sim',np.round(sim,2))
