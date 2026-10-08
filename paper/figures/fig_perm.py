# SCRIPT: IQGT-FIG-PERM
# fig_perm: principal angles between graph(I) and graph(pi) for five swaps vs one 10-cycle (n=20): equal Hamming,
# different geodesic length; and the sandwich d_H/2 <= (1/n) sum theta^2 <= pi^2/8 d_H over random pairs.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
BG='#0c1220'; FG='#b9dcff'; HEAD='#4cc3ff'; MUTED='#7088a8'; LINE='#1a2740'; WARM='#e0b85c'; VIO='#a99bff'; ROSE='#ff7a85'
plt.rcParams.update({'svg.fonttype':'none','font.family':'serif','font.serif':['Newsreader','DejaVu Serif'],'font.size':11,
 'axes.facecolor':BG,'figure.facecolor':BG,'savefig.facecolor':BG,'text.color':FG,'axes.labelcolor':FG,'axes.edgecolor':MUTED,
 'xtick.color':MUTED,'ytick.color':MUTED,'axes.grid':False,'axes.spines.top':False,'axes.spines.right':False})
def perm(p):
    n=len(p); U=np.zeros((n,n)); U[p,np.arange(n)]=1; return U
def angles(U,V):
    s=np.linalg.svd((np.eye(len(U))+U.T@V)/2,compute_uv=False); return np.sort(np.arccos(np.clip(s,0,1)))[::-1]
n=20; sw=list(range(n))
for i in range(0,10,2): sw[i],sw[i+1]=sw[i+1],sw[i]
cy=list(range(n)); cy[:10]=list(range(1,10))+[0]
tA,tB=angles(np.eye(n),perm(sw)),angles(np.eye(n),perm(cy))
fig,axs=plt.subplots(1,3,figsize=(12,4.0),gridspec_kw={'width_ratios':[1,1,1.1]})
for ax,th,col,lab in ((axs[0],tA,ROSE,'five swaps'),(axs[1],tB,HEAD,'one 10-cycle')):
    ax.set_aspect('equal'); ax.axis('off')
    for t in th[:10]:
        ax.plot([0,np.cos(t)],[0,np.sin(t)],color=col,lw=2.2,alpha=0.85)
    a=np.linspace(0,np.pi/2,100); ax.plot(np.cos(a),np.sin(a),color=LINE,lw=1)
    ax.text(0.5,-0.18,lab+('   (5 at 90°, 5 at 0°)' if lab=='five swaps' else ''),ha='center',color=FG,fontsize=11,transform=ax.transData)
    ax.text(0.5,-0.36,'Σ sin²θ = %.1f    √Σθ² = %.3f'%(np.sum(np.sin(th)**2),np.sqrt(np.sum(th**2))),ha='center',color=MUTED,fontsize=9.5)
    ax.set_xlim(-0.1,1.15); ax.set_ylim(-0.45,1.1)
axs[0].text(0.02,1.05,'principal angles between graph(I) and graph(π)',color=MUTED,fontsize=9)
rng=np.random.default_rng(2); X=[];Y=[]
for t in range(1500):
    m=rng.integers(2,13); U=perm(rng.permutation(m)); V=perm(rng.permutation(m))
    dH=np.mean(np.any(U!=V,axis=0)); th=angles(U,V); X.append(dH); Y.append(np.sum(th**2)/m)
ax=axs[2]; ax.plot(X,Y,'o',ms=3,color=VIO,alpha=0.5)
xx=np.linspace(0,1,10); ax.plot(xx,xx/2,color=MUTED,ls=':',lw=1.2); ax.plot(xx,np.pi**2/8*xx,color=MUTED,ls=':',lw=1.2)
ax.text(0.62,0.22,'d_H / 2',color=MUTED,fontsize=9.5); ax.text(0.05,0.36,'π²/8 · d_H',color=MUTED,fontsize=9.5)
ax.plot([0.5],[np.sum(tA**2)/n],'s',color=ROSE,ms=8); ax.plot([0.5],[np.sum(tB**2)/n],'o',color=HEAD,ms=8); ax.text(0.53,np.sum(tA**2)/n+0.02,'swaps',color=ROSE,fontsize=9); ax.text(0.53,np.sum(tB**2)/n-0.06,'cycle',color=HEAD,fontsize=9)
ax.set_xlabel('Hamming distance d_H'); ax.set_ylabel('(1/n) Σ θ²  (squared geodesic, per point)')
fig.tight_layout(); fig.savefig('fig_perm.svg',bbox_inches='tight'); fig.savefig('fig_perm.png',dpi=110,bbox_inches='tight')
