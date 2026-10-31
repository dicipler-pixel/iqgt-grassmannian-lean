# SCRIPT: IQGT-FIG-S12
# fig_cuts: BEC, every cut (saved numbers): single eigen-image vs subspace by rank; measured/predicted vs relative gap.
# fig_corridor: tip and tail against the inner gap (computed); metric speed under the gap ceiling (computed).
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
BG='#0c1220'; FG='#b9dcff'; HEAD='#4cc3ff'; MUTED='#7088a8'; LINE='#1a2740'; WARM='#e0b85c'; VIO='#a99bff'; ROSE='#ff7a85'
plt.rcParams.update({'svg.fonttype':'none','font.family':'serif','font.serif':['Newsreader','DejaVu Serif'],'font.size':11,
 'axes.facecolor':BG,'figure.facecolor':BG,'savefig.facecolor':BG,'text.color':FG,'axes.labelcolor':FG,'axes.edgecolor':MUTED,
 'xtick.color':MUTED,'ytick.color':MUTED,'axes.grid':False,'axes.spines.top':False,'axes.spines.right':False})
d=json.load(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','data','bec_numbers.json')))
fig,axs=plt.subplots(1,2,figsize=(11.2,4.2),gridspec_kw={'width_ratios':[1.25,1]})
ax=axs[0]; R=d['T1'][0]['rows']; k=np.array([x['k'] for x in R])
for c in (2,7,10): ax.axvspan(c-0.5,c+1.5,color=VIO,alpha=0.10)
ax.bar(k-0.18,[x['vector_deg'] for x in R],0.36,color=ROSE,label='single eigen-image')
ax.bar(k+0.18,[x['projector_deg'] for x in R],0.36,color=HEAD,label='rank-k subspace')
for c,lab in ((2,'pair 2,3'),(7,'pair 7,8'),(10,'pair 10,11')): ax.text(c+0.5,86,lab,ha='center',color=VIO,fontsize=9)
ax.set_xticks(k); ax.set_xlabel('rank k (odd against even shots)'); ax.set_ylabel('angle between halves (deg)'); ax.set_ylim(0,92)
ax.legend(frameon=False,fontsize=9.5,loc='upper left',bbox_to_anchor=(0.0,0.93))
ax=axs[1]
for r,col,mk in zip(d['T1'][:3],(HEAD,WARM,VIO),('o','s','^')):
    g=[x['rel_gap'] for x in r['rows']]; q=[x['ratio'] for x in r['rows']]
    ax.semilogy(g,q,mk,color=col,ms=6,label=r['label'])
ax.axvline(0.10,color=MUTED,ls=':',lw=1); ax.axhspan(0.93,1.12,color=HEAD,alpha=0.12)
ax.text(0.105,0.13,'resolved gap',color=MUTED,fontsize=9); ax.text(0.30,1.45,'18 of 18 in 0.93–1.12',color=HEAD,fontsize=9.5)
ax.set_xlabel('relative gap at the cut'); ax.set_ylabel('measured / predicted'); ax.set_ylim(0.08,12); ax.legend(frameon=False,fontsize=8.5,loc='upper right')
fig.tight_layout(); fig.savefig('fig_cuts.svg',bbox_inches='tight'); fig.savefig('fig_cuts.png',dpi=110,bbox_inches='tight'); plt.close(fig)
# corridor
rng=np.random.default_rng(4); dins=np.logspace(-2.2,0,14); tips=[];tails=[]
for din in dins:
    w=np.array([0,din,2*din,2*din+3.0,2*din+3.5,2*din+4.0,2*din+5.0]); n=7; a=[];b=[]
    for t in range(200):
        U,_=np.linalg.qr(rng.normal(size=(n,n))); A=U@np.diag(w)@U.T
        E=rng.normal(size=(n,n)); E=(E+E.T)/2; E*=0.005/np.linalg.norm(E,2)
        w2,U2=np.linalg.eigh(A+E); P=U[:,:3]@U[:,:3].T; P2=U2[:,:3]@U2[:,:3].T
        b.append(np.degrees(np.arcsin(np.linalg.norm((np.eye(n)-P2)@P,2)))); a.append(np.degrees(np.arcsin(min(1,np.sqrt(max(0,1-(U[:,1]@U2[:,1])**2))))))
    tips.append(np.median(a)); tails.append(np.median(b))
fig,axs=plt.subplots(1,2,figsize=(11.0,4.0))
ax=axs[0]; ax.loglog(dins,tips,'o-',color=ROSE,lw=2,label='tip: one vector inside the cluster')
ax.loglog(dins,tails,'s-',color=HEAD,lw=2,label='tail: the cluster subspace')
ax.loglog(dins,np.degrees(0.005/dins),color=ROSE,ls=':',lw=1); ax.axhline(np.degrees(0.005/3.0),color=HEAD,ls=':',lw=1)
ax.text(0.0065,0.115,'‖E‖ / outer gap',color=HEAD,fontsize=9); ax.text(0.03,5.0,'‖E‖ / inner gap',color=ROSE,fontsize=9)
ax.set_xlabel('inner gap (outer gap fixed)'); ax.set_ylabel('angle moved (deg)'); ax.legend(frameon=False,fontsize=9,loc='upper right')
ax=axs[1]
H0=np.diag([0,1.2,3.0]).astype(float); H1=np.array([[0,1,0.3],[1,0,0.8],[0.3,0.8,0]]); s=np.linspace(-2,2,400); sp=[];ce=[]
for x in s:
    f=lambda y:(lambda w,U:(U[:,:1]@U[:,:1].T,w))(*np.linalg.eigh(H0+y*H1)); P,w=f(x); V=(f(x+1e-6)[0]-f(x-1e-6)[0])/2e-6
    sp.append(np.sqrt(0.5*np.trace(V@V))); ce.append(np.linalg.norm(H1,2)/(w[1]-w[0]))
ax.semilogy(s,ce,color=WARM,lw=2,label='ceiling  √k ‖∂Σ‖ / gap'); ax.semilogy(s,sp,color=HEAD,lw=2,label='metric speed  √g')
ax.fill_between(s,sp,ce,color=WARM,alpha=0.08)
ax.set_xlabel('drive parameter'); ax.set_ylabel('speed of the subspace'); ax.legend(frameon=False,fontsize=9.5,loc='upper left')
fig.tight_layout(); fig.savefig('fig_corridor.svg',bbox_inches='tight'); fig.savefig('fig_corridor.png',dpi=110,bbox_inches='tight')
print('ok', np.all(np.array(sp)<=np.array(ce)))
