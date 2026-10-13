# SCRIPT: IQGT-FIG-V7-3
# Figures new in v7: rigidity (equality cases), orders (silent directions by order), shape (three-body fold), perm7.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Circle
BG='#0c1220'; FG='#b9dcff'; HEAD='#4cc3ff'; MUTED='#7088a8'; LINE='#1a2740'; WARM='#e0b85c'; VIO='#a99bff'; ROSE='#ff7a85'; GREEN='#7fe0b0'
plt.rcParams.update({'svg.fonttype':'none','font.family':'serif','font.serif':['Newsreader','DejaVu Serif'],'font.size':11,
 'axes.facecolor':BG,'figure.facecolor':BG,'savefig.facecolor':BG,'text.color':FG,'axes.labelcolor':FG,'axes.edgecolor':MUTED,
 'xtick.color':MUTED,'ytick.color':MUTED,'axes.grid':False,'axes.spines.top':False,'axes.spines.right':False})
sx=np.array([[0,1],[1,0]]);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1])
def projI(H,I): w,U=np.linalg.eigh(H); return U[:,I]@U[:,I].conj().T
def dproj(H,dS,I,h=1e-6): return (projI(H+h*dS,I)-projI(H-h*dS,I))/(2*h)
def gapratio(H,dS,I):
    w=np.linalg.eigvalsh(H); n=len(w); Ic=[j for j in range(n) if j not in I]
    D=min(abs(w[a]-w[b]) for a in I for b in Ic); V=dproj(H,dS,I); k=len(I)
    return np.linalg.norm(V)/(np.sqrt(2*k)*np.linalg.norm(dS,2)/D)
# ---------- rigidity
fig=plt.figure(figsize=(10.5,8.2)); gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],hspace=0.42,wspace=0.28)
axs=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[1,:]),fig.add_subplot(gs[0,1])]
ax=axs[0]; res={'equal':[],'gap':[],'str':[]}
for k in (1,2,3):
    n=2*k; H=np.diag([0.]*k+[1.]*k); dS=np.zeros((n,n),complex)
    for j in range(k): dS[j,k+j]=dS[k+j,j]=1.0
    res['equal'].append(gapratio(H,dS,list(range(k))))
    n2=n+1; H2=np.diag([0.]*k+[1.]*(k-1)+[1.5]+[1.0]); dS2=np.zeros((n2,n2),complex); dS2[:n,:n]=dS
    res['gap'].append(gapratio(H2,dS2,list(range(k))))
    dS3=dS.copy(); dS3[0,k]=dS3[k,0]=0.5; res['str'].append(gapratio(H,dS3,list(range(k))) if k>1 else np.nan)
x=np.arange(3); w=0.26
ax.bar(x-w,res['equal'],w,color=HEAD,label='k equal couplings at one gap'); ax.bar(x,res['gap'],w,color=WARM,label='one partner moved to 1.5Δ'); ax.bar(x+w,res['str'],w,color=ROSE,label='unequal strengths')
ax.axhline(1,color=MUTED,ls=':',lw=1); ax.set_xticks(x); ax.set_xticklabels(['k = 1','k = 2','k = 3']); ax.set_ylim(0,1.55)
ax.set_ylabel('‖V‖ / bound'); ax.set_title('gap bound: met only at one gap, one strength',fontsize=12,color=FG); ax.legend(frameon=False,fontsize=10,loc='upper center',ncol=1,labelcolor=FG)
ax=axs[1]; psi=np.linspace(0,np.pi/2,121); beta=np.linspace(0,2*np.pi,181); P,B=np.meshgrid(beta,psi)
D=np.sin(B)**2
from matplotlib.colors import LinearSegmentedColormap
cm=LinearSegmentedColormap.from_list('dk',['#05080f','#2a2470',VIO,HEAD,WARM])
im=ax.pcolormesh(np.degrees(P),np.degrees(B),D,cmap=cm,shading='auto',vmin=0,vmax=1,rasterized=True)
ax.plot([0,360],[0,0],color=GREEN,lw=4,solid_capstyle='butt'); ax.text(180,5,'kernel: the complex line of V',ha='center',color=GREEN,fontsize=11)
ax.set_xlabel('phase β of W along V  (°)'); ax.set_ylabel('angle ψ of W off the complex line (°)')
ax.set_title('curvature deficit (6.3) = g(V,V) g(W⊥,W⊥): zero only on the complex line of V',fontsize=12,color=FG); ax.set_xticks([0,90,180,270,360]); ax.set_yticks([0,30,60,90])
cb=fig.colorbar(im,ax=ax,fraction=0.05,pad=0.02); cb.set_label('deficit / (g_VV g_WW)',color=MUTED); cb.ax.yaxis.set_tick_params(color=MUTED)
ax=axs[2]; Dl=0.7; om=1.3; ts=np.linspace(0,2,201)
nv=lambda t:np.array([np.sin(om*t),0,np.cos(om*t)]); ndv=lambda t:om*np.array([np.cos(om*t),0,-np.sin(om*t)])
def ser(Sf,Sdf,gapf):
    sp=[];dn=[];gp=[]
    for t in ts:
        V=dproj(Sf(t),Sdf(t),[0]); sp.append(np.sqrt(0.5*np.trace(V@V).real)); dn.append(np.linalg.norm(Sdf(t),2)); gp.append(gapf(t))
    sp=np.array(sp); return sp/(max(dn)/min(gp))
Sr=lambda t:Dl/2*(nv(t)[0]*sx+nv(t)[2]*sz); Srd=lambda t:Dl/2*(ndv(t)[0]*sx+ndv(t)[2]*sz)
rr=lambda t:1+0.5*np.sin(t); rd=lambda t:0.5*np.cos(t)
Sb=lambda t:rr(t)*(nv(t)[0]*sx+nv(t)[2]*sz); Sbd=lambda t:rr(t)*(ndv(t)[0]*sx+ndv(t)[2]*sz)+rd(t)*(nv(t)[0]*sx+nv(t)[2]*sz)
a=ser(Sr,Srd,lambda t:Dl); b=ser(Sb,Sbd,lambda t:2*rr(t))
ax.plot(ts,a,color=HEAD,lw=2.2,label='rigid rotation, fixed gap'); ax.plot(ts,b,color=WARM,lw=2.2,label='same rotation, breathing gap')
ax.axhline(1,color=MUTED,ls=':',lw=1); ax.set_ylim(0,1.15); ax.set_xlabel('time'); ax.set_ylabel('√g  /  (√k sup‖∂Σ‖ / Δ_min)')
ax.set_title('speed limit: met only by a rigid rotation',fontsize=12,color=FG); ax.legend(frameon=False,fontsize=10,loc='lower right',labelcolor=FG)
fig.tight_layout(); fig.savefig('fig_rigidity.svg',bbox_inches='tight'); fig.savefig('fig_rigidity.png',dpi=110,bbox_inches='tight'); plt.close(fig)
# ---------- orders
def bloch(th,ph): return projI(np.cos(ph)*np.sin(th)*sx+np.sin(ph)*np.sin(th)*sy+np.cos(th)*sz,[0])
dist=lambda A,B:np.linalg.norm(A-B)/np.sqrt(2)
eps=np.logspace(-3,-0.5,40)
P0=bloch(0.3+1.1,0); fold=[dist(bloch(0.3+1.1,e*e/2),P0) for e in eps]
Q0=bloch(0.3,0); gen=[dist(bloch(0.3+e,0),Q0) for e in eps]
R0=bloch(0.0,0); pol=[dist(bloch(0.0,3*e),R0) for e in eps]
fig,ax=plt.subplots(1,1,figsize=(6.6,4.4))
ax.loglog(eps,gen,color=HEAD,lw=2.4,label='moving direction: d ∝ ε')
ax.loglog(eps,fold,color=WARM,lw=2.4,label='across a fold: d ∝ ε²  (silent at first order)')
P3=bloch(1.4,0); odd=[dist(bloch(1.4,e**3/6),P3) for e in eps]
sel=eps>=3e-3; ax.loglog(eps[sel],np.array(odd)[sel],color=VIO,lw=2.0,ls='-.',label='odd order (ε³): passes through, removable')
floor=1e-9; ax.loglog(eps,[floor]*len(eps),color=ROSE,lw=2.4,ls='--',label='along a polar stratum: d = 0 at every order')
g2=0.25*np.sin(1.4)**2; ax.loglog(eps,eps**2/2*np.sqrt(g2),color=MUTED,ls=':',lw=1.2); ax.text(1.1e-3,3e-5,'leading term (8.3), m = 2',color=MUTED,fontsize=9)
ax.set_xlabel('step ε along the direction'); ax.set_ylabel('chordal distance d(ε)'); ax.set_ylim(5e-10,1)
ax.text(1.1e-3,1.6e-9,'(drawn at the floor; computed drift at rounding level)',color=ROSE,fontsize=8.5)
ax.legend(frameon=False,fontsize=9,loc='upper left',labelcolor=FG)
fig.tight_layout(); fig.savefig('fig_orders.svg',bbox_inches='tight'); fig.savefig('fig_orders.png',dpi=110,bbox_inches='tight'); plt.close(fig)
# ---------- shape
fig,axs=plt.subplots(1,2,figsize=(11.5,4.6),gridspec_kw={'width_ratios':[1.15,1]})
ax=axs[0]; ax.set_aspect('equal'); ax.axis('off')
tilt=0.32
def proj3(v): return np.array([v[0],v[2]*np.cos(tilt)+v[1]*np.sin(tilt)]), (v[1]*np.cos(tilt)-v[2]*np.sin(tilt))
ax.add_patch(Circle((0,0),1,fill=False,ec=MUTED,lw=1.2))
for lat in (-60,-30,30,60):
    a=np.linspace(0,2*np.pi,200); b=np.radians(lat); pts=[proj3([np.cos(b)*np.cos(t),np.cos(b)*np.sin(t),np.sin(b)]) for t in a]
    xy=np.array([p[0] for p in pts]); dep=np.array([p[1] for p in pts]); ax.plot(xy[dep<0,0],xy[dep<0,1],'.',ms=0.6,color=LINE)
a=np.linspace(0,2*np.pi,400); pts=[proj3([np.cos(t),np.sin(t),0]) for t in a]; xy=np.array([p[0] for p in pts]); dep=np.array([p[1] for p in pts])
ax.plot(np.where(dep<=0,xy[:,0],np.nan),np.where(dep<=0,xy[:,1],np.nan),color=WARM,lw=3)
ax.plot(np.where(dep>0,xy[:,0],np.nan),np.where(dep>0,xy[:,1],np.nan),color=WARM,lw=1.2,ls=':')
ax.text(-1.02,-0.42,'collinear triangles (w₃ = 0)',color=WARM,fontsize=9.5)
up=np.array([0.55,-0.45,0.70]); up/=np.linalg.norm(up); dn=up*np.array([1,1,-1])
pu,_=proj3(up); pd,_=proj3(dn)
ax.plot(*pu,'o',color=HEAD,ms=9); ax.plot(*pd,'o',color=ROSE,ms=9)
ax.text(pu[0]+0.06,pu[1]+0.05,'triangle',color=HEAD,fontsize=9.5); ax.text(pd[0]+0.06,pd[1]-0.1,'its mirror image',color=ROSE,fontsize=9.5)
eq=np.array([np.cos(-0.7),np.sin(-0.7),0]); pe,_=proj3(eq); top,_=proj3(np.array([np.cos(-0.7)*np.cos(0.35),np.sin(-0.7)*np.cos(0.35),np.sin(0.35)])); bot,_=proj3(np.array([np.cos(-0.7)*np.cos(0.35),np.sin(-0.7)*np.cos(0.35),-np.sin(0.35)]))
ax.add_patch(FancyArrowPatch(bot,top,arrowstyle='<->',mutation_scale=14,color=GREEN,lw=2))
ax.text(1.06,0.32,'kernel of the\ndistance metric:\nacross the circle',color=GREEN,fontsize=9.5)
al,_=proj3(np.array([np.cos(-0.45),np.sin(-0.45),0])); ax.add_patch(FancyArrowPatch(pe,al,arrowstyle='->',mutation_scale=13,color=HEAD,lw=1.8))
ax.text(1.06,-0.22,'Π still moves\nalong it',color=HEAD,fontsize=9.5)
ax.set_xlim(-1.15,1.5); ax.set_ylim(-1.15,1.15); ax.set_title('shape sphere = Gr(1, ℂ²)',color=FG,fontsize=11)
ax=axs[1]; ax.set_aspect('equal'); ax.axis('off')
T=np.array([[0,0],[1,0],[0.3,0.8]]); M=T*np.array([1,-1])+np.array([0,-0.15])
ax.add_patch(Polygon(T+np.array([0,0.2]),closed=True,fill=False,ec=HEAD,lw=2.2)); ax.add_patch(Polygon(M-np.array([0,0.05]),closed=True,fill=False,ec=ROSE,lw=2.2))
for P_,c,off in ((T+np.array([0,0.2]),HEAD,1),(M-np.array([0,0.05]),ROSE,-1)):
    for i,j in ((0,1),(0,2),(1,2)):
        mid=(P_[i]+P_[j])/2; ax.text(mid[0]+0.04,mid[1]+0.05*off,'%.3f'%np.linalg.norm(T[i]-T[j]),color=c,fontsize=9)
ax.plot([-0.2,1.3],[0.075,0.075],color=WARM,lw=1,ls=':'); ax.text(1.0,0.11,'mirror',color=WARM,fontsize=9)
ax.text(-0.25,-1.25,'same three distances: the distances see w₁, w₂\nbut not the sign of w₃, so the two sheets\nmeet on the collinear circle — a fold',color=FG,fontsize=9.5)
ax.set_xlim(-0.35,1.45); ax.set_ylim(-1.4,1.15); ax.set_title('a triangle and its mirror image',color=FG,fontsize=11)
fig.tight_layout(); fig.savefig('fig_shape.svg',bbox_inches='tight'); fig.savefig('fig_shape.png',dpi=110,bbox_inches='tight'); plt.close(fig)
# ---------- perm7
def perm(p):
    n=len(p); U=np.zeros((n,n)); U[p,np.arange(n)]=1; return U
def angles(U,V):
    s=np.linalg.svd((np.eye(len(U))+U.T@V)/2,compute_uv=False); return np.sort(np.arccos(np.clip(s,0,1)))[::-1]
n=20; sw=list(range(n))
for i in range(0,10,2): sw[i],sw[i+1]=sw[i+1],sw[i]
cy=list(range(n)); cy[:10]=list(range(1,10))+[0]
db=[(2*r)%9 for r in range(9)]
tA,tB,tC=angles(np.eye(n),perm(sw)),angles(np.eye(n),perm(cy)),angles(np.eye(9),perm(db))
fig=plt.figure(figsize=(11,8.4)); gs=fig.add_gridspec(2,3,height_ratios=[1,1.1],hspace=0.25)
axs=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[0,2]),fig.add_subplot(gs[1,:])]
for ax,th,col,lab,sub in ((axs[0],tA,ROSE,'five swaps (20 points)','5 at 90°, 5 at 0°'),(axs[1],tB,HEAD,'one 10-cycle (20 points)','0° to 90° in 18° steps'),(axs[2],tC,GREEN,'×2 on ℤ/9ℤ','(1 2 4 8 7 5)(3 6)(0)')):
    ax.set_aspect('equal'); ax.axis('off')
    for t in th[:10]:
        if t>1e-9: ax.plot([0,np.cos(t)],[0,np.sin(t)],color=col,lw=2.2,alpha=0.85)
    a=np.linspace(0,np.pi/2,100); ax.plot(np.cos(a),np.sin(a),color=LINE,lw=1)
    ax.text(0.5,-0.18,lab,ha='center',color=FG,fontsize=12); ax.text(0.5,-0.32,sub,ha='center',color=MUTED,fontsize=10.5)
    ax.text(0.5,-0.46,'Σ sin²θ = %.1f    √Σθ² = %.3f'%(np.sum(np.sin(th)**2),np.sqrt(np.sum(th**2))),ha='center',color=MUTED,fontsize=11)
    ax.set_xlim(-0.1,1.15); ax.set_ylim(-0.55,1.1)
axs[0].text(0.0,1.05,'principal angles between graph(I) and graph(π)',color=MUTED,fontsize=10.5)
rng=np.random.default_rng(2); X=[];Y=[]
for t in range(1500):
    m=rng.integers(2,13); U=perm(rng.permutation(m)); V=perm(rng.permutation(m))
    dH=np.mean(np.any(U!=V,axis=0)); th=angles(U,V); X.append(dH); Y.append(np.sum(th**2)/m)
ax=axs[3]; ax.plot(X,Y,'o',ms=3,color=VIO,alpha=0.45)
xx=np.linspace(0,1,10); ax.plot(xx,xx/2,color=MUTED,ls=':',lw=1.2); ax.plot(xx,np.pi**2/8*xx,color=WARM,ls='-',lw=1.3)
ax.text(0.62,0.22,'d_H / 2',color=MUTED,fontsize=9.5); ax.text(0.02,1.2,'π²/8 · d_H : met by involutions only',color=WARM,fontsize=9.2)
ax.plot([0.5],[np.sum(tA**2)/n],'s',color=ROSE,ms=8); ax.plot([0.5],[np.sum(tB**2)/n],'o',color=HEAD,ms=8); ax.plot([8/9],[np.sum(tC**2)/9],'D',color=GREEN,ms=8)
ax.text(0.53,np.sum(tA**2)/n+0.02,'swaps',color=ROSE,fontsize=9); ax.text(0.53,np.sum(tB**2)/n-0.07,'cycle',color=HEAD,fontsize=9); ax.annotate('×2 mod 9',xy=(8/9,np.sum(tC**2)/9),xytext=(0.62,1.08),color=GREEN,fontsize=9,arrowprops=dict(arrowstyle='->',color=GREEN,lw=1))
ax.set_ylim(-0.03,1.3); ax.set_xlabel('Hamming distance d_H'); ax.set_ylabel('(1/n) Σ θ²  (squared geodesic, per point)')
fig.tight_layout(); fig.savefig('fig_perm7.svg',bbox_inches='tight'); fig.savefig('fig_perm7.png',dpi=110,bbox_inches='tight'); plt.close(fig)
print('ok')
# ---------- ring (Example 12.2): light both ways round a ring
r=np.load('../checks/ring_rows.npy')
kap=0.05
fig,axs=plt.subplots(1,2,figsize=(11,4.4),gridspec_kw={'width_ratios':[1.25,1]})
ax=axs[0]
ax.loglog(r[:,0]/(2*kap),r[:,2],'o-',color=ROSE,lw=2,ms=5,label='tip: one mode (median)')
ax.loglog(r[:,0]/(2*kap),r[:,3],'s-',color=HEAD,lw=2,ms=5,label='tail: the forward/backward pair (median)')
ax.loglog(r[:,0]/(2*kap),r[:,4],':',color=WARM,lw=1.6,label='Davis–Kahan bound for the pair, ‖E‖/F')
ax.axvline(1,color=MUTED,ls='--',lw=1); ax.text(1.08,0.27,'Ω = 2κ',color=MUTED,fontsize=10)
ax.set_xlabel('rotation against scattering, Ω / 2κ'); ax.set_ylabel('turn under a perturbation (degrees)')
ax.set_title('the tips swing until rotation separates them; the pair holds',fontsize=11.5,color=FG)
ax.legend(frameon=False,fontsize=9.5,labelcolor=FG,loc='center left')
ax=axs[1]; x=np.geomspace(0.05,20,300); th=0.5*np.degrees(np.arctan(1/x))
ax.semilogx(x,th,color=VIO,lw=2.4); ax.axhline(22.5,color=MUTED,ls=':',lw=1); ax.axvline(1,color=MUTED,ls='--',lw=1)
ax.text(0.06,14,'standing waves\n(scattering wins,\ntwo directions act as one)',color=FG,fontsize=9.5)
ax.text(1.5,30,'travelling waves\n(rotation wins,\neach direction readable)',color=FG,fontsize=9.5)
ax.text(1.1,24,'22.5°',color=MUTED,fontsize=10)
ax.set_xlabel('Ω / 2κ'); ax.set_ylabel('mixing angle inside the pair (degrees)'); ax.set_ylim(0,47)
ax.set_title('which tip the pair shows',fontsize=11.5,color=FG)
fig.tight_layout(); fig.savefig('fig_ring.svg',bbox_inches='tight'); fig.savefig('fig_ring.png',dpi=110,bbox_inches='tight'); plt.close(fig)
print('ring ok')
