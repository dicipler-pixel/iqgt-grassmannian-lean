# SCRIPT: IQGT-FIG-BAND
# Zero-temperature band of Theorem 5.3 for three levels: as the drive shifts its weight from the first to the
# second excited level, zeta_rot/g moves across [min 2 c_n D_n, max 2 c_n D_n]. Computed from the closed form.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
BG='#0c1220'; FG='#b9dcff'; HEAD='#4cc3ff'; MUTED='#7088a8'; LINE='#1a2740'; WARM='#e0b85c'
plt.rcParams.update({'svg.fonttype':'none','font.family':'serif','font.serif':['Newsreader','DejaVu Serif'],'font.size':11,
 'axes.facecolor':BG,'figure.facecolor':BG,'savefig.facecolor':BG,'text.color':FG,'axes.labelcolor':FG,'axes.edgecolor':MUTED,
 'xtick.color':MUTED,'ytick.color':MUTED,'axes.grid':True,'grid.color':LINE,'grid.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False})
D=np.array([1.0,2.6]); fig,ax=plt.subplots(figsize=(7.2,4.0))
for tau,kp,col,lab in ((0.5,0,HEAD,'pure relaxation, τ = 0.5'),(0.5,1,WARM,'relaxation with precession, τ = 0.5')):
    c=tau/(1+kp**2*tau**2*D**2); w=2*c*D; s=np.linspace(0,1,201)
    g1=(1-s)/D[0]**2; g2=s/D[1]**2          # drive weight moved from level 1 to level 2
    ratio=(w[0]*g1+w[1]*g2)/(g1+g2)
    ax.plot(s,ratio,color=col,lw=2.2,label=lab)
    ax.axhline(w.min(),color=col,lw=0.9,ls='--'); ax.axhline(w.max(),color=col,lw=0.9,ls='--')
ax.set_xlabel('share of the drive coupling the ground state to the second excited level')
ax.set_ylabel('friction per unit metric,  ζ_rot / g')
ax.legend(fontsize=9.5,frameon=False,loc='upper left',bbox_to_anchor=(0.0,0.93))
ax.text(0.02,1.06,'dashed: the band ends  min 2c(Δₙ)Δₙ  and  max 2c(Δₙ)Δₙ',color=MUTED,fontsize=9,transform=ax.transAxes)
fig.savefig('fig_band.svg',bbox_inches='tight'); fig.savefig('fig_band.png',dpi=110,bbox_inches='tight')
