# SCRIPT: IQGT-SPIN-J-KAHLER
# Claim: for every spin j, the lowest band of d.S (a spin-j coherent state, the one-particle lowest Landau level
# on the sphere with 2j flux quanta) has g = (j/2) h and |Omega| = j sin(theta) against the round metric h, so it
# meets the curvature bound (6.1) and the trace condition (6.4) with equality, and its Chern number is 2j.
# Kill: any j in 1/2..4 misses equality by more than 1e-6, or the Chern number is not 2j.
# Opposite (control): adding 0.8 Sz^2 must lift the trace condition strictly for j >= 1.
import numpy as np
def spin(j):
    m=np.arange(j,-j-1,-1); n=len(m); Sz=np.diag(m)
    Sp=np.zeros((n,n));
    for k in range(1,n): Sp[k-1,k]=np.sqrt(j*(j+1)-m[k]*(m[k]+1))
    Sx=(Sp+Sp.T)/2; Sy=(Sp-Sp.T)/(2j); return Sx,Sy,Sz
def proj0(H): w,U=np.linalg.eigh(H); return np.outer(U[:,0],U[:,0].conj())
def geom(Hf,th,ph,h=1e-5):
    P=proj0(Hf(th,ph)); Pt=(proj0(Hf(th+h,ph))-proj0(Hf(th-h,ph)))/(2*h); Pp=(proj0(Hf(th,ph+h))-proj0(Hf(th,ph-h)))/(2*h)
    g=lambda A,B:0.5*np.trace(A@B).real
    return g(Pt,Pt),g(Pp,Pp),g(Pt,Pp),(1j*np.trace(P@(Pt@Pp-Pp@Pt))).real
ok=True
for j in (0.5,1,1.5,2,3,4):
    Sx,Sy,Sz=spin(j)
    for extra in (0.0,0.8):
        Hf=lambda th,ph:np.sin(th)*np.cos(ph)*Sx+np.sin(th)*np.sin(ph)*Sy+np.cos(th)*Sz+extra*Sz@Sz
        dev=[];cb=[];gr=[]
        for th in np.linspace(0.3,2.8,7):
            for ph in (0.2,1.9,4.0):
                gtt,gpp,gtp,Om=geom(Hf,th,ph); s=np.sin(th)
                dev.append(gtt+gpp/s**2-abs(Om)/s); cb.append(4*(gtt*gpp-gtp**2)-Om**2); gr.append((gtt,Om/s))
        if extra==0:
            # Chern number by integrating Omega over the sphere
            ths=np.linspace(1e-3,np.pi-1e-3,121); C=np.trapezoid([geom(Hf,t,0.7)[3] for t in ths],ths)*2*np.pi/(2*np.pi)
            good=max(abs(np.array(dev)))<1e-6 and max(abs(np.array(cb)))<1e-6 and abs(abs(C)-2*j)<1e-3 and abs(gr[0][0]-j/2)<1e-6
            ok&=good; print(f'j={j}: trace gap max {max(abs(np.array(dev))):.1e}, curvature deficit max {max(abs(np.array(cb))):.1e}, g_tt={gr[0][0]:.6f} (j/2={j/2}), |Omega|/sin={abs(gr[0][1]):.6f}, Chern {C:.4f}', 'OK' if good else 'FAIL')
        elif j>=1:
            lift=min(dev); ok&=lift>1e-3; print(f'   with 0.8 Sz^2: trace condition lifted, min {lift:.4f}', 'OK' if lift>1e-3 else 'FAIL')
print('RESULT','PASS' if ok else 'FAIL')
