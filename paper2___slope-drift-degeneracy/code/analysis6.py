"""Doubt the drift: (1) within-bin D_L variation biases gamma toward 1 at low z. Detrend fluxes to the bin-centre distance
using only the *shape* of D_L(z) inside the bin (H0-independent; weak Om dependence, tested). (2) group-5 & z>0.7."""
import numpy as np, json, emcee
from scipy.optimize import minimize
from scipy.integrate import quad
rng=np.random.default_rng(66)
rows=[]
for line in open("lusso2020.tsv"):
    if line.startswith("#") or not line.strip(): continue
    p=[c.strip() for c in line.split("\t")]
    try: rows.append([float(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[4]),int(p[5])])
    except ValueError: continue
z,x,ex,y,ey,grp=np.array(rows).T; N=len(z); lz=np.log10(1+z); c=299792.458
def logDL(zz,Om):  # up to an additive constant (H0 drops out)
    return np.array([np.log10((1+v)*quad(lambda t:1/np.sqrt(Om*(1+t)**3+1-Om),0,v)[0]) for v in zz])
w=0.03; e=np.arange(lz.min(),lz.max()+w,w); bid=np.full(N,-1); k=0
for i in range(len(e)-1):
    m=(lz>=e[i])&(lz<e[i+1])
    if m.sum()>=25: bid[m]=k; k+=1
use=bid>=0
def prep(mask,Om=None):
    zu,xu,exu,yu,eyu,bu=z[mask],x[mask],ex[mask],y[mask],ey[mask],bid[mask]
    ub=np.unique(bu); remap={b:i for i,b in enumerate(ub)}; bu=np.array([remap[b] for b in bu]); K=len(ub)
    if Om is not None:  # detrend: move every quasar to its bin's median-z distance
        ld=logDL(zu,Om); zc_bin=np.array([np.median(zu[bu==i]) for i in range(K)]); ldc=logDL(zc_bin,Om)
        shift=2*(ld-ldc[bu]); xu=xu+shift; yu=yu+shift    # flux at bin-centre distance: F' = F (D/Dc)^2
    xp=np.array([np.median(xu[bu==i]) for i in range(K)]); u=xu-xp[bu]
    return dict(u=u,y=yu,ex=exu,ey=eyu,b=bu,K=K,zc=zu-1.3,z=zu)
def run(D,label,nstep=3500,burn=1200):
    K=D["K"]
    def lnp(t):
        g0,s,dl=t[0],t[1],t[2]
        if not(0<g0<1.5 and 0<dl<1): return -np.inf
        g=g0+s*D["zc"]; mu=g*D["u"]+t[3:][D["b"]]; s2=dl**2+D["ey"]**2+(g*D["ex"])**2
        return -0.5*np.sum((D["y"]-mu)**2/s2+np.log(2*np.pi*s2))
    p0=np.concatenate([[0.58,-0.04,0.22],[np.median(D["y"][D["b"]==i]) for i in range(K)]])
    r=minimize(lambda t:-lnp(t),p0,method="L-BFGS-B"); r=minimize(lambda t:-lnp(t),r.x,method="Nelder-Mead",options=dict(maxiter=80000))
    nd=len(r.x); nw=2*nd+4; sam=emcee.EnsembleSampler(nw,nd,lnp); sam.run_mcmc(r.x+1e-3*rng.standard_normal((nw,nd)),nstep,progress=False)
    ch=sam.get_chain(discard=burn,flat=True); q=np.percentile(ch[:,:3],[16,50,84],axis=0)
    # per-bin gammas (plain MLE) for the figure
    gb=[]
    for i in range(K):
        m=D["b"]==i
        if m.sum()<60: gb.append(None); continue
        def nb(t):
            g,b,dl=t
            if not(0<g<1.5 and 0<dl<1): return np.inf
            s2=dl**2+D["ey"][m]**2+(g*D["ex"][m])**2; return 0.5*np.sum((D["y"][m]-g*D["u"][m]-b)**2/s2+np.log(2*np.pi*s2))
        t=minimize(nb,[0.6,np.median(D["y"][m]),0.23],method="Nelder-Mead").x; gb.append([float(np.mean(D["z"][m])),float(t[0]),int(m.sum())])
    out=dict(N=int(len(D["y"])),K=K,g0=q[:,0].tolist(),m=q[:,1].tolist(),delta=q[:,2].tolist(),P_m_ge_0=float(np.mean(ch[:,1]>=0)),bins=[b for b in gb if b])
    print(f"{label}: N={out['N']} g0={q[1,0]:.4f}+-{(q[2,0]-q[0,0])/2:.4f} m={q[1,1]:+.4f}+-{(q[2,1]-q[0,1])/2:.4f} delta={q[1,2]:.4f} P(m>=0)={out['P_m_ge_0']:.4f}")
    return out
R={}
R["raw"]=run(prep(use),"raw flux-flux (as before)")
for om in (0.315,0.25,0.40):
    R[f"detrended_Om{om}"]=run(prep(use,Om=om),f"detrended to bin-centre distance, Om={om}")
R["detrended_z>0.7"]=run(prep(use&(z>0.7),Om=0.315),"detrended, z>0.7")
R["detrended_z<0.7"]=run(prep(use&(z<0.7),Om=0.315),"detrended, z<0.7")
R["detrended_group5"]=run(prep(use&(grp==5),Om=0.315),"detrended, group 5")
R["detrended_group5_z>0.7"]=run(prep(use&(grp==5)&(z>0.7),Om=0.315),"detrended, group 5 & z>0.7")
# size of the within-bin distance spread per bin (to explain the effect)
ld=logDL(z[use],0.315); bu=bid[use]
spread=[float(np.std(2*ld[bu==i])) for i in np.unique(bu)]; zb=[float(np.mean(z[use][bu==i])) for i in np.unique(bu)]
print("within-bin std of 2*logD_L per bin (dex):",[f"{a:.2f}:{b:.3f}" for a,b in zip(zb,spread)])
R["within_bin_2logDL_std"]=dict(z=zb,std=spread)
json.dump(R,open("analysis6_results.json","w"),indent=1); print("saved")
