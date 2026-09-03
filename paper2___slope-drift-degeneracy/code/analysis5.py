"""Honest MCMC errors for the two sub-fits whose bootstrap errors looked too small (group-5 only; Gamma_X-controlled)."""
import numpy as np, json, emcee
from scipy.optimize import minimize
rng=np.random.default_rng(5)
rows=[]
for line in open("lusso2020.tsv"):
    if line.startswith("#") or not line.strip(): continue
    p=[c.strip() for c in line.split("\t")]
    try: rows.append([float(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[4]),int(p[5]),float(p[6])])
    except ValueError: continue
z,x,ex,y,ey,grp,gx=np.array(rows).T; N=len(z); lz=np.log10(1+z)
w=0.03; e=np.arange(lz.min(),lz.max()+w,w); bid=np.full(N,-1); k=0
for i in range(len(e)-1):
    m=(lz>=e[i])&(lz<e[i+1])
    if m.sum()>=25: bid[m]=k; k+=1
use=bid>=0; K=k
def prep(mask):
    zu,xu,exu,yu,eyu,bu,gxu=z[mask],x[mask],ex[mask],y[mask],ey[mask],bid[mask],gx[mask]
    # re-index bins present
    ub=np.unique(bu); remap={b:i for i,b in enumerate(ub)}; bu=np.array([remap[b] for b in bu]); Kl=len(ub)
    xp=np.array([np.median(xu[bu==i]) for i in range(Kl)]); u=xu-xp[bu]
    return dict(u=u,y=yu,ex=exu,ey=eyu,b=bu,K=Kl,zc=zu-1.3,dG=gxu-2.175)
def run(D,label,withG=False,nstep=4000,burn=1500):
    K=D["K"]
    def lnp(t):
        g0,s,dl=t[0],t[1],t[2]
        if not(0<g0<1.5 and 0<dl<1): return -np.inf
        cG=t[3] if withG else 0.0; beta=t[4:] if withG else t[3:]
        g=g0+s*D["zc"]; mu=g*D["u"]+cG*D["dG"]+beta[D["b"]]; s2=dl**2+D["ey"]**2+(g*D["ex"])**2
        return -0.5*np.sum((D["y"]-mu)**2/s2+np.log(2*np.pi*s2))
    p0=np.concatenate([[0.58,-0.04,0.22],[0.0] if withG else [],[np.median(D["y"][D["b"]==i]) for i in range(K)]])
    r=minimize(lambda t:-lnp(t),p0,method="L-BFGS-B"); r=minimize(lambda t:-lnp(t),r.x,method="Nelder-Mead",options=dict(maxiter=80000))
    nd=len(r.x); nw=2*nd+4; sam=emcee.EnsembleSampler(nw,nd,lnp)
    sam.run_mcmc(r.x+1e-3*rng.standard_normal((nw,nd)),nstep,progress=False)
    ch=sam.get_chain(discard=burn,flat=True); q=np.percentile(ch[:,:4 if withG else 3],[16,50,84],axis=0)
    out=dict(N=int(len(D["y"])),K=int(K),g0=q[:,0].tolist(),m=q[:,1].tolist(),delta=q[:,2].tolist(),P_m_ge_0=float(np.mean(ch[:,1]>=0)))
    if withG: out["c"]=q[:,3].tolist()
    print(f"{label}: N={out['N']} K={K} g0={q[1,0]:.4f}+-{(q[2,0]-q[0,0])/2:.4f} m={q[1,1]:+.4f}+-{(q[2,1]-q[0,1])/2:.4f} delta={q[1,2]:.4f} P(m>=0)={out['P_m_ge_0']:.4f}"+(f" c={q[1,3]:+.4f}+-{(q[2,3]-q[0,3])/2:.4f}" if withG else ""))
    return out
R={}
R["group5"]=run(prep(use&(grp==5)),"group5 z-model")
R["photon_index"]=run(prep(use),"all, z-model + Gamma_X term",withG=True)
R["z_gt_0p7"]=run(prep(use&(z>0.7)),"z>0.7 only (Risaliti+26 range)")
json.dump(R,open("analysis5_results.json","w"),indent=1); print("saved")
