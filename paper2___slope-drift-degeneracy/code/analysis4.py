"""Paper-2 core: (A) binning-free joint flux-space likelihood with per-bin distance nuisance and global gamma(z);
(B) model comparison gamma(z) vs gamma(<L>) vs constant; (C) within-bin curvature vs curvature needed to fake the drift;
(D) truncation-corrected (Eddington-bias-aware) per-bin slopes; (E) photon-index-controlled slopes.
Lusso+2020 table3 only. Nothing here is tuned; all thresholds are stated in the paper."""
import numpy as np, json, emcee
from scipy.optimize import minimize
from scipy.integrate import quad
from scipy import stats
from scipy.special import log_ndtr
rng=np.random.default_rng(11)
rows=[]
for line in open("lusso2020.tsv"):
    if line.startswith("#") or not line.strip(): continue
    p=[c.strip() for c in line.split("\t")]
    try: rows.append([float(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[4]),int(p[5]),float(p[6])])
    except ValueError: continue
z,x,ex,y,ey,grp,gx=np.array(rows).T; N=len(z)
c=299792.458
def DLcm(zz,Om=0.315,H0=70.):
    return np.array([(1+v)*c/H0*quad(lambda t:1/np.sqrt(Om*(1+t)**3+1-Om),0,v)[0] for v in zz])*3.0857e24
def lum(Om): d=DLcm(z,Om); return x+np.log10(4*np.pi)+2*np.log10(d)
Luv=lum(0.315)
lz=np.log10(1+z)
R={}

# ---------- bins for the distance nuisance: w=0.03 in log(1+z), N>=25 ----------
w=0.03; e=np.arange(lz.min(),lz.max()+w,w); bid=np.full(N,-1); k=0; binz=[]; binL=[]
for i in range(len(e)-1):
    m=(lz>=e[i])&(lz<e[i+1])
    if m.sum()>=25: bid[m]=k; binz.append(np.mean(z[m])); binL.append(np.mean(Luv[m])); k+=1
use=bid>=0; K=k; binz=np.array(binz); binL=np.array(binL)
print(f"joint fit: {use.sum()} quasars in {K} bins (w={w}); dropped {N-use.sum()} in sparse bins")
zu,xu,exu,yu,eyu,bu,gxu=z[use],x[use],ex[use],y[use],ey[use],bid[use],gx[use]
xpiv=np.zeros(K)
for i in range(K): xpiv[i]=np.median(xu[bu==i])
u=xu-xpiv[bu]                      # flux offset from own-bin median  (cosmology-free)
zc=zu-1.3                          # z centred
Lc=(binL-30.0)[bu]                 # bin-mean luminosity centred (fiducial cosmology; only used for model ML)
R["joint_setup"]=dict(N=int(use.sum()),K=int(K),w=w,binz=binz.tolist(),binL=binL.tolist())
print(f"  across-bin slope d<logL_UV>/dz = {np.polyfit(binz,binL,1)[0]:.3f} dex per unit z; corr={np.corrcoef(binz,binL)[0,1]:.3f}")
R["dLdz"]=float(np.polyfit(binz,binL,1)[0])

def nll_joint(theta, mode):
    # theta = [g0, slope_param, delta, beta_1..beta_K]
    g0, s, dl = theta[0], theta[1], theta[2]; beta=theta[3:]
    if not (0<g0<1.5 and 0<dl<1): return np.inf
    if mode=="const": g=g0
    elif mode=="z": g=g0+s*zc
    elif mode=="L": g=g0+s*Lc
    mu=g*u+beta[bu]; s2=dl**2+eyu**2+(g*exu)**2
    return 0.5*np.sum((yu-mu)**2/s2+np.log(2*np.pi*s2))
def fit_joint(mode, boot=0):
    p0=np.concatenate([[0.58,0.0,0.23],[np.median(yu[bu==i]) for i in range(K)]])
    f=lambda t: nll_joint(t,mode)
    r=minimize(f,p0,method="L-BFGS-B",options=dict(maxiter=5000))
    r=minimize(f,r.x,method="Nelder-Mead",options=dict(maxiter=60000,xatol=1e-8,fatol=1e-8))
    return r.x, r.fun
print("\n[A/B] joint likelihood, 3 models")
fits={}
for mode in ("const","z","L"):
    th,f=fit_joint(mode); npar=3+K-(1 if mode=="const" else 0)
    aic=2*f+2*npar; bic=2*f+npar*np.log(len(yu))
    fits[mode]=dict(theta=th,nll=f,aic=aic,bic=bic,npar=npar)
    print(f"  {mode:5s}: g0={th[0]:.4f} s={th[1]:+.4f} delta={th[2]:.4f}  -lnL={f:.2f} AIC={aic:.2f} BIC={bic:.2f}")
R["joint"]={m:dict(g0=float(v["theta"][0]),s=float(v["theta"][1]),delta=float(v["theta"][2]),nll=float(v["nll"]),aic=float(v["aic"]),bic=float(v["bic"])) for m,v in fits.items()}
LR=2*(fits["const"]["nll"]-fits["z"]["nll"]); print(f"  LR test const vs z: 2dlnL={LR:.2f} -> p={1-stats.chi2.cdf(LR,1):.4f}  ({np.sqrt(max(LR,0)):.1f} sigma-equiv)")
R["joint"]["LR_const_vs_z"]=[float(LR),float(1-stats.chi2.cdf(LR,1))]
R["joint"]["dAIC_z_minus_L"]=float(fits["z"]["aic"]-fits["L"]["aic"])
# emcee on the z model for honest error on m (marginalising the K intercepts)
print("  emcee on gamma(z) model ...")
th0=fits["z"]["theta"]; nd=len(th0); nw=2*nd+4
p0=th0+1e-3*rng.standard_normal((nw,nd))
sam=emcee.EnsembleSampler(nw,nd,lambda t:-nll_joint(t,"z"))
sam.run_mcmc(p0,4000,progress=False); ch=sam.get_chain(discard=1500,flat=True)
q=np.percentile(ch[:,:3],[16,50,84],axis=0)
print(f"  gamma0(z=1.3)={q[1,0]:.4f} +{q[2,0]-q[1,0]:.4f} -{q[1,0]-q[0,0]:.4f} | m={q[1,1]:.4f} +{q[2,1]-q[1,1]:.4f} -{q[1,1]-q[0,1]:.4f} | delta={q[1,2]:.4f}")
print(f"  P(m>=0) = {np.mean(ch[:,1]>=0):.4f}; 95% upper bound on |m| if m<0: {np.percentile(-ch[:,1],95):.4f}")
R["joint"]["mcmc_z"]=dict(g0=q[:,0].tolist(),m=q[:,1].tolist(),delta=q[:,2].tolist(),P_m_ge_0=float(np.mean(ch[:,1]>=0)),m_abs_95=float(np.percentile(-ch[:,1],95)))
# bootstrap over quasars (robustness of m to sampling; keeps bins)
bs=[]
for _ in range(60):
    j=rng.integers(0,len(yu),len(yu));
    def nb(t):
        g=t[0]+t[1]*zc[j]; mu=g*u[j]+t[3:][bu[j]]; s2=t[2]**2+eyu[j]**2+(g*exu[j])**2
        return np.inf if not(0<t[0]<1.5 and 0<t[2]<1) else 0.5*np.sum((yu[j]-mu)**2/s2+np.log(2*np.pi*s2))
    r=minimize(nb,th0,method="L-BFGS-B"); bs.append(r.x[1])
print(f"  bootstrap sigma_m = {np.std(bs):.4f}")
R["joint"]["boot_sigma_m"]=float(np.std(bs))
# same on group 5 only
mask5=(grp[use]==5);
def nll5(t):
    g=t[0]+t[1]*zc[mask5]; mu=g*u[mask5]+t[3:][bu[mask5]]; s2=t[2]**2+eyu[mask5]**2+(g*exu[mask5])**2
    return np.inf if not(0<t[0]<1.5 and 0<t[2]<1) else 0.5*np.sum((yu[mask5]-mu)**2/s2+np.log(2*np.pi*s2))
r5=minimize(nll5,th0,method="L-BFGS-B"); r5=minimize(nll5,r5.x,method="Nelder-Mead",options=dict(maxiter=60000))
bs5=[]
for _ in range(60):
    idx=np.where(mask5)[0]; j=rng.choice(idx,len(idx))
    def nb(t):
        g=t[0]+t[1]*zc[j]; mu=g*u[j]+t[3:][bu[j]]; s2=t[2]**2+eyu[j]**2+(g*exu[j])**2
        return np.inf if not(0<t[0]<1.5 and 0<t[2]<1) else 0.5*np.sum((yu[j]-mu)**2/s2+np.log(2*np.pi*s2))
    bs5.append(minimize(nb,r5.x,method="L-BFGS-B").x[1])
print(f"  group-5 only (N={mask5.sum()}): m={r5.x[1]:.4f} +- {np.std(bs5):.4f}, g0={r5.x[0]:.4f}, delta={r5.x[2]:.4f}")
R["joint"]["group5"]=dict(N=int(mask5.sum()),m=[float(r5.x[1]),float(np.std(bs5))],g0=float(r5.x[0]),delta=float(r5.x[2]))

# ---------- (C) within-bin curvature vs needed ----------
print("\n[C] within-bin quadratic curvature (cosmology-free) vs curvature needed to mimic the drift")
def nll_q(t,xx,yy,exx,eyy):
    g,b,dl,qq=t
    if not(0<g<1.5 and 0<dl<1 and -1<qq<1): return np.inf
    s2=dl**2+eyy**2+((g+2*qq*xx)*exx)**2; r=yy-g*xx-qq*xx*xx-b
    return 0.5*np.sum(r**2/s2+np.log(2*np.pi*s2))
qs=[];qe=[];qz=[]
for i in range(K):
    m=bu==i
    if m.sum()<100: continue
    xx=u[m]; r=minimize(nll_q,[0.6,np.median(yu[m]),0.23,0.0],args=(xx,yu[m],exu[m],eyu[m]),method="Nelder-Mead",options=dict(maxiter=40000))
    bsq=[]
    for _ in range(60):
        j=rng.integers(0,m.sum(),m.sum()); bsq.append(minimize(nll_q,r.x,args=(xx[j],yu[m][j],exu[m][j],eyu[m][j]),method="Nelder-Mead",options=dict(maxiter=20000)).x[3])
    qs.append(r.x[3]); qe.append(np.std(bsq)); qz.append(binz[i])
qs,qe=np.array(qs),np.array(qe); W=1/qe**2; qw=np.sum(W*qs)/np.sum(W); qwe=1/np.sqrt(np.sum(W))
kappa_within=2*qw; kappa_within_e=2*qwe
m_z=R["joint"]["mcmc_z"]["m"][1]; kappa_needed=m_z/R["dLdz"]
print(f"  {len(qs)} bins with N>=100: weighted q = {qw:+.4f} +- {qwe:.4f}  -> d(gamma)/d(logL) within bins = 2q = {kappa_within:+.4f} +- {kappa_within_e:.4f}")
print(f"  to mimic m={m_z:+.4f}/z with d<L>/dz={R['dLdz']:.3f}: need d(gamma)/d(logL) = {kappa_needed:+.4f}")
print(f"  tension between measured and needed curvature: {(kappa_within-kappa_needed)/kappa_within_e:.1f} sigma")
R["curvature"]=dict(q_bins=qs.tolist(),q_err=qe.tolist(),z=qz,q_weighted=[float(qw),float(qwe)],kappa_within=[float(kappa_within),float(kappa_within_e)],kappa_needed=float(kappa_needed),
                    tension_sigma=float((kappa_within-kappa_needed)/kappa_within_e))

# ---------- (D) truncation-corrected per-bin slopes ----------
print("\n[D] truncated-Gaussian likelihood: y > ylim_i (per-bin X-ray flux floor)")
def nll_trunc(t,xx,yy,exx,eyy,ylim):
    g,b,dl=t
    if not(0<g<1.5 and 0<dl<1): return np.inf
    s2=dl**2+eyy**2+(g*exx)**2; mu=g*xx+b
    # p(y|x, y>ylim) = N(y;mu,s)/ (1-Phi((ylim-mu)/s))
    lognorm=log_ndtr((mu-ylim)/np.sqrt(s2))
    return 0.5*np.sum((yy-mu)**2/s2+np.log(2*np.pi*s2))+np.sum(lognorm)
def nll_plain(t,xx,yy,exx,eyy):
    g,b,dl=t
    if not(0<g<1.5 and 0<dl<1): return np.inf
    s2=dl**2+eyy**2+(g*exx)**2; r=yy-g*xx-b
    return 0.5*np.sum(r**2/s2+np.log(2*np.pi*s2))
def wlin(zz,g,s):
    W=1/s**2; A=np.vstack([np.ones(len(zz)),zz]).T; cov=np.linalg.inv(A.T@(W[:,None]*A)); p=cov@A.T@(W*g); return p[1],np.sqrt(cov[1,1])
R["trunc"]={}
for pct in (0.0,1.0,3.0):
    gp,gt,se,zz=[],[],[],[]
    for i in range(K):
        m=bu==i
        if m.sum()<60: continue
        xx,yy,exx,eyy=u[m],yu[m],exu[m],eyu[m]; ylim=np.percentile(yy,pct)-1e-6
        keep=yy>ylim; xx,yy,exx,eyy=xx[keep],yy[keep],exx[keep],eyy[keep]
        p=minimize(nll_plain,[0.6,np.median(yy),0.23],args=(xx,yy,exx,eyy),method="Nelder-Mead",options=dict(maxiter=30000)).x
        t=minimize(nll_trunc,p,args=(xx,yy,exx,eyy,ylim),method="Nelder-Mead",options=dict(maxiter=30000)).x
        bsg=[]
        for _ in range(40):
            j=rng.integers(0,len(yy),len(yy)); bsg.append(minimize(nll_trunc,t,args=(xx[j],yy[j],exx[j],eyy[j],ylim),method="Nelder-Mead",options=dict(maxiter=15000)).x[0])
        gp.append(p[0]); gt.append(t[0]); se.append(np.std(bsg)); zz.append(binz[i])
    gp,gt,se,zz=map(np.array,(gp,gt,se,zz))
    mp,mpe=wlin(zz,gp,se); mt,mte=wlin(zz,gt,se)
    print(f"  floor=p{pct:g}: mean(gamma_trunc - gamma_plain)={np.mean(gt-gp):+.4f}; m_plain={mp:+.4f}+-{mpe:.4f}  m_trunc={mt:+.4f}+-{mte:.4f}")
    R["trunc"][f"p{pct:g}"]=dict(z=zz.tolist(),g_plain=gp.tolist(),g_trunc=gt.tolist(),sig=se.tolist(),m_plain=[float(mp),float(mpe)],m_trunc=[float(mt),float(mte)])
# how hard does the floor bite? fraction within 0.3 dex of the floor per bin
frac=[np.mean(yu[bu==i]<np.min(yu[bu==i])+0.3) for i in range(K)]
print("  fraction of bin within 0.3 dex of its X-ray floor, low->high z:",np.round(frac,2).tolist())
R["trunc"]["frac_near_floor"]=[float(f) for f in frac]

# ---------- (E) photon-index-controlled slopes ----------
print("\n[E] gamma(z) with a Gamma_X term: logF_X = gamma u + c (Gamma_X - 2.175) + beta_i")
gmed=np.nanmedian(gxu); dG=gxu-gmed
def nll_G(t):
    g0,s,dl,cG=t[0],t[1],t[2],t[3]; beta=t[4:]
    if not(0<g0<1.5 and 0<dl<1): return np.inf
    g=g0+s*zc; mu=g*u+cG*dG+beta[bu]; s2=dl**2+eyu**2+(g*exu)**2
    return 0.5*np.sum((yu-mu)**2/s2+np.log(2*np.pi*s2))
p0=np.concatenate([[th0[0],th0[1],th0[2],0.0],th0[3:]])
rG=minimize(nll_G,p0,method="L-BFGS-B"); rG=minimize(nll_G,rG.x,method="Nelder-Mead",options=dict(maxiter=80000))
bsG=[]
for _ in range(40):
    j=rng.integers(0,len(yu),len(yu))
    def nb(t):
        g=t[0]+t[1]*zc[j]; mu=g*u[j]+t[3]*dG[j]+t[4:][bu[j]]; s2=t[2]**2+eyu[j]**2+(g*exu[j])**2
        return np.inf if not(0<t[0]<1.5 and 0<t[2]<1) else 0.5*np.sum((yu[j]-mu)**2/s2+np.log(2*np.pi*s2))
    bsG.append(minimize(nb,rG.x,method="L-BFGS-B").x[:4])
bsG=np.array(bsG)
print(f"  gamma0={rG.x[0]:.4f}  m={rG.x[1]:+.4f}+-{bsG[:,1].std():.4f}  c_Gamma={rG.x[3]:+.4f}+-{bsG[:,3].std():.4f} dex/unit  delta={rG.x[2]:.4f} (was {th0[2]:.4f})")
R["photon_index"]=dict(g0=float(rG.x[0]),m=[float(rG.x[1]),float(bsG[:,1].std())],c=[float(rG.x[3]),float(bsG[:,3].std())],delta=float(rG.x[2]),delta_before=float(th0[2]))
rho,pv=stats.spearmanr(zu,gxu); R["photon_index"]["spearman_z_Gx"]=[float(rho),float(pv)]
# within-bin correlation of Gamma_X with F_UV (would let Gamma_X leak into gamma)
cc=[np.corrcoef(u[bu==i],dG[bu==i])[0,1] for i in range(K) if (bu==i).sum()>=60]
print(f"  within-bin corr(F_UV, Gamma_X): mean={np.mean(cc):+.3f}, range {min(cc):+.2f}..{max(cc):+.2f}")
R["photon_index"]["within_bin_corr_FUV_Gx"]=[float(np.mean(cc)),float(min(cc)),float(max(cc))]

json.dump(R,open("analysis4_results.json","w"),indent=1); print("\nsaved analysis4_results.json")
