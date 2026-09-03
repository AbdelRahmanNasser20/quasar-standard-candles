"""Paper-2 analyses: is the apparent gamma(z) drift cosmic evolution, luminosity-dependence (curvature),
photon-index dependence, or sample heterogeneity?  All on Lusso+2020 table3 (VizieR J/A+A/642/A150).
Every number printed here goes into the paper; nothing else does.
"""
import numpy as np, json
from scipy.optimize import minimize
from scipy.integrate import quad
from scipy import stats
rng = np.random.default_rng(2026)

rows=[]
for line in open("lusso2020.tsv"):
    if line.startswith("#") or not line.strip(): continue
    p=[c.strip() for c in line.split("\t")]
    try: rows.append([float(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[4]),int(p[5]),float(p[6]) if p[6] else np.nan])
    except ValueError: continue
z,xuv,exuv,xx,exx,grp,gx=np.array(rows).T
N=len(z)
H0,Om=70.,0.315; c=299792.458
def DLcm(zz,Om=Om,H0=H0):
    return np.array([(1+v)*c/H0*quad(lambda x:1/np.sqrt(Om*(1+x)**3+1-Om),0,v)[0] for v in zz])*3.0857e24
DL=DLcm(z); Luv=xuv+np.log10(4*np.pi)+2*np.log10(DL); Lx=xx+np.log10(4*np.pi)+2*np.log10(DL)

def nll(t,x,y,ex,ey,xp,quad_=False):
    if quad_: g,b,dl,q=t
    else: g,b,dl=t; q=0.
    if not (0<g<1.5 and -60<b<60 and 0<dl<1 and -1<q<1): return np.inf
    u=x-xp; s2=dl**2+ey**2+((g+2*q*u)*ex)**2; r=y-g*u-q*u*u-b
    return 0.5*np.sum(r**2/s2+np.log(2*np.pi*s2))
def mle(x,y,ex,ey,xp,quad_=False):
    p0=[0.6,np.mean(y-0.6*(x-xp)),0.25]+([0.]if quad_ else [])
    return minimize(nll,p0,args=(x,y,ex,ey,xp,quad_),method="Nelder-Mead",options=dict(xatol=1e-7,fatol=1e-7,maxiter=40000)).x
def boot(x,y,ex,ey,xp,nb=150,quad_=False,idx=0):
    v=[]
    for _ in range(nb):
        j=rng.integers(0,len(x),len(x)); v.append(mle(x[j],y[j],ex[j],ey[j],xp,quad_)[idx])
    return np.std(v)
def linfit(zz,g,s):
    W=1/s**2; A=np.vstack([np.ones(len(zz)),zz]).T; cov=np.linalg.inv(A.T@(W[:,None]*A)); p=cov@A.T@(W*g)
    return p,np.sqrt(np.diag(cov)),float(np.sum(W*(g-p[0]-p[1]*zz)**2))
def constfit(g,s):
    W=1/s**2; gw=np.sum(W*g)/np.sum(W); return gw,1/np.sqrt(np.sum(W)),float(np.sum(W*(g-gw)**2))

def narrow_bins(mask, w=0.05, nmin=60, nb=100):
    lz=np.log10(1+z); e=np.arange(lz[mask].min(),lz[mask].max()+w,w); out=[]
    for i in range(len(e)-1):
        m=mask&(lz>=e[i])&(lz<e[i+1])
        if m.sum()<nmin: continue
        xp=np.median(xuv[m]); t=mle(xuv[m],xx[m],exuv[m],exx[m],xp); s=boot(xuv[m],xx[m],exuv[m],exx[m],xp,nb)
        out.append(dict(z=float(np.mean(z[m])),n=int(m.sum()),gamma=float(t[0]),sig=float(s),delta=float(t[2]),
                        luv_mean=float(np.mean(Luv[m])),luv_std=float(np.std(Luv[m])),gx_mean=float(np.nanmean(gx[m]))))
    return out
def summarize(bins,label):
    zz=np.array([b["z"] for b in bins]); g=np.array([b["gamma"] for b in bins]); s=np.array([b["sig"] for b in bins])
    p,pe,chi2=linfit(zz,g,s); gw,gwe,chi2c=constfit(g,s)
    print(f"  {label}: nbins={len(zz)} <gamma>={gw:.3f}+-{gwe:.3f} chi2c/dof={chi2c:.1f}/{len(zz)-1} (p={1-stats.chi2.cdf(chi2c,len(zz)-1):.3f}) | m={p[1]:.4f}+-{pe[1]:.4f} ({p[1]/pe[1]:.1f}sig)")
    return dict(label=label,bins=bins,gconst=[float(gw),float(gwe)],chi2_const=chi2c,dof_const=int(len(zz)-1),
                p_const=float(1-stats.chi2.cdf(chi2c,len(zz)-1)),m=[float(p[1]),float(pe[1])],g0=[float(p[0]),float(pe[0])],chi2_lin=chi2)

R={}
allm=np.ones(N,bool)
print(f"N={N}; groups:",dict(zip(*np.unique(grp.astype(int),return_counts=True))))

# ---- T1. baseline + robustness to bin width ----
print("\n[T1] baseline narrow-bin drift; robustness to bin width")
R["T1"]={}
for w in (0.03,0.04,0.05,0.07,0.10):
    b=narrow_bins(allm,w=w,nb=80); R["T1"][str(w)]=summarize(b,f"w={w}")
base=R["T1"]["0.05"]

# ---- T2. Spearman of per-bin gamma with z vs with <logL_UV> (they are collinear in a flux-limited sample) ----
print("\n[T2] collinearity: per-bin <logL_UV> vs z")
bz=np.array([b["z"] for b in base["bins"]]); bl=np.array([b["luv_mean"] for b in base["bins"]]); bg=np.array([b["gamma"] for b in base["bins"]])
print(f"  corr(z,<logL_UV>) across bins = {np.corrcoef(bz,bl)[0,1]:.3f}; corr(gamma,z)={np.corrcoef(bg,bz)[0,1]:.3f}; corr(gamma,<L>)={np.corrcoef(bg,bl)[0,1]:.3f}")
R["T2"]=dict(z=bz.tolist(),luv=bl.tolist(),gamma=bg.tolist(),corr_z_L=float(np.corrcoef(bz,bl)[0,1]))

# ---- T3. FIXED-LUMINOSITY-WINDOW test: same narrow z bins, only quasars in a common L_UV window ----
print("\n[T3] fixed L_UV window (needs fiducial cosmology only to define the window)")
R["T3"]={}
for lo,hi in [(30.0,31.0),(29.8,31.2),(30.2,30.8)]:
    wm=(Luv>=lo)&(Luv<hi); print(f"  window {lo}-{hi}: N={wm.sum()}")
    b=narrow_bins(wm,w=0.07,nmin=50,nb=80); R["T3"][f"{lo}-{hi}"]=summarize(b,f"window {lo}-{hi}")

# ---- T4. curvature in the full-sample luminosity relation ----
print("\n[T4] quadratic term q in logLx = g u + q u^2 + b, u=logL_UV-30 (fiducial cosmology; robustness to Om)")
R["T4"]={}
for om in (0.25,0.315,0.40):
    DLo=DLcm(z,Om=om); Lu=xuv+np.log10(4*np.pi)+2*np.log10(DLo); Lxo=xx+np.log10(4*np.pi)+2*np.log10(DLo)
    t=mle(Lu,Lxo,exuv,exx,30.0,quad_=True); sq=boot(Lu,Lxo,exuv,exx,30.0,nb=60,quad_=True,idx=3)
    print(f"  Om={om}: gamma@30={t[0]:.4f}  q={t[3]:.4f}+-{sq:.4f} ({t[3]/sq:.1f}sig)  delta={t[2]:.3f}")
    R["T4"][str(om)]=dict(gamma=float(t[0]),q=[float(t[3]),float(sq)],delta=float(t[2]))
# curvature in flux-flux inside the two best-populated narrow bins (cosmology-free)
print("  cosmology-free curvature in individual narrow bins:")
lz=np.log10(1+z); R["T4"]["flux_bins"]=[]
for lo,hi in [(0.28,0.33),(0.33,0.38),(0.38,0.43),(0.43,0.48)]:
    m=(lz>=lo)&(lz<hi); xp=np.median(xuv[m]); t=mle(xuv[m],xx[m],exuv[m],exx[m],xp,quad_=True); sq=boot(xuv[m],xx[m],exuv[m],exx[m],xp,nb=60,quad_=True,idx=3)
    print(f"    <z>={np.mean(z[m]):.2f} N={m.sum()} gamma={t[0]:.3f} q={t[3]:.3f}+-{sq:.3f}")
    R["T4"]["flux_bins"].append(dict(z=float(np.mean(z[m])),n=int(m.sum()),gamma=float(t[0]),q=[float(t[3]),float(sq)]))

# ---- T5. photon-index dependence ----
print("\n[T5] X-ray photon index: does gamma depend on Gamma_X, and does <Gamma_X> drift with z?")
ok=np.isfinite(gx); med=np.nanmedian(gx)
print(f"  Gamma_X available for {ok.sum()}; median={med:.3f}")
rho,pv=stats.spearmanr(z[ok],gx[ok]); print(f"  Spearman(z,Gamma_X)={rho:.3f} p={pv:.2e}")
R["T5"]=dict(median_gx=float(med),spearman_z_gx=[float(rho),float(pv)])
for lab,m in (("soft Gamma_X>=median",ok&(gx>=med)),("hard Gamma_X<median",ok&(gx<med))):
    b=narrow_bins(m,w=0.07,nmin=50,nb=80); R["T5"][lab]=summarize(b,lab)
# gamma vs Gamma_X at fixed z: residual test. Fit baseline in each narrow bin, regress residual on (Gamma_X - median)
print("  residual slope d(resid)/d(Gamma_X) pooled over narrow bins:")
res=[];gxr=[]
for i,bn in enumerate(base["bins"]): pass
e=np.arange(lz.min(),lz.max()+0.05,0.05)
for i in range(len(e)-1):
    m=(lz>=e[i])&(lz<e[i+1])&ok
    if m.sum()<60: continue
    xp=np.median(xuv[m]); t=mle(xuv[m],xx[m],exuv[m],exx[m],xp)
    res.append(xx[m]-t[0]*(xuv[m]-xp)-t[1]); gxr.append(gx[m]-med)
res=np.concatenate(res); gxr=np.concatenate(gxr); sl=stats.linregress(gxr,res)
print(f"    slope={sl.slope:.4f}+-{sl.stderr:.4f} dex per unit Gamma_X (p={sl.pvalue:.2e}), N={len(res)}")
R["T5"]["resid_vs_gx"]=[float(sl.slope),float(sl.stderr),float(sl.pvalue)]

# ---- T6. sample heterogeneity: main XMM group only ----
print("\n[T6] homogeneity: group 5 (SDSS-4XMM main) only vs all")
g5=grp==5; print(f"  group5 N={g5.sum()}, z range {z[g5].min():.2f}-{z[g5].max():.2f}")
b=narrow_bins(g5,w=0.05,nmin=60,nb=80); R["T6"]=summarize(b,"group5 only")
g56=(grp==5)|(grp==6); b=narrow_bins(g56,w=0.05,nmin=60,nb=80); R["T6b"]=summarize(b,"groups 5+6")

# ---- T7. what the allowed drift does to distances: mag bias per bin from slope error times luminosity lever arm ----
print("\n[T7] distance-modulus bias budget: dmu = 5*dgamma*dB/(2(1-gamma))")
gam0=base["gconst"][0]; Lref=np.mean(Luv[(z>0.4)&(z<0.7)])
out=[]
for bn in base["bins"]:
    dB=bn["luv_mean"]-Lref
    for lab,dg in (("m_best",base["m"][0]*bn["z"]),("m_2sig",(base["m"][0]-2*base["m"][1])*bn["z"])):
        dmu=5*dg*dB/(2*(1-gam0))
        out.append(dict(z=bn["z"],dB=dB,case=lab,dgamma=dg,dmu=dmu,fracDL=10**(dmu/5)-1))
for o in out:
    if o["case"]=="m_best": print(f"  z={o['z']:.2f} dB={o['dB']:+.2f} best dmu={o['dmu']:+.3f} mag ({100*o['fracDL']:+.1f}% D_L)")
R["T7"]=dict(gamma0=gam0,Lref=float(Lref),rows=out)
# requirement: sigma_m needed for 2% H0 at z~2 with the observed lever arm
dB2=[o["dB"] for o in out if abs(o["z"]-2.0)<0.15][0]
sig_m_req=(0.02/np.log(10)*5)*(2*(1-gam0))/(5*2.0*dB2)  # dmu=0.02*5/ln10 -> dgamma -> m
print(f"  lever arm at z~2: dB={dB2:.2f} dex. sigma_m for 2% D_L at z=2: {sig_m_req:.4f}")
R["T7"]["sigma_m_for_2pct_at_z2"]=float(sig_m_req)

json.dump(R,open("analysis3_results.json","w"),indent=1); print("\nsaved analysis3_results.json")
