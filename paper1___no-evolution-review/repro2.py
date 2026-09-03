"""Second pass: proper narrow-bin flux-flux test + luminosity-space bins."""
import numpy as np, json
from scipy.optimize import minimize
from scipy.integrate import quad
rng = np.random.default_rng(7)
rows=[]
for line in open("lusso2020.tsv"):
    if line.startswith("#") or not line.strip(): continue
    p=[c.strip() for c in line.split("\t")]
    try: rows.append([float(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[4])])
    except ValueError: continue
z,xuv,exuv,xx,exx=np.array(rows).T
H0,Om=70.,0.315; c=299792.458
DL=np.array([(1+v)*c/H0*quad(lambda x:1/np.sqrt(Om*(1+x)**3+1-Om),0,v)[0] for v in z])*3.0857e24
Luv=xuv+np.log10(4*np.pi)+2*np.log10(DL); Lx=xx+np.log10(4*np.pi)+2*np.log10(DL)

def nll(t,x,y,ex,ey,xp):
    g,b,dl=t
    if not (0<g<1.5 and -60<b<60 and 0<dl<1): return np.inf
    s2=dl**2+ey**2+(g*ex)**2; r=y-g*(x-xp)-b
    return 0.5*np.sum(r**2/s2+np.log(2*np.pi*s2))
def mle(x,y,ex,ey,xp):
    return minimize(nll,[0.6,np.mean(y-0.6*(x-xp)),0.25],args=(x,y,ex,ey,xp),method="Nelder-Mead",
                    options=dict(xatol=1e-6,fatol=1e-6,maxiter=20000)).x
def mle_boot(x,y,ex,ey,xp,nb=100):
    g=mle(x,y,ex,ey,xp)[0]; bs=[]
    for _ in range(nb):
        j=rng.integers(0,len(x),len(x)); bs.append(mle(x[j],y[j],ex[j],ey[j],xp)[0])
    return g,np.std(bs)
def linfit(zz,g,s):
    W=1/s**2; A=np.vstack([np.ones(len(zz)),zz]).T; cov=np.linalg.inv(A.T@(W[:,None]*A)); p=cov@A.T@(W*g)
    chi2=np.sum(W*(g-p[0]-p[1]*zz)**2); return p, np.sqrt(np.diag(cov)), chi2
def constfit(g,s):
    W=1/s**2; gw=np.sum(W*g)/np.sum(W); return gw,1/np.sqrt(np.sum(W)),np.sum(W*(g-gw)**2)

res={}
# ---- A. narrow bins in log(1+z): width 0.05 dex (=> D_L spread ~0.1 dex), N>=60 ----
print("[A] NARROW BINS flux-flux, width dlog(1+z)=0.05, N>=60")
lz=np.log10(1+z); w=0.05; e=np.arange(lz.min(),lz.max()+w,w)
za,ga,sa,na=[],[],[],[]
for i in range(len(e)-1):
    m=(lz>=e[i])&(lz<e[i+1])
    if m.sum()<60: continue
    xp=np.median(xuv[m]); g,s=mle_boot(xuv[m],xx[m],exuv[m],exx[m],xp,nb=80)
    za.append(np.mean(z[m])); ga.append(g); sa.append(s); na.append(int(m.sum()))
    print(f"   <z>={np.mean(z[m]):.2f} N={m.sum():4d} gamma={g:.3f}+-{s:.3f}")
za,ga,sa=map(np.array,(za,ga,sa))
p,pe,chi2=linfit(za,ga,sa); gw,gwe,chi2c=constfit(ga,sa)
print(f"   linear: g0={p[0]:.3f}+-{pe[0]:.3f}  m={p[1]:.4f}+-{pe[1]:.4f}  chi2/dof={chi2:.1f}/{len(za)-2}")
print(f"   const : gamma={gw:.3f}+-{gwe:.3f}  chi2/dof={chi2c:.1f}/{len(za)-1}")
res["narrow"]=dict(z=za.tolist(),gamma=ga.tolist(),sig=sa.tolist(),n=na,m=[float(p[1]),float(pe[1])],g0=[float(p[0]),float(pe[0])],
                   chi2_lin=float(chi2),gconst=[float(gw),float(gwe)],chi2_const=float(chi2c))
# same but drop z<0.5 where D_L spread is worst even in narrow bins
k=za>0.5; p2,pe2,chi22=linfit(za[k],ga[k],sa[k]); gw2,gwe2,chi2c2=constfit(ga[k],sa[k])
print(f"   z>0.5 only: m={p2[1]:.4f}+-{pe2[1]:.4f} chi2/dof={chi22:.1f}/{k.sum()-2} | const gamma={gw2:.3f}+-{gwe2:.3f} chi2/dof={chi2c2:.1f}/{k.sum()-1}")
res["narrow_z>0.5"]=dict(m=[float(p2[1]),float(pe2[1])],gconst=[float(gw2),float(gwe2)],chi2_const=float(chi2c2),dof=int(k.sum()-1))
# log(1+z) parametrization as in paper's conclusion
X=np.log10(1+za); pk,pke,_=linfit(X,ga,sa)
print(f"   gamma=g0+k log10(1+z): k={pk[1]:.3f}+-{pke[1]:.3f}")
res["k_log1pz"]=[float(pk[1]),float(pke[1])]

# ---- B. 4 quartile bins in LUMINOSITY space (fiducial cosmology; what Table 1 beta~6 implies) ----
print("\n[B] FOUR QUARTILE BINS, luminosity space (H0=70,Om=0.315)")
ed=np.quantile(z,[0,.25,.5,.75,1]); zb,gb,sb=[],[],[]
for i in range(4):
    m=(z>=ed[i])&(z<=ed[i+1]) if i==3 else (z>=ed[i])&(z<ed[i+1])
    t=mle(Luv[m],Lx[m],exuv[m],exx[m],0.0); g,s=mle_boot(Luv[m],Lx[m],exuv[m],exx[m],0.0,nb=80)
    zb.append(np.mean(z[m])); gb.append(g); sb.append(s)
    print(f"   {ed[i]:.2f}-{ed[i+1]:.2f} N={m.sum()} <z>={np.mean(z[m]):.2f} gamma={g:.3f}+-{s:.3f} beta={t[1]:.2f} delta={t[2]:.3f}")
zb,gb,sb=map(np.array,(zb,gb,sb)); p,pe,chi2=linfit(zb,gb,sb)
print(f"   linear m={p[1]:.4f}+-{pe[1]:.4f} chi2/dof={chi2:.1f}/2")
res["lum4"]=dict(z=zb.tolist(),gamma=gb.tolist(),sig=sb.tolist(),m=[float(p[1]),float(pe[1])])

# ---- C. paper's Table-1 bins in luminosity space ----
print("\n[C] PAPER TABLE-1 BINS (0.1-0.9-1.5-2.3-5.0), luminosity space")
for lo,hi in [(0.1,0.9),(0.9,1.5),(1.5,2.3),(2.3,5.0)]:
    m=(z>=lo)&(z<hi); t=mle(Luv[m],Lx[m],exuv[m],exx[m],0.0)
    print(f"   {lo}-{hi} N={m.sum()} gamma={t[0]:.3f} beta={t[1]:.2f} delta={t[2]:.3f}")
json.dump(res,open("repro2_results.json","w"),indent=1)
