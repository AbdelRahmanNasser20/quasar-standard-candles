"""Reproduce the core claims of 'No Evidence for Redshift Evolution of the LX-LUV relation'
on the public Lusso et al. 2020 catalog (VizieR J/A+A/642/A150, table3).
"""
import numpy as np, emcee, json, sys
from scipy.optimize import minimize
from scipy.integrate import quad

rng = np.random.default_rng(42)

# ---------- load ----------
rows = []
with open("lusso2020.tsv") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = [c.strip() for c in line.split("\t")]
        if p[0] in ("z", "") or p[0].startswith("-"):
            continue
        try:
            rows.append([float(p[0]), float(p[1]), float(p[2]), float(p[3]), float(p[4]), int(p[5])])
        except ValueError:
            continue
d = np.array(rows)
z, xuv, exuv, xx, exx, grp = d.T
N = len(z)
print(f"N={N}  z: min={z.min():.3f} median={np.median(z):.3f} max={z.max():.3f}")
print(f"N with 0.1<=z<=5.0: {np.sum((z>=0.1)&(z<=5.0))}")
print(f"groups: {dict(zip(*np.unique(grp.astype(int), return_counts=True)))}")

# ---------- cosmology (only for luminosity-space fit / plotting) ----------
H0, Om = 70.0, 0.315
c_kms = 299792.458
def DL_Mpc(zz):
    I = quad(lambda x: 1/np.sqrt(Om*(1+x)**3 + 1-Om), 0, zz)[0]
    return (1+zz)*c_kms/H0*I
DL = np.array([DL_Mpc(v) for v in z]) * 3.0857e24  # cm
logL_uv = xuv + np.log10(4*np.pi) + 2*np.log10(DL)
logL_x  = xx  + np.log10(4*np.pi) + 2*np.log10(DL)
print(f"logL_UV range: {logL_uv.min():.2f}..{logL_uv.max():.2f}")

# ---------- D'Agostini likelihood ----------
def nll(theta, x, y, ex, ey, xpiv):
    g, b, dl = theta
    if not (0 < g < 1.5 and -50 < b < 50 and 0 < dl < 1): return np.inf
    s2 = dl**2 + ey**2 + (g*ex)**2
    r = y - g*(x - xpiv) - b
    return 0.5*np.sum(r**2/s2 + np.log(2*np.pi*s2))

def fit(x, y, ex, ey, xpiv=0.0, nwalk=32, nstep=3000, burn=800, mcmc=True):
    best = minimize(nll, [0.6, np.mean(y - 0.6*(x-xpiv)), 0.25], args=(x, y, ex, ey, xpiv), method="Nelder-Mead",
                    options=dict(xatol=1e-6, fatol=1e-6, maxiter=20000))
    if not mcmc:
        return best.x, None
    p0 = best.x + 1e-3*rng.standard_normal((nwalk, 3))
    sam = emcee.EnsembleSampler(nwalk, 3, lambda t: -nll(t, x, y, ex, ey, xpiv))
    sam.run_mcmc(p0, nstep, progress=False)
    ch = sam.get_chain(discard=burn, flat=True)
    med = np.median(ch, axis=0); lo = np.percentile(ch, 16, axis=0); hi = np.percentile(ch, 84, axis=0)
    try: tau = sam.get_autocorr_time(quiet=True)
    except Exception: tau = None
    return med, dict(lo=lo, hi=hi, sig=(hi-lo)/2, tau=tau)

out = {}

# ---------- 1. full-sample luminosity-space fit (what Fig.1/Fig.2 claim) ----------
med, st = fit(logL_uv, logL_x, exuv, exx, xpiv=0.0)
print("\n[1] FULL SAMPLE, luminosity space, log L_X = g*log L_UV + b  (H0=70, Om=0.315)")
print(f"    gamma = {med[0]:.4f} +- {st['sig'][0]:.4f}   beta = {med[1]:.3f} +- {st['sig'][1]:.3f}   delta = {med[2]:.4f} +- {st['sig'][2]:.4f}")
out["full_lum"] = dict(gamma=[med[0], st['sig'][0]], beta=[med[1], st['sig'][1]], delta=[med[2], st['sig'][2]])
# pivoted (removes gamma-beta degeneracy)
piv = 30.0
medp, stp = fit(logL_uv, logL_x, exuv, exx, xpiv=piv)
print(f"    pivoted at logL_UV=30: gamma = {medp[0]:.4f} +- {stp['sig'][0]:.4f}  delta = {medp[2]:.4f}")

# ---------- 2. four equal-occupancy z bins, FLUX-FLUX space ----------
edges = np.quantile(z, [0, .25, .5, .75, 1.0])
print("\n[2] FOUR EQUAL-COUNT BINS, flux-flux space (cosmology-free)")
print("    edges:", np.round(edges, 3))
gam, gsig, zmid, binrows = [], [], [], []
for i in range(4):
    m_ = (z >= edges[i]) & (z <= edges[i+1]) if i == 3 else (z >= edges[i]) & (z < edges[i+1])
    xb, yb = xuv[m_], xx[m_]
    xp = np.median(xb)
    med, st = fit(xb, yb, exuv[m_], exx[m_], xpiv=xp, nstep=2500, burn=600)
    gam.append(med[0]); gsig.append(st['sig'][0]); zmid.append(np.mean(z[m_]))
    binrows.append(dict(zlo=edges[i], zhi=edges[i+1], n=int(m_.sum()), zmean=float(np.mean(z[m_])),
                        gamma=float(med[0]), gsig=float(st['sig'][0]), delta=float(med[2]), dsig=float(st['sig'][2])))
    print(f"    {edges[i]:.2f}-{edges[i+1]:.2f}  N={m_.sum():4d}  <z>={np.mean(z[m_]):.2f}  gamma={med[0]:.3f}+-{st['sig'][0]:.3f}  delta={med[2]:.3f}+-{st['sig'][2]:.3f}")
gam, gsig, zmid = map(np.array, (gam, gsig, zmid))
W = 1/gsig**2
A = np.vstack([np.ones(4), zmid]).T
cov = np.linalg.inv(A.T @ (W[:, None]*A))
g0, m = cov @ A.T @ (W*gam)
chi2 = np.sum(W*(gam - g0 - m*zmid)**2)
print(f"    gamma(z)=g0+m z :  g0={g0:.4f}+-{np.sqrt(cov[0,0]):.4f}   m={m:.4f}+-{np.sqrt(cov[1,1]):.4f}   chi2/dof={chi2:.2f}/2")
# constant-gamma test
gw = np.sum(W*gam)/np.sum(W); chi2c = np.sum(W*(gam-gw)**2)
print(f"    weighted-mean gamma = {gw:.4f}+-{1/np.sqrt(np.sum(W)):.4f}  chi2(const)/dof={chi2c:.2f}/3")
out["bins4"] = binrows; out["m_lin"] = [float(m), float(np.sqrt(cov[1,1]))]; out["g0"] = [float(g0), float(np.sqrt(cov[0,0]))]

# ---------- 2b. paper's other bin scheme (Sec 3.2: 0.1-0.7-1.2-2.0-5.0) ----------
print("\n[2b] PAPER'S SEC-3.2 BINS (0.1,0.7,1.2,2.0,5.0), flux-flux")
for lo_, hi_ in [(0.1,0.7),(0.7,1.2),(1.2,2.0),(2.0,5.0)]:
    m_ = (z>=lo_)&(z<hi_)
    xp = np.median(xuv[m_])
    med, _ = fit(xuv[m_], xx[m_], exuv[m_], exx[m_], xpiv=xp, mcmc=False)
    print(f"    {lo_}-{hi_}  N={m_.sum():4d}  gamma(MLE)={med[0]:.3f}")

# ---------- 3. finer bins (Risaliti&Lusso-style), 10 equal-count ----------
print("\n[3] TEN EQUAL-COUNT BINS, flux-flux (MLE + bootstrap err)")
e10 = np.quantile(z, np.linspace(0,1,11))
g10, s10, z10 = [], [], []
for i in range(10):
    m_ = (z>=e10[i])&(z<=e10[i+1]) if i==9 else (z>=e10[i])&(z<e10[i+1])
    xb, yb, exb, eyb = xuv[m_], xx[m_], exuv[m_], exx[m_]; xp=np.median(xb)
    g = fit(xb, yb, exb, eyb, xpiv=xp, mcmc=False)[0][0]
    bs = []
    for _ in range(60):
        j = rng.integers(0, len(xb), len(xb))
        bs.append(fit(xb[j], yb[j], exb[j], eyb[j], xpiv=xp, mcmc=False)[0][0])
    g10.append(g); s10.append(np.std(bs)); z10.append(np.mean(z[m_]))
    print(f"    <z>={np.mean(z[m_]):.2f} N={m_.sum()} gamma={g:.3f}+-{np.std(bs):.3f}")
g10, s10, z10 = map(np.array,(g10,s10,z10))
W=1/s10**2; A=np.vstack([np.ones(10), z10]).T; cov=np.linalg.inv(A.T@(W[:,None]*A)); g0b, mb = cov@A.T@(W*g10)
print(f"    10-bin: g0={g0b:.4f}+-{np.sqrt(cov[0,0]):.4f}  m={mb:.4f}+-{np.sqrt(cov[1,1]):.4f}")
out["bins10"] = dict(z=z10.tolist(), gamma=g10.tolist(), sig=s10.tolist(), m=[float(mb), float(np.sqrt(cov[1,1]))])

# ---------- 4. mock test: does an X-ray flux limit fake an evolution? ----------
print("\n[4] MOCKS: truth gamma=0.60 no evolution, per-bin X-ray flux floor = observed 2nd percentile")
gtrue = 0.60
nm = 300
mres = []
for k in range(nm):
    gm_bins = []
    for i in range(4):
        m_ = (z>=edges[i])&(z<=edges[i+1]) if i==3 else (z>=edges[i])&(z<edges[i+1])
        xb = xuv[m_]; xp=np.median(xb); ex=exuv[m_]; ey=exx[m_]
        # generate 3x sample, truncate at observed floor to mimic selection, keep first N_bin
        xm = np.tile(xb,3); exm=np.tile(ex,3); eym=np.tile(ey,3)
        btrue = np.median(xx[m_]) - gtrue*0  # intercept at pivot
        ym = btrue + gtrue*(xm-xp) + rng.normal(0,0.25,len(xm)) + rng.normal(0,1,len(xm))*eym
        xo = xm + rng.normal(0,1,len(xm))*exm
        floor = np.percentile(xx[m_], 2)
        keep = ym > floor
        xo, ym, exm, eym = xo[keep][:m_.sum()], ym[keep][:m_.sum()], exm[keep][:m_.sum()], eym[keep][:m_.sum()]
        gm_bins.append(fit(xo, ym, exm, eym, xpiv=xp, mcmc=False)[0][0])
    gm_bins=np.array(gm_bins)
    A=np.vstack([np.ones(4), zmid]).T
    mres.append(np.linalg.lstsq(A, gm_bins, rcond=None)[0][1])
mres=np.array(mres)
print(f"    recovered m over {nm} mocks: mean={mres.mean():.4f}  std={mres.std():.4f}  (bias in m = {mres.mean():.4f})")
out["mock_m"] = [float(mres.mean()), float(mres.std())]

json.dump(out, open("repro_results.json","w"), indent=1)
print("\nsaved repro_results.json")
