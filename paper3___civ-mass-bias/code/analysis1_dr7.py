#!/usr/bin/env python
"""
Paper 3, analysis 1: do Hβ, Mg II and C IV line widths / virial masses agree?
Data: Shen et al. (2011) SDSS DR7 quasar catalog, dr7_bh_Nov19_2013.fits.gz
      (105,783 quasars; FWHM + single-epoch masses for Hβ, Mg II, C IV).

Every number in results/dr7_results.json is produced here. Nothing is typed by hand.
Run from paper3___civ-mass-bias/code:  ../../qenv/bin/python analysis1_dr7.py
"""
import json, sys, os
import numpy as np
from astropy.io import fits
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "dr7_bh_Nov19_2013.fits.gz")
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

d = fits.open(DATA)[1].data
N_all = len(d)
z = d["REDSHIFT"].astype(float)
bal = d["BAL_FLAG"].astype(int)

def good(fwhm_col, err_col, snr_min=3.0):
    f = d[fwhm_col].astype(float); e = d[err_col].astype(float)
    ok = np.isfinite(f) & np.isfinite(e) & (f > 0) & (e > 0) & (f / e >= snr_min) & (f >= 1000.0)
    return ok, f

ok_hb, f_hb = good("FWHM_BROAD_HB", "FWHM_BROAD_HB_ERR")
ok_mg, f_mg = good("FWHM_BROAD_MGII", "FWHM_BROAD_MGII_ERR")
ok_civ, f_civ = good("FWHM_CIV", "FWHM_CIV_ERR")
ok_civ &= (bal == 0)              # BAL quasars have absorbed C IV profiles

m_hb = d["LOGBH_HB_VP06"].astype(float)     # Vestergaard & Peterson 2006
m_mg = d["LOGBH_MGII_VO09"].astype(float)   # Vestergaard & Osmer 2009
m_mg_s10 = d["LOGBH_MGII_S10"].astype(float)  # Shen et al. 2011 recalibration
m_civ = d["LOGBH_CIV_VP06"].astype(float)   # Vestergaard & Peterson 2006
voff_civ = d["VOFF_CIV_PEAK"].astype(float)  # km/s; sign convention checked below

def pct(a, q): return float(np.percentile(a, q))

def zwin(ok):
    zz = z[ok]
    return {"N": int(ok.sum()), "z_min": round(pct(zz, 0), 3), "z_p1": round(pct(zz, 1), 3),
            "z_median": round(pct(zz, 50), 3), "z_p99": round(pct(zz, 99), 3), "z_max": round(pct(zz, 100), 3)}

def corr(x, y):
    """x, y linear FWHM; stats on log10 and on linear values."""
    lx, ly = np.log10(x), np.log10(y)
    sp = stats.spearmanr(lx, ly)
    pe_log = stats.pearsonr(lx, ly)
    pe_lin = stats.pearsonr(x, y)
    # OLS y|x in log space, and bisector of the two OLS lines
    b1, a1 = np.polyfit(lx, ly, 1)
    b2i, a2i = np.polyfit(ly, lx, 1); b2 = 1 / b2i; a2 = -a2i / b2i
    bis = ((b1 * b2 - 1 + np.sqrt((1 + b1**2) * (1 + b2**2))) / (b1 + b2))
    resid = ly - (a1 + b1 * lx)
    return {"N": int(len(x)),
            "spearman_rho": round(float(sp.correlation), 3), "spearman_p": float(sp.pvalue),
            "pearson_r_log": round(float(pe_log[0]), 3), "pearson_r_linear": round(float(pe_lin[0]), 3),
            "ols_slope_log": round(float(b1), 3), "ols_intercept_log": round(float(a1), 3),
            "bisector_slope_log": round(float(bis), 3),
            "ols_resid_scatter_dex": round(float(np.std(resid)), 3),
            "median_ratio_y_over_x": round(float(np.median(y / x)), 3)}

def mass_cmp(mx, my, label_x, label_y):
    dm = my - mx
    ok = np.isfinite(dm)
    dm = dm[ok]
    sp = stats.spearmanr(mx[ok], my[ok])
    b, a = np.polyfit(mx[ok], my[ok], 1)
    return {"N": int(ok.sum()), "x": label_x, "y": label_y,
            "median_offset_y_minus_x_dex": round(float(np.median(dm)), 3),
            "mean_offset_dex": round(float(np.mean(dm)), 3),
            "std_offset_dex": round(float(np.std(dm)), 3),
            "nmad_offset_dex": round(float(1.4826 * np.median(np.abs(dm - np.median(dm)))), 3),
            "p16_offset": round(pct(dm, 16), 3), "p84_offset": round(pct(dm, 84), 3),
            "frac_y_exceeds_x_by_0p5dex": round(float(np.mean(dm > 0.5)), 4),
            "frac_y_below_x_by_0p5dex": round(float(np.mean(dm < -0.5)), 4),
            "spearman_rho_masses": round(float(sp.correlation), 3),
            "ols_slope": round(float(b), 3), "ols_intercept": round(float(a), 3)}

out = {"catalog": "Shen et al. 2011, ApJS 194, 45 (dr7_bh_Nov19_2013.fits.gz)", "N_catalog": N_all,
       "quality_cut": "FWHM>=1000 km/s, FWHM/err>=3, finite; C IV additionally BAL_FLAG==0",
       "coverage": {"Hbeta": zwin(ok_hb), "MgII": zwin(ok_mg), "CIV": zwin(ok_civ)}}

# ---- which pairs can exist in one SDSS spectrum? ----
both_hb_mg = ok_hb & ok_mg
both_mg_civ = ok_mg & ok_civ
both_hb_civ = ok_hb & ok_civ
all3 = ok_hb & ok_mg & ok_civ
out["pair_counts"] = {"Hbeta_and_MgII": int(both_hb_mg.sum()), "MgII_and_CIV": int(both_mg_civ.sum()),
                      "Hbeta_and_CIV": int(both_hb_civ.sum()), "all_three": int(all3.sum())}
out["pair_zwindows"] = {"Hbeta_and_MgII": zwin(both_hb_mg) if both_hb_mg.any() else None,
                        "MgII_and_CIV": zwin(both_mg_civ) if both_mg_civ.any() else None,
                        "Hbeta_and_CIV": zwin(both_hb_civ) if both_hb_civ.any() else None}

# ---- FWHM correlations ----
out["fwhm_Hbeta_vs_MgII"] = corr(f_mg[both_hb_mg], f_hb[both_hb_mg])   # x = Mg II, y = Hβ
out["fwhm_CIV_vs_MgII"] = corr(f_mg[both_mg_civ], f_civ[both_mg_civ])  # x = Mg II, y = C IV
if both_hb_civ.sum() > 10:
    out["fwhm_CIV_vs_Hbeta"] = corr(f_hb[both_hb_civ], f_civ[both_hb_civ])

# ---- C IV blueshift sign convention (data-driven) ----
v = voff_civ[both_mg_civ]
vok = np.isfinite(v) & (np.abs(v) < 6000)
out["civ_voff_sign_check"] = {"median_VOFF_CIV_PEAK_kms": round(float(np.median(v[vok])), 1),
    "note": "C IV is blueshifted on average in luminous quasars; the sign of the median tells which "
            "sign the catalog uses for blueshift."}

# ---- virial mass comparisons ----
sel = both_hb_mg & np.isfinite(m_hb) & np.isfinite(m_mg) & (m_hb > 0) & (m_mg > 0)
out["mass_Hbeta_vs_MgII_VO09"] = mass_cmp(m_mg[sel], m_hb[sel], "logM(MgII,VO09)", "logM(Hbeta,VP06)")
sel2 = both_hb_mg & np.isfinite(m_hb) & np.isfinite(m_mg_s10) & (m_hb > 0) & (m_mg_s10 > 0)
out["mass_Hbeta_vs_MgII_S10"] = mass_cmp(m_mg_s10[sel2], m_hb[sel2], "logM(MgII,S10)", "logM(Hbeta,VP06)")
selc = both_mg_civ & np.isfinite(m_civ) & np.isfinite(m_mg) & (m_civ > 0) & (m_mg > 0)
out["mass_CIV_vs_MgII_VO09"] = mass_cmp(m_mg[selc], m_civ[selc], "logM(MgII,VO09)", "logM(CIV,VP06)")
selc2 = both_mg_civ & np.isfinite(m_civ) & np.isfinite(m_mg_s10) & (m_civ > 0) & (m_mg_s10 > 0)
out["mass_CIV_vs_MgII_S10"] = mass_cmp(m_mg_s10[selc2], m_civ[selc2], "logM(MgII,S10)", "logM(CIV,VP06)")

# ---- C IV − Mg II mass offset as a function of C IV blueshift ----
dm = (m_civ - m_mg)[selc]; vv = voff_civ[selc]
blue = vv if out["civ_voff_sign_check"]["median_VOFF_CIV_PEAK_kms"] > 0 else -vv  # make blueshift positive
edges = [-5000, -500, 0, 500, 1000, 1500, 2000, 3000, 6000]
rows = []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (blue >= lo) & (blue < hi) & np.isfinite(blue)
    if m.sum() >= 30:
        rows.append({"blueshift_lo": lo, "blueshift_hi": hi, "N": int(m.sum()),
                     "median_dlogM_CIV_minus_MgII": round(float(np.median(dm[m])), 3),
                     "nmad": round(float(1.4826 * np.median(np.abs(dm[m] - np.median(dm[m])))), 3),
                     "median_FWHM_CIV": round(float(np.median(f_civ[selc][m])), 0),
                     "median_FWHM_MgII": round(float(np.median(f_mg[selc][m])), 0)})
out["civ_offset_vs_blueshift"] = rows
okb = np.isfinite(blue)
out["civ_offset_blueshift_spearman"] = {"rho": round(float(stats.spearmanr(blue[okb], dm[okb]).correlation), 3),
                                         "N": int(okb.sum())}

# ---- same blueshift table against the catalog's fiducial Mg II calibration (S10) ----
dm_s10 = (m_civ - m_mg_s10)[selc2]; vv2 = voff_civ[selc2]
blue2 = vv2 if out["civ_voff_sign_check"]["median_VOFF_CIV_PEAK_kms"] > 0 else -vv2
rows_s10 = []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (blue2 >= lo) & (blue2 < hi) & np.isfinite(blue2)
    if m.sum() >= 30:
        rows_s10.append({"blueshift_lo": lo, "blueshift_hi": hi, "N": int(m.sum()),
                         "median_dlogM_CIV_minus_MgII_S10": round(float(np.median(dm_s10[m])), 3),
                         "nmad": round(float(1.4826 * np.median(np.abs(dm_s10[m] - np.median(dm_s10[m])))), 3)})
out["civ_offset_vs_blueshift_S10"] = rows_s10

# ---- C IV blueshift relative to the broad Mg II line (VOFF_CIV_PEAK - VOFF_BROAD_MGII), fiducial Mg II mass ----
voff_mg = d["VOFF_BROAD_MGII"].astype(float)
rel = voff_civ - voff_mg
okr = selc2 & np.isfinite(rel) & (np.abs(rel) < 8000)
dmr = (m_civ - m_mg_s10)[okr]; rr = rel[okr]
out["civ_rel_mgii_blueshift_stats"] = {"N": int(okr.sum()), "median_kms": round(float(np.median(rr)), 1),
    "frac_gt_1000": round(float(np.mean(rr > 1000)), 4), "frac_gt_1500": round(float(np.mean(rr > 1500)), 4),
    "frac_gt_2000": round(float(np.mean(rr > 2000)), 4),
    "spearman_rho_offset_vs_blueshift": round(float(stats.spearmanr(rr, dmr).correlation), 3)}
rows_rel = []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (rr >= lo) & (rr < hi)
    if m.sum() >= 30:
        rows_rel.append({"blueshift_lo": lo, "blueshift_hi": hi, "N": int(m.sum()),
                         "median_dlogM_CIV_minus_MgII_S10": round(float(np.median(dmr[m])), 3),
                         "nmad": round(float(1.4826 * np.median(np.abs(dmr[m] - np.median(dmr[m])))), 3)})
out["civ_offset_vs_rel_blueshift_S10"] = rows_rel

# ---- conditional spread: at fixed Mg II mass, what range of C IV masses is assigned? ----
cond = []
for lo, hi in [(8.0, 8.5), (8.5, 9.0), (9.0, 9.5), (9.5, 10.0)]:
    m = (m_mg[selc] >= lo) & (m_mg[selc] < hi)
    if m.sum() >= 50:
        c = m_civ[selc][m]
        cond.append({"logM_MgII_lo": lo, "logM_MgII_hi": hi, "N": int(m.sum()), "CIV_p16": round(pct(c, 16), 2),
                     "CIV_p50": round(pct(c, 50), 2), "CIV_p84": round(pct(c, 84), 2),
                     "CIV_p2.5": round(pct(c, 2.5), 2), "CIV_p97.5": round(pct(c, 97.5), 2)})
out["civ_mass_given_mgii_mass"] = cond
# fraction of the Mg II–C IV sample with large blueshift (where the bias exceeds ~0.4 dex)
out["civ_blueshift_fractions"] = {"frac_blueshift_gt_1000": round(float(np.mean(blue[okb] > 1000)), 4),
                                  "frac_blueshift_gt_1500": round(float(np.mean(blue[okb] > 1500)), 4),
                                  "frac_blueshift_gt_2000": round(float(np.mean(blue[okb] > 2000)), 4),
                                  "median_blueshift_kms": round(float(np.median(blue[okb])), 1)}
# luminosity dependence of the C IV offset (z>6 quasars are the most luminous objects)
l1350 = d["LOGL1350"].astype(float)[selc]
lrows = []
for lo, hi in [(44.5, 45.5), (45.5, 46.0), (46.0, 46.5), (46.5, 47.5)]:
    m = (l1350 >= lo) & (l1350 < hi)
    if m.sum() >= 50:
        lrows.append({"logL1350_lo": lo, "logL1350_hi": hi, "N": int(m.sum()),
                      "median_dlogM_CIV_minus_MgII": round(float(np.median(dm[m])), 3),
                      "median_blueshift_kms": round(float(np.nanmedian(blue[m])), 0)})
out["civ_offset_vs_L1350"] = lrows

# ---- FWHM distributions per redshift bin (replaces the draft's unsourced claims) ----
def zbins(ok, f, edges):
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = ok & (z >= lo) & (z < hi)
        if m.sum() >= 20:
            rows.append({"z_lo": lo, "z_hi": hi, "N": int(m.sum()), "median_FWHM": round(float(np.median(f[m])), 0),
                         "p16": round(pct(f[m], 16), 0), "p84": round(pct(f[m], 84), 0), "p99": round(pct(f[m], 99), 0)})
    return rows
out["fwhm_vs_z"] = {"Hbeta": zbins(ok_hb, f_hb, [0, 0.25, 0.5, 0.75, 1.0]),
                    "MgII": zbins(ok_mg, f_mg, [0.3, 0.7, 1.0, 1.5, 2.0, 2.5]),
                    "CIV": zbins(ok_civ, f_civ, [1.5, 2.0, 2.5, 3.0, 3.5, 5.0])}

# ---- bimodality test for the draft's "double-peaked Mg II at 1<z<2" ----
m = ok_mg & (z >= 1) & (z < 2)
hist, be = np.histogram(f_mg[m], bins=np.arange(1000, 15001, 500))
out["mgii_1lt_z_lt2_hist"] = {"bin_edges_kms": be.tolist(), "counts": hist.tolist(),
                               "n_local_maxima_above_half_peak": int(np.sum((hist[1:-1] > hist[:-2]) & (hist[1:-1] > hist[2:]) & (hist[1:-1] > 0.5 * hist.max())))}

with open(os.path.join(RES, "dr7_results.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps({k: out[k] for k in ("coverage", "pair_counts", "fwhm_Hbeta_vs_MgII", "fwhm_CIV_vs_MgII",
                                      "mass_Hbeta_vs_MgII_VO09", "mass_CIV_vs_MgII_VO09", "civ_voff_sign_check",
                                      "civ_offset_blueshift_spearman")}, indent=1))

# ---- figures ----
def hexpanel(ax, x, y, xl, yl, title, lim):
    ax.hexbin(np.log10(x), np.log10(y), gridsize=60, bins="log", cmap="Greys", mincnt=1)
    ax.plot(lim, lim, "r--", lw=1, label="1:1")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(title, fontsize=10)
    ax.legend(loc="upper left", fontsize=8)

fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
r = out["fwhm_Hbeta_vs_MgII"]
hexpanel(axs[0], f_mg[both_hb_mg], f_hb[both_hb_mg], "log FWHM(Mg II) [km/s]", "log FWHM(Hβ) [km/s]",
         f"DR7, {r['N']} quasars, 0.4<z<0.9: ρ={r['spearman_rho']}", [3.0, 4.2])
r = out["fwhm_CIV_vs_MgII"]
hexpanel(axs[1], f_mg[both_mg_civ], f_civ[both_mg_civ], "log FWHM(Mg II) [km/s]", "log FWHM(C IV) [km/s]",
         f"DR7, {r['N']} non-BAL quasars, 1.5<z<2.2: ρ={r['spearman_rho']}", [3.0, 4.2])
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig1_fwhm_pairs_dr7.png"), dpi=160); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
r = out["mass_Hbeta_vs_MgII_VO09"]
hexpanel(axs[0], 10**m_mg[sel], 10**m_hb[sel], "log M(Mg II, VO09)", "log M(Hβ, VP06)",
         f"median offset {r['median_offset_y_minus_x_dex']:+.2f} dex, NMAD {r['nmad_offset_dex']:.2f}", [7, 10.5])
r = out["mass_CIV_vs_MgII_VO09"]
hexpanel(axs[1], 10**m_mg[selc], 10**m_civ[selc], "log M(Mg II, VO09)", "log M(C IV, VP06)",
         f"median offset {r['median_offset_y_minus_x_dex']:+.2f} dex, NMAD {r['nmad_offset_dex']:.2f}", [7, 10.5])
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2_mass_pairs_dr7.png"), dpi=160); plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 4.2))
xs = [(rw["blueshift_lo"] + rw["blueshift_hi"]) / 2 for rw in rows]
ys = [rw["median_dlogM_CIV_minus_MgII"] for rw in rows]
es = [rw["nmad"] / np.sqrt(rw["N"]) for rw in rows]
ax.errorbar(xs, ys, yerr=es, fmt="o-", color="k")
ax.axhline(0, color="r", ls="--", lw=1)
ax.set_xlabel("C IV blueshift [km/s] (positive = blue)"); ax.set_ylabel("median log M(C IV) − log M(Mg II) [dex]")
ax.set_title("DR7: C IV mass excess grows with C IV blueshift", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_civ_offset_vs_blueshift_dr7.png"), dpi=160); plt.close(fig)
print("figures written")
