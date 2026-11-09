#!/usr/bin/env python
"""
Paper 3, analysis 2: same tests as analysis1_dr7.py on the SDSS DR16Q catalog of
Wu & Shen (2022, ApJS 263, 42), dr16q_prop_May01_2024.fits.gz (750,414 quasars).

Line columns are 6-element arrays [peak λ, 50% centroid λ, flux, logL, FWHM, rest EW];
FWHM is element 4. Masses: LOGMBH_HB, LOGMBH_MGII, LOGMBH_CIV (Shen et al. 2011 recipes:
Hβ VP06, Mg II S10, C IV VP06).

Every number in results/dr16q_results.json is produced here.
Run from paper3___civ-mass-bias/code:  ../../qenv/bin/python analysis2_dr16q.py
"""
import json, os
import numpy as np
from astropy.io import fits
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "dr16q_prop_May01_2024.fits")
RES = os.path.join(HERE, "..", "results"); FIG = os.path.join(HERE, "..", "figures")
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

hdu = fits.open(DATA, memmap=True)[1]
names = hdu.columns.names
d = hdu.data
N_all = len(d)
print("rows", N_all)

def pick(*cands):
    for c in cands:
        if c in names: return c
    raise KeyError(cands)

zcol = pick("Z_DR16Q", "Z", "Z_FIT", "Z_SYS")
z = np.asarray(d[zcol], float)
balcol = None  # the properties catalog carries no BAL flag (BAL_PROB lives in the parent DR16Q catalog)
bal = np.zeros(N_all)
print("z column", zcol, "bal column", balcol, "names with FWHM-bearing arrays:",
      [n for n in names if n in ("HBETA_BR", "MGII_BR", "MGII", "CIV_BR", "CIV", "HBETA")])

def line_fwhm(col):
    arr = np.asarray(d[col], float)
    err = np.asarray(d[col + "_ERR"], float)
    return arr[:, 4], err[:, 4]

hb_col = pick("HBETA_BR"); mg_col = pick("MGII_BR", "MGII"); civ_col = pick("CIV_BR", "CIV")
f_hb, e_hb = line_fwhm(hb_col); f_mg, e_mg = line_fwhm(mg_col); f_civ, e_civ = line_fwhm(civ_col)

def good(f, e, snr_min=3.0):
    return np.isfinite(f) & np.isfinite(e) & (f >= 1000) & (e > 0) & (f / e >= snr_min)

ok_hb = good(f_hb, e_hb); ok_mg = good(f_mg, e_mg); ok_civ = good(f_civ, e_civ)
# restrict each line to the redshift window where line + local continuum sit inside the BOSS bandpass
# (3600-10400 A); outside it the catalog still carries edge fits that are not usable.
ok_hb &= (z >= 0.0) & (z <= 1.0); ok_mg &= (z >= 0.35) & (z <= 2.5); ok_civ &= (z >= 1.5) & (z <= 5.0)
sn_all = np.asarray(d["SN_MEDIAN_ALL"], float)
# BAL_PROB in DR16Q is a probability (0–1); legacy BAL_FLAG is 0/1. Exclude likely BALs from C IV.
# no BAL exclusion possible here; DR7 (analysis1) is the BAL-cleaned check

m_hb = np.asarray(d["LOGMBH_HB"], float); m_mg = np.asarray(d["LOGMBH_MGII"], float); m_civ = np.asarray(d["LOGMBH_CIV"], float)
for m in (m_hb, m_mg, m_civ):
    m[~np.isfinite(m) | (m <= 0)] = np.nan

# C IV blueshift: peak wavelength (element 0) vs rest 1549.06 Å, in the quasar rest frame (positive = blue)
civ_arr = np.asarray(d[civ_col], float)
peak = civ_arr[:, 0]                       # rest-frame peak wavelength in the Z_FIT frame
z_fit = np.asarray(d["Z_FIT"], float); z_sys = np.asarray(d["Z_SYS"], float)
c_kms = 299792.458
lam_sys = peak * (1 + z_fit) / (1 + z_sys)  # peak wavelength in the systemic (Z_SYS) frame
blue = np.where(np.isfinite(lam_sys) & (peak > 0) & (z_sys > 0), (1549.06 - lam_sys) / 1549.06 * c_kms, np.nan)

def pct(a, q): return float(np.percentile(a, q))
def zwin(ok):
    zz = z[ok]
    return {"N": int(ok.sum()), "z_min": round(pct(zz, 0), 3), "z_p1": round(pct(zz, 1), 3),
            "z_median": round(pct(zz, 50), 3), "z_p99": round(pct(zz, 99), 3), "z_max": round(pct(zz, 100), 3)}

def corr(x, y):
    lx, ly = np.log10(x), np.log10(y)
    sp = stats.spearmanr(lx, ly); pl = stats.pearsonr(lx, ly); pn = stats.pearsonr(x, y)
    b1, a1 = np.polyfit(lx, ly, 1)
    resid = ly - (a1 + b1 * lx)
    return {"N": int(len(x)), "spearman_rho": round(float(sp.correlation), 3), "spearman_p": float(sp.pvalue),
            "pearson_r_log": round(float(pl[0]), 3), "pearson_r_linear": round(float(pn[0]), 3),
            "ols_slope_log": round(float(b1), 3), "ols_intercept_log": round(float(a1), 3),
            "ols_resid_scatter_dex": round(float(np.std(resid)), 3),
            "median_ratio_y_over_x": round(float(np.median(y / x)), 3)}

def mass_cmp(mx, my, lx, ly):
    ok = np.isfinite(mx) & np.isfinite(my); dm = (my - mx)[ok]
    sp = stats.spearmanr(mx[ok], my[ok]); b, a = np.polyfit(mx[ok], my[ok], 1)
    return {"N": int(ok.sum()), "x": lx, "y": ly,
            "median_offset_y_minus_x_dex": round(float(np.median(dm)), 3), "mean_offset_dex": round(float(np.mean(dm)), 3),
            "std_offset_dex": round(float(np.std(dm)), 3),
            "nmad_offset_dex": round(float(1.4826 * np.median(np.abs(dm - np.median(dm)))), 3),
            "p16_offset": round(pct(dm, 16), 3), "p84_offset": round(pct(dm, 84), 3),
            "frac_y_exceeds_x_by_0p5dex": round(float(np.mean(dm > 0.5)), 4),
            "frac_y_below_x_by_0p5dex": round(float(np.mean(dm < -0.5)), 4),
            "spearman_rho_masses": round(float(sp.correlation), 3), "ols_slope": round(float(b), 3), "ols_intercept": round(float(a), 3)}

out = {"catalog": "Wu & Shen 2022, ApJS 263, 42 (dr16q_prop_May01_2024.fits.gz)", "N_catalog": N_all,
       "columns_used": {"z": zcol, "bal": "none available", "blueshift": "CIV peak (Z_FIT frame) -> Z_SYS frame", "Hbeta": hb_col, "MgII": mg_col, "CIV": civ_col,
                        "masses": ["LOGMBH_HB", "LOGMBH_MGII", "LOGMBH_CIV"]},
       "quality_cut": "FWHM>=1000 km/s, FWHM/err>=3, finite; windows Hb z<=1.0, MgII 0.35-2.5, CIV 1.5-5.0; NO BAL exclusion (flag not in this catalog)",
       "coverage": {"Hbeta": zwin(ok_hb), "MgII": zwin(ok_mg), "CIV": zwin(ok_civ)}}

both_hb_mg = ok_hb & ok_mg; both_mg_civ = ok_mg & ok_civ; both_hb_civ = ok_hb & ok_civ
out["mass_availability"] = {"Hbeta_pairs_with_both_masses": int((both_hb_mg & np.isfinite(m_hb) & np.isfinite(m_mg)).sum()),
                            "CIV_pairs_with_both_masses": int((both_mg_civ & np.isfinite(m_civ) & np.isfinite(m_mg)).sum()),
                            "note": "LOGMBH_* is NaN/0 where Wu & Shen 2022 did not compute a mass (e.g. continuum luminosity not covered)"}
out["pair_counts"] = {"Hbeta_and_MgII": int(both_hb_mg.sum()), "MgII_and_CIV": int(both_mg_civ.sum()),
                      "Hbeta_and_CIV": int(both_hb_civ.sum()), "all_three": int((ok_hb & ok_mg & ok_civ).sum())}
out["pair_zwindows"] = {"Hbeta_and_MgII": zwin(both_hb_mg), "MgII_and_CIV": zwin(both_mg_civ),
                        "Hbeta_and_CIV": zwin(both_hb_civ) if both_hb_civ.any() else None}
out["fwhm_Hbeta_vs_MgII"] = corr(f_mg[both_hb_mg], f_hb[both_hb_mg])
out["fwhm_CIV_vs_MgII"] = corr(f_mg[both_mg_civ], f_civ[both_mg_civ])
if both_hb_civ.sum() > 10:
    out["fwhm_CIV_vs_Hbeta"] = corr(f_hb[both_hb_civ], f_civ[both_hb_civ])
hs = sn_all >= 10
out["fwhm_Hbeta_vs_MgII_SN10"] = corr(f_mg[both_hb_mg & hs], f_hb[both_hb_mg & hs])
out["fwhm_CIV_vs_MgII_SN10"] = corr(f_mg[both_mg_civ & hs], f_civ[both_mg_civ & hs])

sel = both_hb_mg; out["mass_Hbeta_vs_MgII"] = mass_cmp(m_mg[sel], m_hb[sel], "logM(MgII,S10)", "logM(Hbeta,VP06)")
selc = both_mg_civ; out["mass_CIV_vs_MgII"] = mass_cmp(m_mg[selc], m_civ[selc], "logM(MgII,S10)", "logM(CIV,VP06)")

# C IV − Mg II offset vs blueshift
okd = selc & np.isfinite(m_civ) & np.isfinite(m_mg) & np.isfinite(blue) & (np.abs(blue) < 8000)
dm = (m_civ - m_mg)[okd]; bb = blue[okd]
out["civ_blueshift_stats"] = {"N": int(okd.sum()), "median_blueshift_kms": round(float(np.median(bb)), 1),
                              "frac_gt_1000": round(float(np.mean(bb > 1000)), 4), "frac_gt_1500": round(float(np.mean(bb > 1500)), 4),
                              "frac_gt_2000": round(float(np.mean(bb > 2000)), 4),
                              "spearman_rho_offset_vs_blueshift": round(float(stats.spearmanr(bb, dm).correlation), 3)}
rows = []
for lo, hi in zip([-5000, -500, 0, 500, 1000, 1500, 2000, 3000], [-500, 0, 500, 1000, 1500, 2000, 3000, 8000]):
    m = (bb >= lo) & (bb < hi)
    if m.sum() >= 30:
        rows.append({"blueshift_lo": lo, "blueshift_hi": hi, "N": int(m.sum()),
                     "median_dlogM_CIV_minus_MgII": round(float(np.median(dm[m])), 3),
                     "nmad": round(float(1.4826 * np.median(np.abs(dm[m] - np.median(dm[m])))), 3)})
out["civ_offset_vs_blueshift"] = rows

# C IV blueshift relative to the Mg II peak (same quantity as Schindler et al. 2020), no systemic z needed
mg_arr = np.asarray(d[mg_col], float); mg_peak = mg_arr[:, 0]
rel = np.where(np.isfinite(mg_peak) & (mg_peak > 0) & np.isfinite(peak) & (peak > 0),
               c_kms * ((1549.06 - peak) / 1549.06 - (2798.75 - mg_peak) / 2798.75), np.nan)
okr = selc & np.isfinite(m_civ) & np.isfinite(m_mg) & np.isfinite(rel) & (np.abs(rel) < 8000)
dmr = (m_civ - m_mg)[okr]; rr = rel[okr]
out["civ_rel_mgii_blueshift_stats"] = {"N": int(okr.sum()), "median_kms": round(float(np.median(rr)), 1),
    "frac_gt_1000": round(float(np.mean(rr > 1000)), 4), "frac_gt_1500": round(float(np.mean(rr > 1500)), 4),
    "frac_gt_2000": round(float(np.mean(rr > 2000)), 4),
    "spearman_rho_offset_vs_blueshift": round(float(stats.spearmanr(rr, dmr).correlation), 3)}
rows_rel = []
for lo, hi in zip([-5000, -500, 0, 500, 1000, 1500, 2000, 3000], [-500, 0, 500, 1000, 1500, 2000, 3000, 8000]):
    m = (rr >= lo) & (rr < hi)
    if m.sum() >= 30:
        rows_rel.append({"blueshift_lo": lo, "blueshift_hi": hi, "N": int(m.sum()),
                         "median_dlogM_CIV_minus_MgII": round(float(np.median(dmr[m])), 3),
                         "nmad": round(float(1.4826 * np.median(np.abs(dmr[m] - np.median(dmr[m])))), 3)})
out["civ_offset_vs_rel_blueshift"] = rows_rel

# ---- why might DR16Q differ from DR7? (a) S/N: redo the relative-blueshift table at SN_MEDIAN_ALL>=10;
#      (b) luminosity: redo it for log L1350 >= 46 (DR7-like); (c) host light in L5100 for the Hβ masses.
def rel_table(mask, key):
    rows_ = []
    for lo, hi in zip([-5000, -500, 0, 500, 1000, 1500, 2000, 3000], [-500, 0, 500, 1000, 1500, 2000, 3000, 8000]):
        m = mask & (rr >= lo) & (rr < hi)
        if m.sum() >= 30:
            rows_.append({"blueshift_lo": lo, "blueshift_hi": hi, "N": int(m.sum()),
                          "median_dlogM_CIV_minus_MgII": round(float(np.median(dmr[m])), 3),
                          "nmad": round(float(1.4826 * np.median(np.abs(dmr[m] - np.median(dmr[m])))), 3)})
    out[key] = rows_
sn_r = sn_all[okr]; l_r = np.asarray(d["LOGL1350"], float)[okr]
rel_table(sn_r >= 10, "civ_offset_vs_rel_blueshift_SN10")
rel_table(l_r >= 46.0, "civ_offset_vs_rel_blueshift_L46")
rel_table((l_r >= 46.0) & (sn_r >= 10), "civ_offset_vs_rel_blueshift_L46_SN10")
out["civ_offset_subsamples"] = {
    "all": {"N": int(len(dmr)), "median": round(float(np.median(dmr)), 3), "nmad": round(float(1.4826*np.median(np.abs(dmr-np.median(dmr)))), 3),
            "median_rel_blueshift": round(float(np.median(rr)), 0), "spearman_offset_vs_blueshift": round(float(stats.spearmanr(rr, dmr).correlation), 3)},
    "SN10": {"N": int((sn_r >= 10).sum()), "median": round(float(np.median(dmr[sn_r >= 10])), 3),
             "spearman_offset_vs_blueshift": round(float(stats.spearmanr(rr[sn_r >= 10], dmr[sn_r >= 10]).correlation), 3)},
    "L46": {"N": int((l_r >= 46).sum()), "median": round(float(np.median(dmr[l_r >= 46])), 3),
            "median_rel_blueshift": round(float(np.median(rr[l_r >= 46])), 0),
            "spearman_offset_vs_blueshift": round(float(stats.spearmanr(rr[l_r >= 46], dmr[l_r >= 46]).correlation), 3)},
    "L46_SN10": {"N": int(((l_r >= 46) & (sn_r >= 10)).sum()), "median": round(float(np.median(dmr[(l_r >= 46) & (sn_r >= 10)])), 3),
                 "spearman_offset_vs_blueshift": round(float(stats.spearmanr(rr[(l_r >= 46) & (sn_r >= 10)], dmr[(l_r >= 46) & (sn_r >= 10)]).correlation), 3)}}

# (c) host-galaxy light: Wu & Shen 2022 give FHOST_5100 (host-to-total flux ratio at 5100 A, z<1 only).
# The Hβ mass scales as L5100^0.5, so the AGN-only mass is lower by 0.5*log10(1 - fhost).
fh = np.asarray(d["FHOST_5100"], float) if "FHOST_5100" in names else np.full(N_all, np.nan)
s_hb = sel & np.isfinite(m_hb) & np.isfinite(m_mg)
fh_ok = s_hb & np.isfinite(fh) & (fh >= 0) & (fh < 0.95)
corr_hb = m_hb + 0.5 * np.log10(1 - np.clip(fh, 0, 0.95))
dm_raw = (m_hb - m_mg)[fh_ok]; dm_corr = (corr_hb - m_mg)[fh_ok]
out["hbeta_host_correction"] = {"N_with_fhost": int(fh_ok.sum()), "median_fhost": round(float(np.median(fh[fh_ok])), 3),
    "median_dlogM_Hb_minus_MgII_raw": round(float(np.median(dm_raw)), 3),
    "median_dlogM_Hb_minus_MgII_host_corrected": round(float(np.median(dm_corr)), 3),
    "nmad_raw": round(float(1.4826*np.median(np.abs(dm_raw-np.median(dm_raw)))), 3),
    "nmad_host_corrected": round(float(1.4826*np.median(np.abs(dm_corr-np.median(dm_corr)))), 3),
    "note": "correction = 0.5*log10(1 - FHOST_5100) applied to log M(Hbeta); Mg II mass untouched"}

# luminosity dependence
l1350 = np.asarray(d["LOGL1350"], float)[okd]
lrows = []
for lo, hi in [(44.5, 45.5), (45.5, 46.0), (46.0, 46.5), (46.5, 47.0), (47.0, 48.0)]:
    m = (l1350 >= lo) & (l1350 < hi)
    if m.sum() >= 50:
        lrows.append({"logL1350_lo": lo, "logL1350_hi": hi, "N": int(m.sum()),
                      "median_dlogM_CIV_minus_MgII": round(float(np.median(dm[m])), 3),
                      "median_blueshift_kms": round(float(np.median(bb[m])), 0)})
out["civ_offset_vs_L1350"] = lrows

# conditional spread of C IV mass at fixed Mg II mass
cond = []
for lo, hi in [(8.0, 8.5), (8.5, 9.0), (9.0, 9.5), (9.5, 10.0)]:
    m = (m_mg[selc] >= lo) & (m_mg[selc] < hi) & np.isfinite(m_civ[selc])
    if m.sum() >= 50:
        c = m_civ[selc][m]
        cond.append({"logM_MgII_lo": lo, "logM_MgII_hi": hi, "N": int(m.sum()), "CIV_p16": round(pct(c, 16), 2),
                     "CIV_p50": round(pct(c, 50), 2), "CIV_p84": round(pct(c, 84), 2)})
out["civ_mass_given_mgii_mass"] = cond

# FWHM vs z per line
def zbins(ok, f, edges):
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = ok & (z >= lo) & (z < hi)
        if m.sum() >= 20:
            rows.append({"z_lo": lo, "z_hi": hi, "N": int(m.sum()), "median_FWHM": round(float(np.median(f[m])), 0),
                         "p16": round(pct(f[m], 16), 0), "p84": round(pct(f[m], 84), 0)})
    return rows
out["fwhm_vs_z"] = {"Hbeta": zbins(ok_hb, f_hb, [0, 0.25, 0.5, 0.75, 1.0, 1.2]),
                    "MgII": zbins(ok_mg, f_mg, [0.3, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0]),
                    "CIV": zbins(ok_civ, f_civ, [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0])}
m = ok_mg & (z >= 1) & (z < 2)
hist, be = np.histogram(f_mg[m], bins=np.arange(1000, 15001, 500))
out["mgii_1lt_z_lt2_hist"] = {"bin_edges_kms": be.tolist(), "counts": hist.tolist(),
    "n_local_maxima_above_half_peak": int(np.sum((hist[1:-1] > hist[:-2]) & (hist[1:-1] > hist[2:]) & (hist[1:-1] > 0.5 * hist.max())))}

with open(os.path.join(RES, "dr16q_results.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps({k: out[k] for k in ("coverage", "pair_counts", "mass_availability", "fwhm_Hbeta_vs_MgII", "fwhm_CIV_vs_MgII",
                                      "fwhm_Hbeta_vs_MgII_SN10", "fwhm_CIV_vs_MgII_SN10",
                                      "mass_Hbeta_vs_MgII", "mass_CIV_vs_MgII", "civ_blueshift_stats",
                                      "civ_rel_mgii_blueshift_stats", "civ_offset_vs_rel_blueshift", "civ_offset_vs_blueshift",
                                      "civ_offset_vs_L1350", "civ_mass_given_mgii_mass", "fwhm_vs_z", "civ_offset_vs_rel_blueshift_SN10", "civ_offset_vs_rel_blueshift_L46", "civ_offset_vs_rel_blueshift_L46_SN10", "civ_offset_subsamples", "hbeta_host_correction")}, indent=0))

# figures
def hexpanel(ax, x, y, xl, yl, title, lim):
    ax.hexbin(np.log10(x), np.log10(y), gridsize=70, bins="log", cmap="Greys", mincnt=1)
    ax.plot(lim, lim, "r--", lw=1, label="1:1"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(title, fontsize=10); ax.legend(loc="upper left", fontsize=8)
fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
r = out["fwhm_Hbeta_vs_MgII"]
hexpanel(axs[0], f_mg[both_hb_mg], f_hb[both_hb_mg], "log FWHM(Mg II) [km/s]", "log FWHM(Hβ) [km/s]",
         f"DR16Q, N={r['N']}: Spearman ρ={r['spearman_rho']}", [3.0, 4.2])
r = out["fwhm_CIV_vs_MgII"]
hexpanel(axs[1], f_mg[both_mg_civ], f_civ[both_mg_civ], "log FWHM(Mg II) [km/s]", "log FWHM(C IV) [km/s]",
         f"DR16Q (BALs included), N={r['N']}: Spearman ρ={r['spearman_rho']}", [3.0, 4.2])
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig1_fwhm_pairs_dr16q.png"), dpi=160); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
r = out["mass_Hbeta_vs_MgII"]; s = sel & np.isfinite(m_hb) & np.isfinite(m_mg)
hexpanel(axs[0], 10**m_mg[s], 10**m_hb[s], "log M(Mg II)", "log M(Hβ)",
         f"median offset {r['median_offset_y_minus_x_dex']:+.2f} dex, NMAD {r['nmad_offset_dex']:.2f}", [7, 10.5])
r = out["mass_CIV_vs_MgII"]; s = selc & np.isfinite(m_civ) & np.isfinite(m_mg)
hexpanel(axs[1], 10**m_mg[s], 10**m_civ[s], "log M(Mg II)", "log M(C IV)",
         f"median offset {r['median_offset_y_minus_x_dex']:+.2f} dex, NMAD {r['nmad_offset_dex']:.2f}", [7, 10.5])
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2_mass_pairs_dr16q.png"), dpi=160); plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 4.2))
for rws, lab, col in ((rows, "vs Z_SYS", "0.6"), (rows_rel, "relative to Mg II peak", "k")):
    xs = [(rw["blueshift_lo"] + rw["blueshift_hi"]) / 2 for rw in rws]; xs[-1] = 3500
    ys = [rw["median_dlogM_CIV_minus_MgII"] for rw in rws]; es = [rw["nmad"] / np.sqrt(rw["N"]) for rw in rws]
    ax.errorbar(xs, ys, yerr=es, fmt="o-", color=col, label=lab)
ax.axhline(0, color="r", ls="--", lw=1); ax.legend(fontsize=8)
ax.set_xlabel("C IV blueshift [km/s] (positive = blue)"); ax.set_ylabel("median log M(C IV) − log M(Mg II) [dex]")
ax.set_title("DR16Q: C IV mass excess vs C IV blueshift", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_civ_offset_vs_blueshift_dr16q.png"), dpi=160); plt.close(fig)
print("figures written")
