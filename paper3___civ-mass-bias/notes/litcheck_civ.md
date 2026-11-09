# Literature check: C IV virial masses, SDSS coverage, catalogs (agent report, 2026-11-09)

## SDSS coverage
- SDSS-I/II spectrograph 3800–9200 Å (SkyServer DR7; Wu & Shen 2022 text); BOSS/eBOSS 3600–10400 Å (Smee et al. 2013, AJ 146, 32; sdss4.org).
- Hβ 4861 and C IV 1549 need λmax/λmin ≥ 3.14; SDSS gives 2.42, BOSS 2.89 → no single SDSS/BOSS spectrum can contain both lines.
- Shen et al. 2011 fit windows: Hα z≤0.39, Hβ z≤0.89, Mg II 0.35≤z≤2.25, C IV 1.5≤z≤4.95 (ar5iv 1006.5178). Velocity shifts: "Positive values indicate blueshift and negative values indicate redshift."

## Catalogs
- Wu & Shen 2022, ApJS 263, 42 (arXiv:2209.03987): 750,414 broad-line quasars from DR16Q; line arrays (HBETA_BR, MGII, CIV …) with 6 elements [peak λ, 50% centroid λ, flux, logL, FWHM, rest EW]; LOGMBH_HB / _MGII / _CIV / LOGMBH; fiducial recipes = Shen et al. 2011 (Hβ VP06; Mg II S10; C IV VP06); fiducial Hβ z<0.7, Mg II 0.7≤z<1.9 (legacy) / <2.0 (BOSS), C IV above. L_Edd = 1.3e38 (M/M☉). File dr16q_prop_May01_2024.fits.gz, 1,391,052,577 bytes, http://quasar.astro.illinois.edu/paper_data/DR16Q/.
- Shen et al. 2011, ApJS 194, 45 (arXiv:1006.5178): 105,783 quasars. Recipes log M = a + b log(λLλ/1e44) + 2 log FWHM: Hβ VP06 a=0.910 b=0.50; Mg II S10 a=0.740 b=0.62; Mg II VO09 a=0.860 b=0.50; C IV VP06 a=0.660 b=0.53. Mg II(S10)–Hβ mean offset <0.01 dex by construction. File dr7_bh_Nov19_2013.fits.gz (62,749,630 bytes).

## C IV bias literature (verified numbers)
- Baskin & Laor 2005, MNRAS 356, 1029: C IV broader than Hβ when FWHM(Hβ)<4000, narrower when >4000; C IV–Hβ mass dispersion 0.505 dex; "significantly less accurate, and possibly biased".
- Shen et al. 2008, ApJ 680, 169: ~60,000 quasars; log FWHM(Hβ)/FWHM(Mg II) mean 0.006, dispersion 0.11 dex; mass ratio mean 0.034, dispersion 0.22 dex; C IV–Mg II FWHM correlation "very weak"; log(M_MgII/M_CIV) mean −0.06, dispersion 0.34 dex, systematic with blueshift.
- Shen & Liu 2012, ApJ 753, 125: 60 quasars z~1.5–2.2; r_s FWHM(C IV)–FWHM(Hβ) = 0.11 (P=0.39); r_s FWHM(Mg II)–FWHM(Hβ) = 0.64 (P=2.8e-8); C IV mass scatter vs Hβ 0.26–0.28 dex, Mg II 0.15–0.16 dex.
- Denney 2012, ApJ 759, 44: non-variable C IV component; shape correction halves scatter to ~0.2 dex.
- Runnoe et al. 2013, MNRAS 434, 848: Hβ–C IV mass scatter 0.43 → 0.33 dex with 1400 Å/C IV peak-ratio correction.
- Coatman et al. 2016, MNRAS 461, 647: FWHM masses overestimated by factor ~5 at blueshift >2000 km/s; C IV masses LOWER at blueshift <1200 km/s.
- Coatman et al. 2017, MNRAS 465, 2120: 230 quasars 1.5<z<4; C IV masses ×5 at 3000 km/s blueshift, ~×10 at >5000; scatter 0.40 → ~0.2 dex after correction; FWHM_corr = FWHM/[(0.41±0.02)(Δv/1e3)+(0.62±0.04)] (Hα-based; Hβ-based 0.36, 0.61); M = 10^6.71 (FWHM_corr/1e3)^2 (λL1350/1e44)^0.53.
- Park et al. 2017, ApJ 839, 93: C IV FWHM calibration slope on FWHM 0.50 (not 2), scatter 0.37 dex; recommends σ_line.
- Mejía-Restrepo et al. 2018, Nat. Astron. 2, 63: virial factor f ∝ FWHM^−1.2; Mg II correction lowers brightest z≳6 masses by up to ×2.
- Marziani et al. 2019, A&A 627, A88: C IV FWHM "not immediately offering a virial broadening estimator"; corrected C IV gives unbiased masses with ~0.3 dex scatter.
- Richards et al. 2011, AJ 141, 167: C IV blueshift mean 810 km/s (radio-quiet), 360 (radio-loud); RM sample covers only part of C IV parameter space.
- Wang et al. 2009, ApJ 707, 1334: log FWHM(Mg II) = 0.81 log FWHM(Hβ) + const, scatter 0.08 dex.

## Calibrations
- VP06 (ApJ 641, 689): C IV log M = log[(FWHM/1e3)^2 (λL1350/1e44)^0.53] + 6.66±0.01, scatter 0.36 dex; Hβ zero point 6.91, 0.43 dex; "factor of 4" absolute uncertainty. Zero point = weighted mean vs RM masses → mean offset zero on the calibration sample by construction.
- VO09 (ApJ 699, 800): Mg II zp 6.86 (3000 Å), scatter 0.55 dex; calibrated to SDSS Hβ and C IV masses, consistent within ~0.1 dex.

## Γ – λ_Edd
- Shemmer et al. 2008, ApJ 682, 81: Γ = 0.31 log λ + 2.11; r_S=0.55; scatter in predicted log λ_Edd ~0.35 dex; masses to factor ~2–3.
- Risaliti, Young & Elvis 2009, ApJ 700, L6: 343 quasars; slope 0.31±0.06, dispersion 0.33–0.4 dex, r=0.32; NO correlation for C IV-mass subsample (−0.11±0.17).
- Brightman et al. 2013, MNRAS 433, 2485: Γ = 0.32 log λ + 2.27, r_S=0.60; "not viable for single sources".
- Trakhtenbrot et al. 2017, MNRAS 470, 800: 228 BASS AGN; r_s=0.31 (weak); no correlation for direct/RM masses; inverse relation 1σ range 1.37 dex; intrinsic scatter >0.6 dex.
- L_Edd: SDSS catalogs use 1.3e38 erg/s per M☉; 1.26e38 is the pure-hydrogen Thomson value.

## Fan, Bañados & Simcoe 2023, ARA&A 61, 373 (§4.1)
- Reviews "only ... results based on Mg II measurements, which are thought to be the most reliable until Hβ-derived masses are possible"; Mg II relation intrinsic scatter 0.55 dex; C IV "uncertain", blueshifts >2000 km/s more common at z>6.
