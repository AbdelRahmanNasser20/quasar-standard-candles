# Audit of the draft "The JWST Mass Mirage: C IV Outflows as a Driver of Overestimated Black Hole Masses at High Redshift" (Abu Seif, Nov 2026)

**Verdict in three lines.**
The one solid thing in the draft is real and already well known: C IV line widths track Mg II/Hβ widths poorly, and C IV virial masses are biased high when C IV is blueshifted. Everything built on top of it is wrong: the "same-object" three-line sample cannot exist, the headline 0.41 dex bias is not in the data, the z > 6 quasars and the Little Red Dots never had C IV-based masses, and the proposed X-ray estimator fails exactly for the objects it targets. The paper has been rebuilt (`paper.md` / `paper.docx`) around what the data actually show.

Every number below comes from `results/dr7_results.json` and `results/dr16q_results.json` (scripts in `code/`) or from a paper that was opened and read (`notes/litcheck_*.md`).

---

## A. Problems with the data and figures

### A1. A three-line "same-object sample of 12,423 quasars" is physically impossible from SDSS
Hβ 4861 Å and C IV 1549 Å are a factor 3.14 apart in wavelength. The SDSS spectrograph covers 3800–9200 Å (ratio 2.42); BOSS covers 3600–10400 Å (ratio 2.89). No SDSS or BOSS spectrum can contain both lines. In Shen et al. (2011) Hβ is fitted at z ≤ 0.89 and C IV at 1.5 ≤ z ≤ 4.95; the windows do not touch.

Real counts (DR7, quality cut FWHM ≥ 1000 km/s, FWHM/err ≥ 3):

| pair | N | redshift window |
|---|---|---|
| Hβ + Mg II | 12,514 | 0.35–0.89 |
| Mg II + C IV (non-BAL) | 21,810 | 1.50–2.25 |
| Hβ + C IV | **0** | — |
| all three | **0** | — |

The draft's 12,423 is within 1 % of the DR7 Hβ+Mg II overlap. So the draft measured Hβ–Mg II on one set of objects (z ≈ 0.35–0.9) and Mg II–C IV on a different set (z ≈ 1.5–2.25), then described both as one sample. The text must say so.

### A2. The "Hβ vs C IV" mass panel pairs different objects
Figure "Comparison between the Different Virial Estimators", lower-right panel, plots M(C IV) against M(Hβ) for several thousand points. Per A1 that set is empty. The only way to draw it is to place the Hβ column of one object next to the C IV column of another (e.g. two independently compacted spreadsheet columns). Its fit (y = −0.0756x + 9.86) is therefore meaningless, and the flat, structureless "C IV vs Mg II" clouds in the same figure and in the Excel scatter plot are suspect for the same reason. On real matched objects the C IV–Mg II relation is weak but clearly positive (DR7: Spearman ρ = 0.14, p ≈ 10⁻⁹²; the 1:1 line runs through the density ridge, see `figures/fig1_fwhm_pairs_dr7.png`). "No correlation" is not what the data say; "weak correlation, scatter-dominated" is.

### A3. The two-panel "Panel A / Panel B" figure (image9) is not a plot of catalog data
Panel (a) shows a curved ridge with uniform scatter reaching FWHM(Hβ) = 14,000 km/s at FWHM(Mg II) = 10,000 km/s, and panel (b) shows a uniformly filled rectangle. Neither resembles the catalog distributions (compare `figures/fig1_fwhm_pairs_dr7.png`). It looks like an illustration generated to match the quoted r-values and must not appear in a paper as data.

### A4. The headline "0.41 ± 0.12 dex overestimate" is not supported, and is calibration-dependent
Measured median offsets, log M(C IV, VP06) − log M(Mg II), DR7 matched objects (N = 21,323):

| Mg II calibration | median offset | NMAD scatter |
|---|---|---|
| Vestergaard & Osmer 2009 | +0.16 dex | 0.35 dex |
| Shen et al. 2011 (the catalog's fiducial) | −0.02 dex | 0.37 dex |

Nothing in the data produces 0.41 dex as a sample-wide offset, and the sign of the mean offset flips with the choice of Mg II zero point. What *is* robust is the dependence on C IV blueshift (A5) and the 0.35–0.40 dex object-to-object scatter, which matches Coatman et al. (2017: 0.40 dex) and Shen et al. (2008: 0.34 dex). The draft's statement also quotes "± 0.12" without saying what it is (standard error? scatter?).

### A5. The real result the draft missed: the C IV excess is a function of C IV blueshift
DR7, non-BAL, median log M(C IV) − log M(Mg II, VO09) by C IV blueshift (positive = blue):

| blueshift [km/s] | N | median Δlog M |
|---|---|---|
| −5000 to −500 (redshifted) | 909 | +0.09 |
| −500 to 0 | 2,689 | +0.08 |
| 0–500 | 7,321 | +0.11 |
| 500–1000 | 6,973 | +0.18 |
| 1000–1500 | 2,708 | +0.28 |
| 1500–2000 | 503 | +0.48 |
| 2000–3000 | 159 | +0.74 |
| > 3000 | 61 | +0.78 |

This is the Coatman et al. (2016, 2017) effect reproduced on 21,000 objects. Only 3.4 % of DR7 quasars have blueshift > 1500 km/s, so the *population* bias is small; the bias on an *individual* high-blueshift object is large. That distinction is the whole story and the draft inverts it.

### A6. The "Hβ versus Mg II" fit is misquoted and misread
Figure caption shows y = 0.5568x + 3.7637 with x = M(Hβ), y = M(Mg II). The text reports "log M(Hβ) = 0.59 log M(Mg II) + 3.76, scatter 0.18 dex" (axes swapped, coefficients changed). A slope of 0.56 is not a "tight linear relation"; it is what regression dilution does to two noisy estimators of the same quantity. On DR7 the ordinary least-squares slope is 0.83 (bisector slope 1.18), the median offset is +0.01 dex against the fiducial Mg II calibration, and the scatter is 0.24–0.25 dex (NMAD), not 0.18.

### A7. Redshift bins that SDSS cannot contain
- "FWHM(Hβ) at z > 1 … mean shifting to ~5,200 km/s": SDSS Hβ fits stop at z = 0.89 (DR7) / ≈1.1 (BOSS). The draft's z > 1 panel has ~120 objects and a histogram peak of 52; it is noise at the edge of the spectrograph, not a population.
- "Mg II … across all redshift bins up to z > 3": Mg II leaves the BOSS window at z ≈ 2.7. The z > 3 Mg II panel (~100 objects) cannot be clean.
- "C IV 3.5 < z < 6.5 bin contains ~350 objects": C IV leaves BOSS at z ≈ 5.7. DR7 has 2,003 non-BAL quasars with good C IV at 3.5 < z < 5.0.
- The Hβ caption says the peak is at ~6,000 km/s; the plotted peak is at ~4,500 km/s. Real DR7 medians: Hβ 4,100–4,500 km/s in every bin, Mg II 4,150–4,970 km/s, C IV 5,200–5,500 km/s.

### A8. "Double-peaked Mg II at 1 < z < 2" is a smoothing artefact
The DR7 histogram of FWHM(Mg II) at 1 < z < 2 (500 km/s bins, N ≈ 47,000) has one maximum (4,000–4,500 km/s) and declines monotonically; there is no second peak at 7,500 km/s. The draft's curve uses ~1,000 km/s points joined by a smoothing spline in Excel, which manufactures bumps.

### A9. "Mass evolution with redshift" figure
Log M versus z with a logarithmic fit (R² = 0.42) shows the flux limit of a magnitude-limited survey (Malmquist bias), as the draft itself says elsewhere. Calling it "evolution" contradicts the draft's own text. Drop it or label it as a selection envelope.

### A10. Methodology describes work that was not done
Section "Dataset Selection and Line Fitting" says the authors performed multi-component spectral fits with a power law, Fe II template and Gauss–Hermite profiles on > 10,000 S/N > 10 spectra. The figures are Excel charts of catalog columns. The paper must cite the catalog whose fits were used (Shen et al. 2011 for DR7; Wu & Shen 2022 for DR16Q) and describe *their* fitting, not claim it.

### A11. Statistics named inconsistently
Abstract and Results: "Pearson r ~ 0.65 / 0.08". Interpretation sections: "Spearman ρ = 0.65 / 0.08". Pick one and state what was correlated (linear FWHM or log FWHM). Real DR7 values: Hβ–Mg II Spearman 0.73, Pearson (linear) 0.69; C IV–Mg II Spearman 0.14, Pearson (linear) 0.17.

---

## B. Problems with the astrophysical argument

### B1. The z > 6 quasar masses were never C IV-based
The standard single-epoch estimator at z > 6 is Mg II from ground-based near-infrared spectroscopy, not C IV. Fan, Bañados & Simcoe (2023, ARA&A 61, 373, §4.1): "we only focus on results based on Mg II measurements, which are thought to be the most reliable until Hβ-derived masses are possible"; C IV estimates "are uncertain". Replacing C IV by Balmer lines therefore cannot "deflate" the published z > 6 masses, because those masses do not rest on C IV. (Object-by-object check: see §C.)

### B2. Little Red Dots never had C IV masses; most do not even show C IV
Every LRD black-hole mass in the literature is from broad Hα or Hβ measured by JWST itself (Kocevski et al. 2023; Matthee et al. 2024; Greene et al. 2024; Furtak et al. 2024; Kokorev et al. 2023). LRDs are a JWST discovery; there are no "legacy UV spectra" to compare against. C IV is typically absent (Rusakov et al. 2026: absent in 5/5 high-S/N LRDs with rest-UV coverage). "UNCOVER 29255" does not exist in any UNCOVER paper (real IDs with masses: 4286, 13123, 13821, 20466, 23608, 35488, 38108, 41225, 45924).

### B3. The "overmassive" tension is about M_BH/M_*, not about the mass scale being inflated by C IV
Pacucci et al. (2023) and Maiolino et al. (2024) argue the JWST black holes are overmassive *relative to their host stellar mass*, with masses from Balmer lines. If anything deflates LRD masses it is electron scattering of the Balmer lines (Rusakov et al. 2026, Nature 649, 574: intrinsic cores imply 10⁵–10⁷ M☉, ~2 dex lower), with a counter-example where a dynamical mass confirms the Balmer virial value (Abell 2744-QSO1: Juodžbalis et al. 2026, ≈5 × 10⁷ M☉). This is the live debate; C IV has no part in it.

### B4. The X-ray "outflow-free" estimator does not work for these objects
- The Γ–λ_Edd relation has an intrinsic scatter of 0.35–0.6 dex in λ_Edd (Shemmer et al. 2008; Risaliti et al. 2009; Trakhtenbrot et al. 2017 find it weak, r_s = 0.31, and absent for AGN with direct masses). Brightman et al. (2013): "not viable for single sources".
- Risaliti et al. (2009) find *no* Γ–λ_Edd correlation when λ_Edd is computed with C IV masses, which is itself evidence against C IV.
- LRDs are X-ray weak by 1–2 dex and undetected even in deep stacks (Yue et al. 2024; Ananna et al. 2024; Maiolino et al. 2025). There is no hard X-ray continuum to measure Γ from. The estimator fails precisely for the population the draft applies it to.
- "ε = 1.26 × 10³⁸ erg s⁻¹ M☉⁻¹" is the Eddington luminosity per solar mass (SDSS catalogs use 1.3 × 10³⁸), not an efficiency; ε is the standard symbol for radiative efficiency and will confuse readers.

### B5. Claims of "definitive", "complete lack of correlation", "mirage"
The data support "weak correlation" (ρ ≈ 0.14, highly significant) and "large, blueshift-dependent bias". They do not support "complete decoupling", "definitive evidence", or that the early-universe mass crisis is "largely a mirage".

---

## C. Object-by-object check of the ten "validation" targets
Sources in `notes/litcheck_z6_quasars.md` and `notes/litcheck_lrd.md`.

| draft object | exists? | original mass line | JWST mass | agree? |
|---|---|---|---|---|
| J0100+2802 | yes | Mg II (Wu+15) 1.2×10¹⁰ | Hβ 10.2 ± 0.1 dex (Eilers+23); 10.06 (Yue+24) | yes |
| J1120+0641 | yes | Mg II (Mortlock+11) 2.0×10⁹ | Hα 1.5×10⁹ (Bosman+24); Hβ 1.9×10⁹ (Marshall+25) | yes |
| J1148+5251 | yes | Mg II (Willott+03) ~3×10⁹ | Hβ 9.64 dex (Yue+24) | yes, 0.3 dex low |
| HSC J2239+0207 | yes | Mg II (Onoue+19) | none published | — |
| HSC J1205−0000 | yes | Mg II (Onoue+19) | none published | — |
| J2346+0529 | **no such quasar** | — | — | — |
| J0218+0007 | yes | Mg II (Yang+21) 0.61×10⁹ | Hβ 0.61×10⁹ (Yang+23) | identical |
| Abell 2744-QSO1 | yes | **Hβ, JWST** (Furtak+24) 4×10⁷ | dynamical 5×10⁷ (Juodžbalis+26) | yes |
| UNCOVER 29255 | **no such ID** | — | — | — |
| CEERS 1670 | yes | **Hα, JWST** (Kocevski+23) 1.3×10⁷ | — | — |

Not one of the ten had a C IV-based mass. The draft also cites "Marshall et al. 2023" and "Loiacono et al. 2024" for J1120; the first is about two other z ≈ 6.8 quasars, the second about PJ308−21 at z = 6.23.

---

## D. Editorial defects
- Empty "Reference" section; "[citation needed]" left in the text; no citations anywhere.
- Broken fragments: "pile a validation sample" (should be "We compile"), "reconciling early black holes wit / h standard seeding", "ABSTRACT: 1", "Introduction;", figure numbers missing ("Figure shows", "Figure 4/5" never defined), a stray "." paragraph.
- Results and Interpretation sections repeat the same numbers three times with different statistics (Pearson/Spearman) and different sample sizes (10,000 / 12,423).
- Units: "M☉" and "M_sun" and "Msun" mixed; "H-Beta"/"Hβ"/"H beta" mixed.

---

## E. What survives, and what the rebuilt paper says
1. Hβ and Mg II widths and masses agree (ρ ≈ 0.73; offset ≈ 0 dex; scatter ≈ 0.25 dex). **Kept, with real numbers.**
2. C IV widths correlate only weakly with Mg II widths (ρ ≈ 0.14). **Kept, re-worded from "none" to "weak".**
3. C IV masses scatter by 0.35–0.40 dex against Mg II and are biased high in proportion to C IV blueshift (+0.1 dex at zero blueshift, +0.5 dex at 1,500–2,000 km/s, +0.75 dex above 2,000 km/s). **New, replaces the 0.41 dex claim.**
4. Because z > 6 quasars have C IV blueshifts ~2.5× larger than at low z (Meyer et al. 2019; Schindler et al. 2020, median 1,848 km/s), a C IV mass at z > 6 *would* be biased by ~0.5 dex. **New, correctly framed: a warning, not a correction to published masses, which use Mg II.**
5. JWST Hβ masses agree with the earlier Mg II masses for the z > 6 quasars that have both (details §C). **Replaces the draft's unsupported claim of systematic deflation.**
6. The LRD mass question is about Balmer-line broadening physics and X-ray weakness, not C IV. **Draft's LRD section replaced.**
7. X-ray Γ–λ_Edd is a population tool with ≥ 0.4 dex scatter and no data for LRDs. **Proposal withdrawn; discussed as a limitation.**

---

## F. Self-review of the rebuilt paper

**Novelty: 4/10.** The C IV–blueshift bias is known (Coatman et al. 2016, 2017; Shen et al. 2008). What is new here is (i) the measurement on 21,000 (DR7) and 66,000 (DR16Q) matched objects instead of a few hundred, with the blueshift defined relative to Mg II exactly as in the z > 6 samples; (ii) the explicit demonstration that the two Mg II zero points (S10 vs VO09) move the "mean C IV bias" across zero, so the mean is not a physical quantity; (iii) the object-by-object audit showing that no published z > 6 or LRD mass is C IV-based, which closes the "mirage" hypothesis.

**What a referee will ask**
1. The S/N ≥ 3 and FWHM ≥ 1000 km s⁻¹ cuts: how do the correlations change at S/N ≥ 5 or 10? (Partly answered: DR16Q at SN_MEDIAN_ALL ≥ 10; see Table 2–3.)
2. DR16Q contains BAL quasars; cross-match with the DR16Q `BAL_PROB` column and show the C IV scatter with BALs removed.
3. The DR7 and DR16Q C IV–Mg II mass offsets differ (−0.02 vs −0.19 dex). The paper attributes this to the fainter DR16Q population via the luminosity trend; show it directly by matching the two samples in λL₁₃₅₀.
4. The blueshift relation is measured on the FWHM from the catalogs' multi-Gaussian fits; Park et al. (2017) prefer the line dispersion. Can the catalogs' Gaussian parameters be used to build σ_line?
5. Mg II itself has a mild FWHM-dependent virial factor (Mejía-Restrepo et al. 2018) and a luminosity-dependent bias; the Mg II reference is not perfect. Quantify how much of the C IV–Mg II scatter is Mg II's.
6. Table 5 is a compilation; a homogeneous re-fit of the public JWST spectra with one recipe would be stronger.

**Ranked next steps**
1. Cross-match DR16Q with the parent DR16Q catalog for `BAL_PROB` and `SN`, redo Tables 3–4 BAL-free (one evening; data public).
2. Luminosity-matched DR7 vs DR16Q comparison to pin the mean-offset difference on λL₁₃₅₀.
3. Fit the blueshift–bias relation with a two-parameter form and compare with the Coatman et al. (2017) correction factor; publish the SDSS-based correction with its uncertainty.
4. Apply that correction to the C IV masses of the Schindler et al. (2020) / Farina et al. (2022) z > 6 sample and compare with their Mg II masses and with the JWST Hβ masses now available (Yang et al. 2023; Yue et al. 2024; Marshall et al. 2025).
5. Drop the X-ray section to a paragraph unless a Γ–λ_Edd test on the DR16Q × eROSITA/XMM overlap is added.

**Not reproducible from the original draft**: the 12,423 three-line sample, the 0.41 ± 0.12 dex bias, the r = 0.08 figure, the Hβ–C IV panel, the Mg II double peak, the z > 1 Hβ and z > 3 Mg II bins, the 350-object 3.5 < z < 6.5 C IV bin, the quasar "J2346+0529", the LRD "UNCOVER 29255", and the three JWST citations for J1120+0641.
