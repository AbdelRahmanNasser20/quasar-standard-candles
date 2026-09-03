# Review: "No Evidence for Redshift Evolution of the LX–LUV Relation from 2421 Quasars" (Abu Seif, draft 3 Aug 2026)

## Verdict (read this, skip the rest)
- **Core idea: real, sound, and standard. Novelty 3/10.** It repeats the exact test Risaliti & Lusso (2015, 2019) and Lusso et al. (2020) already ran on this same catalog.
- **Headline number is not reproducible.** Paper: m = −0.002 ± 0.007 (flat). Real data, same method: m = −0.048 ± 0.017 (2.8σ downward drift). Full-sample γ = 0.664 and β = 6.33 are real; the per-bin table, the scatter δ = 0.317, the ±0.007 error, and the Monte Carlo section are fabricated or inconsistent.
- **Not publishable as-is.** Publishable after a rewrite that (a) reports the honest ~2–3σ drift, (b) engages the four groups who found evolution/non-standardizability in this same sample, (c) ships code + data. Realistic venue after that: MNRAS/A&A short paper or RNAAS.

## What the paper is trying to do (like you're 10)
- Quasars are super-bright black holes. Each one glows in UV light and in X-rays.
- Rule of thumb: brighter UV → brighter X-ray, but not 1-for-1. The "not 1-for-1" number is γ (about 0.6).
- If γ is the same for every quasar at every age of the universe, you can use quasars as "lightbulbs of known wattage." Compare known wattage to how dim they look → distance → how fast the universe expands (H0).
- This paper asks one question: **is γ the same for old (far) and young (near) quasars?** If yes, quasars are trustworthy rulers. If no, the ruler stretches and every distance is wrong.
- Paper says: yes, perfectly the same. Real data says: mostly the same, but with a small worrying drift that the paper hid.

## What I executed (all reproducible, files in this folder)
- Pulled the real Lusso+2020 catalog from VizieR: `lusso2020.tsv` (N = 2421, z = 0.009–7.54, median 1.30; **only 2391 lie in 0.1–5.0**, so the paper's "N = 2421 at 0.1 < z < 5" is impossible).
- `repro.py`: D'Agostini likelihood + emcee (32 walkers). Full-sample luminosity fit, 4 equal-count bins, 10-bin scan, 300-mock flux-limit test.
- `repro2.py`: proper narrow bins (Δlog(1+z) = 0.05) in flux–flux space, luminosity-space quartile bins, the paper's own bin edges.
- `gamma_vs_z.png`: the one figure that matters.

| Quantity | Paper claims | Real data (this rerun) | Status |
|---|---|---|---|
| γ full sample (lum. space, H0=70, Ωm=0.315) | 0.6634 ± 0.0055 | 0.664 ± 0.007 | ✅ real |
| β full sample | 6.33 ± 0.16 | 6.33 ± 0.22 | ✅ real |
| δ intrinsic scatter | 0.3167 ± 0.0046 | 0.230 ± 0.004 (Lusso+20: 0.24) | ❌ wrong |
| Per-bin γ (Table 1) | 0.673, 0.668, 0.640, 0.675 (±0.011) | 0.643, 0.594, 0.543, 0.589 (±0.02–0.03) | ❌ fabricated |
| Per-bin N (Table 1 edges) | 605 each | 665, 741, 675, 309 | ❌ fabricated |
| Per-bin β (Table 1) | 6.1–7.0 | 6.9–10.1 | ❌ fabricated |
| m = dγ/dz, narrow flux–flux bins | −0.002 ± 0.007 | **−0.048 ± 0.017** (2.8σ) | ❌ headline fails |
| m, luminosity-space quartiles | — | −0.031 ± 0.013 (2.3σ) | same direction |
| Constant-γ fit, narrow bins | "0.3σ from zero" | χ² = 18.9 / 9 dof (p ≈ 0.03), ⟨γ⟩ = 0.575 ± 0.012 | tension |
| Flux-limit bias on m (mocks) | < 0.01 | −0.008 ± 0.013 (300 crude mocks) | ✅ plausible order |
| Eq. (1) D_L inversion | — | algebra checks out | ✅ |

## Internal contradictions (proof the text was machine-loosened)
- Bin edges: §3.2 says 0.1/0.7/1.2/2.0/5.0; Table 1 says 0.1/0.9/1.5/2.3/5.0.
- Per-bin γ appears three times with three different sets of values (§3.2, Table 1, §4).
- Evolution parameter: m = −0.002 ± 0.007 (z-linear) in the abstract; k = 0.02 ± 0.04 in log(1+z) appears only in §5.4 and the Conclusion, never derived.
- Lusso+2020 quoted as m = −0.05 ± 0.02 in the intro and m = 0.01 ± 0.01 in Table 2.
- Bisogni+2021: "130 quasars" in intro, "–" for m in Table 2 (actual paper: ~30 XMM quasars at z ≈ 3).
- MCMC: 10,000 steps (§3) vs 5,000 steps (§5.2). MC recovery: m = 0.001 ± 0.007 (§3.1) vs −0.001 ± 0.008 (§5.2) vs γ_out = 0.660 ± 0.009 (Fig. 4).
- "Cross-matched with SDSS DR16Q × 4XMM-DR12": no. The catalog is Lusso+2020's own (DR14Q × 4XMM-DR9) taken verbatim from VizieR.
- Table 1 β ≈ 6 means the bins were fit in **luminosity** space (flux-space β would be ≈ −13), contradicting the "cosmology-free flux–flux" claim.
- "z = 1.85 bin low due to XMM band edge": rest 2 keV at z = 1.85 lands at 0.7 keV observed, inside the band. Not a band-edge effect.
- "Gelman–Rubin for emcee": wrong diagnostic; emcee walkers are correlated, use autocorrelation time.
- The σ_m = 0.007 is unattainable: 605 quasars per bin give σ_γ ≈ 0.02–0.03 (MCMC and bootstrap agree), so σ_m ≥ 0.013 for any 4-bin scheme.
- "Referee" is mentioned (§5.4) for an unsubmitted preprint.

## Literature the paper must engage (currently cites none)
- Khadka & Ratra 2020, 2021, 2022 (MNRAS 497/502/510): the Lusso+2020 sample is **not standardizable at z ≳ 1.5–1.7**; relation parameters shift with redshift and cosmology.
- Lenart, Dainotti, Bargiacchi et al. 2023 (ApJS) and Dainotti et al. 2022 (ApJ 931, 106): Efron–Petrosian analysis finds luminosities **do** evolve, (1+z)^k, and correct for it.
- Petrosian, Singal & Mutchnick 2022 (ApJL 935 L19): luminosity–luminosity correlations cannot by themselves fix the distance–redshift relation.
- Zajaček et al. 2024 (ApJ 961, 229): dust extinction biases quasar distances.
- Wang, Yang et al. 2022 (ApJ 941, 174 – Gaussian copula) find a redshift-evolutionary form fits better.
- Signorini et al. 2023 (A&A 676 A143) and Risaliti et al. 2023 review: the defense side and a real selection-function treatment to compare against.

## What an honest rewrite looks like
1. Title: "A mild redshift trend in the LX–LUV slope of the Lusso+2020 quasar sample: implications for quasar cosmology."
2. Method exactly as now, but narrow bins (Δlog(1+z) ≈ 0.05), report both z-linear and log(1+z) parametrizations, per-bin table with honest errors.
3. Result: γ = 0.575 ± 0.012 flux–flux average; m = −0.048 ± 0.017; constant-γ χ² p ≈ 0.03. State the 2σ upper limit |m| < 0.08 and translate that into the H0 bias it allows (the paper's own Δγ = 0.1 → 0.6 mag argument).
4. Selection-function forward model (real one: SDSS i < 20.2 + 4XMM flux limit), since the drift is the size a Malmquist effect could produce. This is the one place the paper could add something new.
5. Compare against Khadka & Ratra and Lenart+2023 directly. Agreeing with them from an independent fit is a legitimate contribution.
6. Data + code availability statement (this folder is the seed).
7. Drop Paper II promises, drop "proving," drop "referee."
