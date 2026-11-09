# Where we are, and where to go next

Written for Nasser Abu Seif, 2026-09-03. Two pages, plain words.

## What happened
1. **Paper 1** (your draft, "No Evidence for Redshift Evolution…", 3 Aug 2026) was audited against the real Lusso+2020 catalog.
   - Real: γ = 0.664, β = 6.33 for the full sample.
   - Not real: δ = 0.317 (true 0.230), the four-bin table, N = 605 per bin, m = −0.002 ± 0.007. Three different sets of per-bin γ appear in the text. Details: `paper1___no-evolution-review/REVIEW.md`.
   - The honest measurement on the same data is a **2.7σ downward drift**, m = −0.042 ± 0.017.
2. **Paper 2** was built around what the data actually show.
   - Point 1: in flux-limited samples redshift and luminosity are the same axis (r = 0.98). "γ is flat in z" cannot be separated from "γ is flat in L" by binning. We show the cosmology-free test that can (within-shell curvature), and it rules out the luminosity explanation at 3.7σ.
   - Point 2: the drift survives every systematic we could test (selection floor, photon index, distance detrending, homogeneous subsample) but **vanishes above z = 0.7**. It is a low-z vs high-z contrast. This reconciles Risaliti/Lusso ("no evolution", samples start at z ≈ 0.7) with Khadka & Ratra and Li et al. ("inconsistent", anchors reach z ≈ 0.1).
   - Point 3: the cost. Best-fit drift = 23% distance bias at z = 2. For 2% distances you need σ_m ≈ 0.003, thirty times today's sample. A flat γ(z) at ±0.02 is necessary but not sufficient.
   - Self-review, weaknesses, referee questions: `paper2___slope-drift-degeneracy/REVIEW.md`.

## Update 2026-11-09 — Paper 3 (the C IV / "JWST Mass Mirage" draft)
3. **Paper 3** started from your November draft "The JWST Mass Mirage". It was audited the same way as Paper 1, against the real SDSS catalogs (Shen et al. 2011 DR7; Wu & Shen 2022 DR16Q) and the real JWST papers.
   - Real: Hβ and Mg II widths agree (ρ ≈ 0.7, masses equal to 0.01 dex, 0.25 dex scatter); C IV widths barely track Mg II (ρ ≈ 0.14–0.22); C IV masses scatter by 0.36 dex.
   - Not real: a three-line "same-object" SDSS sample (Hβ and C IV never share a spectrum; 12,423 ≈ the DR7 Hβ+Mg II count only); the 0.41 dex constant bias (the mean offset is −0.02 to +0.16 dex depending on the Mg II zero point); the Hβ–C IV panel (pairs different objects); the Mg II double peak; Hβ at z > 1 and Mg II at z > 3; the quasar "J2346+0529"; the LRD "UNCOVER 29255". Details: `paper3___civ-mass-bias/REVIEW.md`.
   - The honest result: the C IV bias is a *function of C IV blueshift* — negative at zero blueshift, +0.3 dex at 1,500–2,000 km/s, +0.4–0.6 dex above 2,000 km/s (DR7; DR16Q same trend, smaller amplitude for known reasons). That reproduces Coatman et al. (2017) on 100× more objects.
   - The big correction: no z > 6 quasar and no Little Red Dot ever had a C IV mass. z > 6 masses are Mg II, and JWST Hβ/Hα masses confirm them (J0100, J1120, J1148, J0218). LRD masses are JWST Balmer masses; the live debate there is electron scattering (Rusakov et al. 2026), not C IV. The X-ray Γ–λ_Edd estimator has ≥ 0.4 dex scatter and no data for X-ray-weak LRDs, so it was withdrawn.
   - Next steps (ranked in REVIEW.md §F): BAL-free DR16Q via the parent catalog; luminosity-matched DR7 vs DR16Q; publish the SDSS blueshift correction and apply it to the z > 6 C IV masses of Schindler et al. (2020).

## Where to go next (one line each; ranked in paper 2's REVIEW.md)
1. Host-galaxy light test at z < 0.7 — could kill or confirm the drift in an afternoon.
2. Forward model of the SDSS × XMM selection — the referee's first question.
3. Rerun on Risaliti+2026 and eROSITA catalogs when public.
4. Joint quasar + Pantheon+ fit with γ(z) free — your original Paper II, now with a real prerequisite.

## How to use this repo with your AI
Tell it:
> Clone https://github.com/AbdelRahmanNasser20/quasar-standard-candles, read AGENT.md and manifest.json, then read paper3___civ-mass-bias/paper.md and REVIEW.md and summarise what was found and what to do next.

Every number in both papers is in `results/*.json`. Nothing was typed by hand. If a number in a PDF disagrees with the JSON, the JSON wins.
