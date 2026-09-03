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

## Where to go next (one line each; ranked in paper 2's REVIEW.md)
1. Host-galaxy light test at z < 0.7 — could kill or confirm the drift in an afternoon.
2. Forward model of the SDSS × XMM selection — the referee's first question.
3. Rerun on Risaliti+2026 and eROSITA catalogs when public.
4. Joint quasar + Pantheon+ fit with γ(z) free — your original Paper II, now with a real prerequisite.

## How to use this repo with your AI
Tell it:
> Clone https://github.com/AbdelRahmanNasser20/quasar-standard-candles, read AGENT.md and manifest.json, then read paper2___slope-drift-degeneracy/paper.md and REVIEW.md and summarise what was found and what to do next.

Every number in both papers is in `results/*.json`. Nothing was typed by hand. If a number in a PDF disagrees with the JSON, the JSON wins.
