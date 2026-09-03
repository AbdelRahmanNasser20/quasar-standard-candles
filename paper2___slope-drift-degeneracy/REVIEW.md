# Review of Paper 2 — "Is a flat γ(z) enough?" (draft v1, 2026-09-03)

Self-review written the same day the draft was produced. Read this before trusting the paper.

## Verdict
- **Solid:** the degeneracy argument (z–L collinearity r = 0.98), the binning-free joint likelihood, the within-shell curvature test, the distance-bias formula and the σ_m ≈ 0.003 requirement. These are the parts a referee will not be able to knock down.
- **Moderate:** the drift itself. m = −0.042 ± 0.017 is a 2.7σ effect. AIC likes it, BIC does not. Binned significance swings from 2.2σ to 4.8σ with bin width. Say "evidence for", never "detection".
- **Fragile:** the drift is a contrast between z < 0.7 (364 objects, γ ≈ 0.61–0.69) and z > 0.7 (γ ≈ 0.57). Above z = 0.7 it is 1σ. Any story that removes the low-z shells removes the result.
- **Novelty:** 5/10. The framing and the requirement number are new as far as we found; the drift agrees with Khadka & Ratra (2021, 2022) and Li, Keeley & Shafieloo (2025).
- **Publishable?** After the fixes below and one pass by N.A.S., yes, as a letter-length paper (MNRAS Letters, A&A Letters, or RNAAS if cut to the requirement argument only).

## What a referee will ask, and the honest answer today
| Question | Where we stand |
|---|---|
| Is the truncated-Gaussian floor a real selection model? | No. It is a sharp floor at the 0th/1st/3rd percentile of each shell. It excludes only the simplest truncation bias. A forward model of the SDSS i < 20.2 and 4XMM flux limits is the definitive test and is not done. |
| Could the low-z excess in γ be host-galaxy light at 2500 Å? | Plausible and untested. Host light raises F_UV for faint low-z objects and steepens the apparent slope. Test: cut to L_UV > 10^30 in the low-z shells, or use L20's host-corrected fluxes if available. |
| Is the Γ_X term physics? | Partly methodological. L20 derive the 2 keV flux from band fluxes using Γ_X, so a residual Γ_X dependence is built in. We say so; we do not interpret c. |
| Does truncation bias the curvature test? | Yes, toward positive q, which is the direction that rejects curvature. The test is conservative but not clean. |
| Model L uses a cosmology. | Only to order the shells along the luminosity axis; tested for Ωm = 0.25–0.40. Still, say it louder. |
| The bias budget assumes a single slope extrapolated from a low-z anchor. | True. Risaliti–Lusso fit cosmology and the relation jointly; the bias then leaks into cosmological parameters rather than distances one-to-one. The order of magnitude stands; the mapping to H0 needs a proper joint fit. |
| Why not the Risaliti et al. 2026 sample? | Not public as of 2026-09-03. Rerun the moment it is. |
| Errors. | Early bootstrap errors on sub-fits were 2× too small (L-BFGS stopping early). Every quoted error in the paper is now MCMC or a full Nelder-Mead bootstrap. `analysis4.py` still prints the small bootstrap values; ignore them, use `analysis5/6` numbers. |
| Li et al. 2026 (A&A 706, A337) | Cited from its abstract only; the full text returned 403. Read it before submission. |
| Authorship/affiliations | Placeholder. N.A.S. must confirm affiliation, ORCID, and whether A.N. is an author or an acknowledgement. |

## Internal consistency checks done
- Full-sample γ = 0.664, β = 6.33, δ = 0.230 match L20's published values → data ingestion correct.
- Detrending to shell-centre distance changes m by < 0.004 for three Ωm → within-shell distance spread is not the cause.
- Two independent MCMC runs (different seeds) give m = −0.042 and −0.046, both ± 0.017.
- Group 5 alone (homogeneous SDSS × 4XMM) reproduces the drift (−0.060 ± 0.021).
- Distance-bias formula checked against the Δγ = 0.1 → ~0.6–0.75 mag rule of thumb used by L20.

## Where to go next (ranked)
1. **Host-light test at z < 0.7.** Cheapest, and it is the one thing that could kill or confirm the result. One afternoon.
2. **Forward-model the selection.** Simulate a parent population with the L20 relation, apply SDSS and 4XMM flux limits as functions of z, refit. If mocks with m = 0 reproduce the observed z < 0.7 excess, the paper's conclusion changes. Two to three days.
3. **Rerun on Risaliti et al. 2026 and eROSITA (Sacchi et al. 2025) catalogs** when public. If the z > 0.7 flatness holds there with 10× the objects, the σ_m requirement becomes the paper's main point.
4. **Joint QSO + Pantheon+ fit** with γ(z) free, to turn the bias budget into an actual H0 shift. This is the "Paper II" N.A.S. originally planned, now with a defensible prerequisite.
5. **Efron–Petrosian cross-check** against Dainotti/Lenart's (1+z)^k corrections on the same shells.

## Files
- `paper.pdf` / `paper.tex` / `paper.md` — the draft.
- `results/*.json` — every number. `results/*.log` — raw run output.
- `code/analysis3.py … analysis6.py` — reproduce everything (~25 min).
