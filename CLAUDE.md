# quasar-standard-candles — CLAUDE.md

## How to answer
1-3 lines first, details nested under. Plain words. No preambles, no closing summaries.

## What this is
Papers by Nasser Abu Seif (NRIAG Cairo, Abdel's dad) + Abdel on whether quasars are standard candles via log L_X = γ log L_UV + β. Public repo; his AI clones it. Story: `WHERE-WE-ARE.md`. Index: `manifest.json`. Agent path: `AGENT.md`.

## HARD RULES
- **No hand-typed numbers.** Every value in a paper comes from `results/*.json` written by `code/*.py`. JSON beats PDF.
- **`data/lusso2020.tsv` is a verbatim VizieR download** (J/A+A/642/A150 table3). Never edit; re-download from the URL in `manifest.json`.
- **Folder rule `paperN___short-title/`** with `paper.pdf`, `paper.tex`, `paper.md`, `REVIEW.md`, `code/`, `data/`, `results/`, `figures/`. New paper = `paper3___…`, add to `manifest.json` and `README.md`.
- **Every paper ships with a `REVIEW.md`**: verdict, novelty /10, referee questions, ranked next steps. Doubt everything; state what is not reproducible.
- **Errors from MCMC or full bootstrap only.** L-BFGS bootstrap under-estimated σ_m by 2× (see paper2 REVIEW). Cosmology-free claims must be tested at Ωm = 0.25/0.315/0.40.
- **Cite only verified refs** (arXiv API / ADS). Li et al. 2026 A&A 706 A337 is cited from abstract only.
- Never reuse the per-bin table from dad's original draft (`paper1…/paper.pdf`); it does not reproduce.

## Key facts (2026-09-03)
- Paper 1 = audit of dad's draft: γ = 0.664, β = 6.33 real; δ = 0.317 and m = −0.002 ± 0.007 not. Real drift m = −0.042 ± 0.017 (2.7σ), vanishes for z > 0.7. Novelty 3/10.
- Paper 2 draft v1: z–L collinearity r = 0.98; within-shell curvature rules out luminosity explanation 3.7σ; 2% distance at z = 2 needs σ_m ≈ 0.003. Novelty 5/10. Next: host-light test at z < 0.7, selection forward model, rerun on Risaliti+26 / eROSITA when public.

## Commands
- Env: `python3 -m venv qenv && ./qenv/bin/pip install numpy scipy emcee matplotlib`
- Reproduce paper 2 (~25 min): `cd paper2___slope-drift-degeneracy/data && for s in 3 4 5 6; do ../../qenv/bin/python ../code/analysis$s.py; done`
- Build PDF: `cd paper2___slope-drift-degeneracy && pdflatex paper.tex && pdflatex paper.tex && pandoc paper.tex -t gfm --wrap=none -o paper.md`
- Local venv with deps lives in the Claude job tmp dir, not here; `qenv/` is gitignored.
- Send to dad: Messenger thread "Nasser Abu-Seif" (browser). Commits end with the Claude co-author trailer.
