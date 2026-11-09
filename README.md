# quasar-standard-candles

Research papers on quasars: whether they can be used as standard candles for the Hubble constant, via the X-ray to UV luminosity relation
(log L_X = γ log L_UV + β), and how reliable their single-epoch black-hole masses are. Papers 1–2 are reproducible from the public Lusso et al. (2020) catalog included here; paper 3 from the public Shen et al. (2011) and Wu & Shen (2022) SDSS catalogs (download commands in `manifest.json`).

| Folder | What | Status |
|---|---|---|
| `paper1___no-evolution-review/` | Audit of the original draft "No Evidence for Redshift Evolution of the LX–LUV Relation" (Abu Seif 2026). Includes the original draft (`paper.pdf`) and the audit (`REVIEW.md`). | done 2026-09-03 |
| `paper3___civ-mass-bias/` | Rebuilt from the Nov 2026 draft "The JWST Mass Mirage": how much C IV λ1549 biases single-epoch black-hole masses (SDSS DR7 + DR16Q, 21k and 230k matched quasars), and why it does not explain JWST's overmassive early black holes. Includes the audit of the original draft (`REVIEW.md`). | draft v1, 2026-11-09 |
| `paper2___slope-drift-degeneracy/` | New paper: a cosmology-free, binning-free test of whether the slope γ drifts with redshift, and why "γ(z) is flat" cannot by itself prove standardizability in a flux-limited sample. | draft v1, 2026-09-03 |

Start with `WHERE-WE-ARE.md` for the story, `AGENT.md` if you are an AI agent or want the fastest path. `manifest.json` indexes every file.

Authors: Nasser Abu Seif (NRIAG, Cairo) and Abdel-Rahman Nasser. Analysis code and drafting assisted by Claude (Anthropic); every number is script output, nothing is hand-typed.

License: code MIT; text CC-BY-4.0; data © Lusso et al. 2020 via CDS/VizieR (cite the original).
