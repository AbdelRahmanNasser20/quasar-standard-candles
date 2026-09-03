# Instructions for an AI agent (or a person) pulling this research

This repository holds the Abu Seif / Nasser quasar-cosmology papers. Everything an agent needs is in plain files. No database, no login.

## One command to get everything
```bash
git clone https://github.com/AbdelRahmanNasser20/quasar-standard-candles.git
cd quasar-standard-candles
cat manifest.json          # machine-readable index of every paper, figure, dataset, script
```

## Layout (stable — do not rename)
```
paperN___short-title/       one folder per paper, N = paper number
  paper.pdf                 the compiled paper (read this first)
  paper.tex / paper.md      source
  figures/                  every figure, PDF + PNG
  results/*.json            every number quoted in the paper, machine-readable
  results/*.log             raw stdout of the analysis runs
  code/*.py                 the scripts that produced results/ from data/
  data/                     the input catalog(s), downloaded from public archives
  REVIEW.md                 (paper1 only) the audit of the original draft
manifest.json               index of all of the above
AGENT.md                    this file
```

## How to read a paper folder as an agent
1. Read `manifest.json` → pick the paper by `id` (e.g. `paper2`).
2. Read `paper.pdf` (or `paper.md`, same content) for the narrative.
3. Every number in the text is in `results/*.json`. Grep for it there rather than trusting the PDF.
4. To reproduce: `python3 -m venv qenv && ./qenv/bin/pip install numpy scipy emcee matplotlib && cd paperN___*/code && ../../qenv/bin/python analysisX.py` (run from the `data/` directory so the scripts find `lusso2020.tsv`, or symlink it).
5. Runtime for all scripts together is about 25 minutes on a laptop.

## Ground rules baked into these papers
- Nothing in `results/` is typed by hand. Every value is script output.
- The data file `data/lusso2020.tsv` is a verbatim VizieR download of Lusso et al. (2020), table 3 (`J/A+A/642/A150`). Do not edit it. Re-download with the URL in `manifest.json` if in doubt.
- If a number in the PDF and a number in `results/` disagree, the JSON wins and the PDF has a typo. Report it as an issue.

## For Nasser (dad): shortest path
- Open `paper2___slope-drift-degeneracy/paper.pdf`.
- Tell your AI: "Clone https://github.com/AbdelRahmanNasser20/quasar-standard-candles and read AGENT.md, then summarise paper2."
