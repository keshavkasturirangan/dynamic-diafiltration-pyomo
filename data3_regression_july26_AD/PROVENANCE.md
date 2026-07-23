# Provenance — data3_regression_july26_AD

Snapshot of Alex Dowling's (PI) DATA3 regression-analysis repo, pulled into this
branch so his model/statistics work travels with our codebase.

| What | Value |
|---|---|
| Source | https://github.com/dowlinglab/data3_regression |
| Branch | `main` |
| Commit | `5b8b1b7c11e281277d1e7d515d382e8a4690d8e3` |
| Commit date | 2026-07-17 ("Added feedback from meeting with collaborators") |
| Imported | 2026-07-22 (`.git` stripped — this is a plain snapshot, not a submodule) |

## What it is
An extended **Nernst–Planck / Donnan-equilibrium** membrane transport model whose
regressed quantities are the membrane **fixed charge χ** and a thermodynamic
**partition factor δ\*** — physics NOT in our DATA2-style Spiegler–Kedem (Lp, B, σ)
workflow — plus the statistical machinery (autocorrelation-corrected confidence
regions, profile likelihood, bootstrap, raw-measurement regression).

## Key entry points
- `DATA3_regression_notes_AD_2026-07-07.pdf` — Alex's 62-page write-up
  (was `~/Downloads/main.pdf`; compiled from `docs/reports/main.tex`, which is
  10 days older than the repo HEAD — the repo may contain post-meeting updates
  not reflected in the PDF).
- `README.md` — repo orientation.
- `analysis/` — the regression scripts; `data/` — processed inputs;
  `references/` — DATA1/DATA2 + related papers.

## To refresh the snapshot
```bash
git clone --depth 1 --branch main https://github.com/dowlinglab/data3_regression /tmp/d3r \
  && rsync -a --delete --exclude=.git --exclude=PROVENANCE.md --exclude=DATA3_regression_notes_AD_2026-07-07.pdf /tmp/d3r/ data3_regression_july26_AD/
```
(then update the commit SHA above).
