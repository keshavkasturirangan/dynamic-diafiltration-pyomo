# Overleaf LaTeX Workspace

This folder stores the local source of truth for the Overleaf-linked document assets in this repository.

## What `main.tex` addresses

`main.tex` is the SIAM UQ 2026 presentation-template document in `article` format (not Beamer). It currently includes:

- A title section for the intrusive UQ + MBDoE membrane characterization topic.
- A structured sequence of slide-like sections (title, motivation, forward model, noise model, inverse problem, sensitivities, FIM, identifiability, MBDoE objective, sequential workflow, case study, results, and limitations/future work).
- Core equations for ODE-constrained modeling, WLS/MLE estimation, sensitivity ODEs, Fisher Information Matrix, and sequential design criteria.
- Automatic slide numbering via a LaTeX counter macro so section numbering stays consistent when adding/removing sections.

## Current Setup

- Repo root: `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo`
- Overleaf remote name: `overleaf`
- Overleaf Git URL: `https://git.overleaf.com/699d176f9b4bfb2fe06e101f`
- GitHub branch: `codex/unified_code`

## First-time authentication

Overleaf requires Git authentication tokens.

```bash
git remote set-url overleaf https://git@git.overleaf.com/699d176f9b4bfb2fe06e101f
git config --global credential.helper osxkeychain
```

When prompted:
- Username: `git`
- Password: Overleaf Git authentication token

## A) Local changes -> Overleaf + GitHub

Use this when you edited files in `overleaf-latex/` locally.

```bash
cd /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo
git checkout codex/unified_code

# Commit local changes in repo
git add overleaf-latex
git commit -m "Update Overleaf document locally"
git push origin codex/unified_code

# Push latest files to Overleaf via fast-forward branch flow
git fetch overleaf master
git checkout -b overleaf-update overleaf/master
git show codex/unified_code:overleaf-latex/main.tex > main.tex
git show codex/unified_code:overleaf-latex/README.md > README.md
git add main.tex README.md
git commit -m "Sync local updates to Overleaf"
git push overleaf HEAD:master
git checkout codex/unified_code
git branch -D overleaf-update
```

## B) Overleaf changes -> Local repo + GitHub

Use this when you edited `main.tex` directly in Overleaf.

```bash
cd /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo
git checkout codex/unified_code
git fetch overleaf master

# Pull Overleaf versions into repo folder
git show overleaf/master:main.tex > overleaf-latex/main.tex
git show overleaf/master:README.md > overleaf-latex/README.md

# Commit and publish to GitHub
git add overleaf-latex/main.tex overleaf-latex/README.md
git commit -m "Sync updates from Overleaf"
git push origin codex/unified_code
```

If Overleaf adds new files (for example `.bib`, images, extra `.tex` files), also copy them with:

```bash
git show overleaf/master:<filename> > overleaf-latex/<filename>
```
