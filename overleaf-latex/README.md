# Overleaf LaTeX Workspace

This folder contains a starter LaTeX project (`main.tex`) intended to sync with Overleaf using Git.

## Current Setup

- Overleaf remote name: `overleaf`
- Overleaf project URL: `https://git.overleaf.com/699d176f9b4bfb2fe06e101f`

## First-time authentication

Overleaf now uses Git authentication tokens.

1. Get your token from Overleaf Git settings.
2. Use remote URL format with `git@` user:

```bash
git remote set-url overleaf https://git@git.overleaf.com/699d176f9b4bfb2fe06e101f
```

3. Optional (recommended on macOS): store credentials in keychain.

```bash
git config --global credential.helper osxkeychain
```

## Recommended push workflow (works with existing Overleaf history)

Run from repo root:

```bash
cd /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo
git fetch overleaf master
git checkout -b overleaf-update overleaf/master
git show codex/unified_code:overleaf-latex/main.tex > main.tex
git show codex/unified_code:overleaf-latex/README.md > README.md
git add main.tex README.md
git commit -m "Update Overleaf files from local overleaf-latex"
git push overleaf HEAD:master
git checkout codex/unified_code
git branch -D overleaf-update
```

## Pull Overleaf changes back into `overleaf-latex/`

Run from repo root:

```bash
cd /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo
git fetch overleaf master
git show overleaf/master:main.tex > overleaf-latex/main.tex
git show overleaf/master:README.md > overleaf-latex/README.md
git add overleaf-latex/main.tex overleaf-latex/README.md
git commit -m "Sync Overleaf updates into overleaf-latex"
```

If Overleaf adds more files, copy those too during sync.
