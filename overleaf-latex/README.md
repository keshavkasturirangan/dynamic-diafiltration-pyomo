# Overleaf LaTeX Workspace

This folder contains a starter LaTeX project (`main.tex`) intended to sync with Overleaf using Git.

## Current Setup

- Overleaf remote name: `overleaf`
- Overleaf project URL: `https://git.overleaf.com/699d176f9b4bfb2fe06e101f`

## First-time authentication

If prompted, use your Overleaf Git credentials from **Overleaf -> Menu -> Git**.

```bash
git config --global credential.helper osxkeychain
```

## Push local changes to Overleaf

Run from repo root:

```bash
cd /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo
git subtree push --prefix overleaf-latex overleaf master
```

## Pull Overleaf changes back into this folder

Run from repo root:

```bash
cd /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo
git subtree pull --prefix overleaf-latex overleaf master --squash
```

This imports the Overleaf project into `overleaf-latex/` and records it as a single merge commit.
