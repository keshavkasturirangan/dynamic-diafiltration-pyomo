# Overleaf LaTeX Workspace

This folder contains a starter LaTeX project (`main.tex`) intended to sync with Overleaf using Git.

## Connect this folder to Overleaf (Git)

1. In Overleaf, open your project and enable Git access.
2. Copy the Overleaf Git URL from **Menu -> Git**.
3. In this repository, run:

```bash
cd overleaf-latex
git init
# Add your remote (replace URL)
git remote add overleaf <OVERLEAF_GIT_URL>
# First push
git add .
git commit -m "Initialize LaTeX project"
git push -u overleaf master
```

If this repository is already a Git repo (it is), you can instead add the Overleaf remote from the repo root and push only this folder.
