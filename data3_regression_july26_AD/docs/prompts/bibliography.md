# Prompt for Sonnet — add a source-material bibliography to the modeling document

You are refining the LaTeX report `docs/reports/main.tex` in the repo
`/Users/adowling/DowlingLab/Membranes/data3_regression`. Add a proper
bibliography that catalogs **all source materials** — journal papers, group
manuscripts/reports, presentations, datasets, and prior code — and cite them from
the text. For any material that lives in the shared Box folder, put its
**relative file path in the entry's `note` field** so the path renders in the
bibliography (this is an explicit PI requirement: the team wants to see, in the
rendered PDF, exactly which file each source is).

## Toolchain (already verified on this machine)
- Use **biblatex + biber** (`biber` 2.19 and `pdflatex` are installed; `latexmk`
  auto-runs biber). Add to the preamble:
  `\usepackage[backend=biber, style=numeric, sorting=none]{biblatex}` and
  `\addbibresource{refs.bib}`. Put `\printbibliography` before `\end{document}`.
- Build with `latexmk -pdf main.tex` from `docs/reports/` and confirm a clean
  compile with **no undefined citations** and no biber errors. Run
  `latexmk -c` after to clean aux files (the repo `.gitignore` already ignores
  LaTeX build artifacts; do not commit `.bbl/.bcf/.aux`).

## `note`-field convention (do this for every Box/repo source)
Add the location as a `note`, e.g.:
`note = {Box: \path{.../Manuscript/Supplementary Information.docx}}` or
`note = {repo: \path{prior_analysis/B_regression_concentration_NAClanal.m}}`.
Box paths are relative to the shared-folder root
`Diafiltration Modeling and Experiments (Phillip and Dowling groups)/`; repo paths
relative to `data3_regression/`. Use `\path{}` (load `url`/`hyperref`, already
loaded) so underscores/spaces render safely.

## Entries to create in `refs.bib`

### Journal papers (in `references/` — open each PDF to confirm exact citation details)
1. **DATA1** — Ouimet, Brown, ... , Muetzel, Dowling, Phillip, "DATA: Diafiltration
   Apparatus for high-Throughput Analysis," *J. Membrane Sci.* **641** (2022)
   119743. File `references/ouimet_2022_jmembranesci.pdf`.
2. **DATA2** — Liu, Estrada, Ouimet, McClure, Latulippe, Phillip, Dowling,
   "Characterizing Transport Properties of Surface Charged Nanofiltration
   Membranes via Model-Based Data Analytics," *Ind. Eng. Chem. Res.* **64** (2025)
   12111–12130. File `references/liu_2025_iecr.pdf`.
3. **Soft Sensors** — Lilonfe, Estrada, Singh, Ouimet, Phillip, Dowling, "Soft
   Sensors Enable Real-Time Ion Concentration Measurements," SSRN preprint 6661382
   (not yet peer-reviewed). File `references/lilonfe_2025_ssrn.pdf`. Use `@misc`/`@unpublished`.
4. **Ink formulation** — Liu, Ouimet, Hoffman, Xu, Phillip, Dowling, "Optimization
   of Reactive Ink Formulation for Controlled Additive Manufacturing of Copolymer
   Membrane Functionalization," *ACS Appl. Mater. Interfaces* **16** (2024)
   59216–59233. File `references/liu_2024_acsami.pdf`.

### Group theory documents / reports (use `@unpublished` or `@techreport`)
Pull the exact set and their relative paths from **Table `tab:sources`** (keys
[S], [M], [SI], [Dn], [Tr], [BoE], [MS]) already in `main.tex`, plus:
- Active manuscript outline — Box `.../Manuscript/Outline.docx`.
- Process-model draft — Box `.../Manuscript/Diafiltration_Process_Model_v10.pdf`.
- `Theory - Multicomponent Donnan Exclusion.docx` — Box `.../Analyses/`.

### Datasets (`@dataset` or `@misc`)
From **Table `tab:nacldata`**: the analysis workbook (Box
`.../Analyses/BoE Analysis.xlsx`; repo twin `data/Rejection_Analysis.xlsx`)
and the raw campaign workbooks (Box `.../Diafiltration Campaign and Raw Data/
NF270_MC{2,3,4,5}.xlsx`). One entry per workbook is fine; note which sheets are used.

### Presentations (`@misc`, include date)
- `DATA3_single_salt_analysis_v2.pptx` (2026-05-13) — Box
  `.../Presentations/Dowling Lab Updates/`.
- `ICON_poster_final_v2.pdf` (2026-01-28) — same folder.
- Most recent Phillip-lab update: `Multicomponent calculations 3.18.26.pptx` —
  Box `.../Presentations/Phillip Lab Updates/`.
(These paths are under `Multicomponent Diafiltration/Diafiltration Campaign and Raw
Data/`.)

## Wire up citations in `main.tex`
Replace the informal parenthetical references with `\cite`s, at least:
- DATA1/DATA2 where the transport framework and workflow are introduced
  (§Scope, §Transport, §Statistical methods).
- The Soft-Sensor paper in §2.3 (conductivity map).
- The group theory docs where their equations are used (§Membrane transport model,
  §Donnan partitioning) — you may add a `\cite` column or footnote to
  Table `tab:sources` rather than disturbing its layout.
- The analyzed dataset in §2 and Table `tab:nacldata`.
Keep Tables `tab:sources`/`tab:nacldata` (they are the reconciliation/inventory
aids); the bibliography complements them as the full catalog.

## Guardrails
- **Only edit** `docs/reports/main.tex` and create `docs/reports/refs.bib`. Do
  not change any analysis code, `prior_analysis/`, or other docs.
- Do not alter the technical content, equations, tables, or numbers — this is a
  citations-only refinement.
- Preserve the existing `\todo`/`\note` red/blue flags.

## Report back
- The list of bib keys created and which sources still lack full citation details
  (so a human can fill them).
- Confirmation of a clean `latexmk` build with zero undefined citations, and the
  final PDF page count.
