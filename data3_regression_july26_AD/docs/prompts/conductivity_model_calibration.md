# Prompt for Sonnet — First-principles conductivity calibration curves (Shedlovsky + MSA) to rule out the nonlinearity claim and adjudicate the calibration discrepancy

You are adding a **physics-based conductivity calibration** to the DATA3
diafiltration analysis, using the two validated conductivity models from the
Lilonfe/Estrada/Singh/Ouimet/Phillip/Dowling soft-sensor manuscript. Work in
`/Users/adowling/DowlingLab/Membranes/data3_regression`, env `data3-regression`
(`/opt/anaconda3/envs/data3-regression/bin/python`).

## Why this is decisive (read first)
The manuscript **"Soft Sensors Enable Real-Time Ion Concentration Measurements"**
(Lilonfe, Estrada, Singh, Ouimet, Phillip, Dowling; *Journal of Membrane Science*
submission MEMSCI-S-26-01544-2, 2026; code at
`github.com/dowlinglab/salts-conductivity-modelling`, local copy at
`/Users/adowling/DowlingLab/Membranes/salts-conductivity-modelling/`) is the ideal
authority for the DATA3 calibration question because it studies **the exact same
salts** (NaCl, CaCl2, LaCl3), with the **same conductivity sensor** (Innovative
Sensor Technology LFS 1107 — the same model as the DATA3 diafiltration inline
probes), and is **co-authored by our collaborators Phillip and Dowling**. Its
central results settle two open DATA3 items:

1. **The "nonlinear conductivity–concentration relationship" claim is refuted by the
   manuscript itself.** Its Figure 3 shows κ(c) is essentially **linear** for
   NaCl/CaCl2/LaCl3 over 0–100 mM, and both physics models match the measured data
   to **MAPE 1.0% (NaCl), 2.2% (CaCl2), 1.9% (LaCl3)** (Table 3). Reproducing κ(c)
   over the DATA3 ranges puts a first-principles, peer-reviewed number on this.
2. **The calibration-source discrepancy (discussion point 9)** can be adjudicated
   against a physics-based κ(c) from the same sensor + same salts. Crucially, the
   sensor obeys **χ = Φ·(I/V)** with Φ a *purely multiplicative* cell constant (Eq. 3
   of the manuscript) — so a probe miscalibration is a **scale error with NO additive
   offset**. Physical κ(c) passes through the origin; the DATA3 authors' large
   *negative* intercepts (CaCl2 −132.6, LaCl3 −302.6 µS/cm) therefore cannot come
   from the cell constant and are physically unjustifiable.

## The models (vendor verbatim; cite the manuscript)
Source: `salts-conductivity-modelling/conductivity.py`:
- `shedlovsky(conc, temp, epsilon, eta, lambda_0, a, z_1, z_2, lambda_0_cation, lambda_0_anion)`
  — single-salt Shedlovsky (manuscript **Eqs. 4–9**); `conc` in **M**, `eta` in
  **poise**, returns **mS/cm**.
- `msa(valency, diameters, diff_coeff, temp, eta, epsilon, lambda_0, salt_1_conc, salt_2_conc=, salt_3_conc=)`
  — MSA (manuscript **Eqs. 11–16**); `salt_1_conc` in **mM** (1 mM = 1 mol/m³),
  `eta` in **Pa·s**, per-ion lists in **[cation, anion]** order, returns **mS/cm**.

**Model validity ranges (drives model choice for DATA3):**
- **Shedlovsky is validated only to 0.1 N.** DATA3 exceeds this at the high end:
  NaCl →180 mM = **0.18 N**; CaCl2 →108 mM = **0.216 N** (c_eq = z·c). LaCl3 →28 mM
  = 0.084 N is within range.
- **MSA is validated to 4 M ionic strength** — covers every DATA3 range (max is
  CaCl2's I = 3·108 = 324 mM = 0.32 M). MSA is also the model the manuscript shows
  is superior for binary/ternary mixtures (Table 4: CaCl2–LaCl3 MSA 2.0% vs
  Shedlovsky 5.0%; Table 5 ternary: MSA 4.4% vs Shedlovsky 8.8%).
- **→ Use MSA as the primary DATA3 calibration; report Shedlovsky alongside but flag
  it as extrapolated above 0.1 N for NaCl and CaCl2.** (For single salts the two
  agree to ~0.1–0.3% MAPE, Table 3, so this mainly matters for rigor and for the
  eventual multicomponent step.)

**Setup:** copy `conductivity.py` **verbatim** into
`analysis/conductivity_models/conductivity.py` with a provenance header (manuscript
title/authors, JMS submission no., repo URL, "run unmodified"). Do not edit the model
code. Read `single_salt_conductivity.ipynb` for the call patterns/units and mirror
them.

## Parameters — manuscript Tables 1 & 2 (authoritative; T = 298.15 K, ε_r = 78.43)
These match the notebook; use the manuscript values and cite the tables.

**Shedlovsky (Table 1; η = 8.9e-3 poise; z_Cl = −1; conc in M):**
| Salt | Λ°_salt [cm²·S·equiv⁻¹] | Λ°_cation | Λ°_Cl | a [cm] | z_cation |
|---|---|---|---|---|---|
| NaCl | 126.45 | 50.10 | 76.35 | 4.0e-8 | 1 |
| CaCl2 | 135.85 | 59.50 | 76.35 | 4.3e-8 | 2 |
| LaCl3 | 145.90 | 69.70 | 76.35 | 4.9e-8 | 3 |

**MSA (Table 2; η = 0.89e-3 Pa·s; conc in mM; order [cation, Cl]):**
| Ion | λ°_mol [S·m²/mol] | σ (diameter) [m] | D° [m²/s] | z |
|---|---|---|---|---|
| Na⁺ | 50.08e-4 | 2.04e-10 | 1.33e-9 | +1 |
| Ca²⁺ | 59.47e-4 | 2.00e-10 | 0.79e-9 | +2 |
| La³⁺ | 69.70e-4 | 2.36e-10 | 0.62e-9 | +3 |
| Cl⁻ | 76.31e-4 | 3.62e-10 | 2.03e-9 | −1 |

## Tasks (driver: `analysis/conductivity_models/calibration_curves.py`)

1. **Generate κ(c)** from MSA (primary) and Shedlovsky (secondary) for NaCl, CaCl2,
   LaCl3 over each salt's diafiltration $\cint$ range (NaCl to ~180 mM, CaCl2 to
   ~108 mM, LaCl3 ~2–28 mM; also sample to c→0). Output µS/cm (probe unit) and mS/cm.
   Run at the experimental retentate temperature (~294–296 K from the raw sheets) for
   the Task-3 comparison and at 298.15 K to match the manuscript. Flag the points
   where Shedlovsky is extrapolated beyond 0.1 N.

2. **Linearity test (primary deliverable).** For each salt × model, fit **linear** and
   **quadratic** κ(c) over the diafiltration range; report slope, intercept, $R^2$,
   the quadratic curvature coefficient, and the **max deviation of the linear fit from
   the model curve** (absolute + % of range). State plainly whether κ(c) is linear to
   within the sensor's precision (**the manuscript's probe MAPE is ~0.93%**, and model
   MAPE is 1.0–2.3%). Expect: yes — smooth, very slightly concave-down, curvature
   negligible. **This, combined with the manuscript's Figure 3 for the same salts and
   sensor, is the definitive rebuttal of the nonlinearity guess.**

3. **Adjudicate the calibration discrepancy (discussion point 9)**, per salt. Overlay
   model κ(c) with (a) the DATA3 authors' inline calibration (NaCl κ=−63.71+76.69·c;
   CaCl2 κ=−132.6+152.1·c; LaCl3 κ=−302.6+216.3·c, µS/cm vs mM) and (b) the Phase-1
   vial-ICP fit (e.g. CaCl2 MC3 +96.4+187.6·c). Report:
   - **Slope** agreement (a cell-constant Φ difference is a pure multiplicative scale,
     so slope *can* legitimately differ between probes — compare in relative terms).
   - **Intercept:** since χ = Φ·(I/V) has no additive term and physical κ(c) → 0 as
     c → 0, the physically admissible intercept is ≈0 (a concave-down κ(c) fit over a
     high-c window gives a *small positive* intercept, never a large negative one).
     Use this to argue the authors' **negative** intercepts are a fitting artifact and
     that a through-origin (or small-positive-intercept) linear calibration is the
     physically correct form.
   - **Caveat — known multivalent underprediction:** the manuscript reports the models
     slightly *underpredict* CaCl2 and (more) LaCl3 conductivity (~2% MAPE; dielectric
     decrement / hydration effects on multivalent ions), so treat the model as the
     *shape/linearity/intercept* authority, not a to-the-percent absolute calibration
     for the multivalent salts. Say which empirical calibration the physics supports
     and by how much, with this caveat stated.

4. **Reproduce the manuscript's single-salt result** as an implementation check: using
   the vendored code + Table 1/2 params + the manuscript's own experimental data
   (`salts-conductivity-modelling/Data/{NaCl,CaCl2,LaCl3}.csv`), reproduce the Table 3
   MAPEs (NaCl 1.0%, CaCl2 2.2/2.3%, LaCl3 1.9/2.2% for Shedlovsky/MSA). This confirms
   the copy + parameters before drawing conclusions, and gives an **experimental,
   same-sensor κ(c) reference** for the three salts to anchor Task 3.

5. **Figure** (`docs/reports/figures/`, plot_style compliant — boxed legends, bold
   math, PNG+PDF, rasterize-verify): model κ(c) (MSA primary, Shedlovsky secondary)
   for the three salts over the DATA3 ranges, with the authors'-inline and vial-ICP
   calibrations overlaid, plus a residual/linearity panel. One clear multi-panel
   figure.

## Report write-up
Add a compact subsection to `sec:preproc` — e.g. "First-principles conductivity models
(Shedlovsky & MSA)" — that: (i) introduces the two models with their manuscript
equation numbers (Shedlovsky Eqs. 4–9, extended Eq. 10; MSA Eqs. 11–16) and validation
pedigree (same salts + same LFS 1107 sensor, MAPE 1.0–2.3% vs measured, ~0.93% probe
precision — **cite the manuscript**); (ii) presents the **linearity verdict** as the
definitive, first-principles + peer-reviewed rebuttal of the nonlinearity guess; (iii)
adjudicates discussion point 9 using the cell-constant/zero-intercept argument, with
the multivalent-underprediction caveat. Then update the Phase-1 validation paragraph
and **discussion point 9** to cite this as first-principles confirmation. Note that the
**MSA multi-salt model (Eqs. 11–16, 20–23) is the tool for the eventual multicomponent
step (Step 14)** — vendoring it here is a building block, and the manuscript's
Hunter–Reiner model-discrimination design (Eqs. 20–23, its Fig. 5) is directly
reusable for planning multi-salt diafiltration experiments.

Add **one** `refs.bib` entry for the manuscript (this is the exception to "don't touch
refs.bib"); if a soft-sensor citation already exists, reuse/update it.

## Guardrails
- Vendor `conductivity.py` verbatim (provenance header only); do not modify model math.
- New code confined to `analysis/conductivity_models/`; do **not** change the
  regression drivers, results sections/numbers, `data/`, or `prior_analysis/`.
- Edit `main.tex` only in `sec:preproc` (new subsection + validation/point-9 updates);
  one allowed `refs.bib` addition.
- Compile clean (`latexmk -pdf`: 0 undefined, 0 overfull; `latexmk -c` after); report
  page count. Rasterize + eyeball the new figure.

## Report back
- The linearity test table (per salt × model: slope, intercept, $R^2$, curvature, max
  linear-fit deviation) and the plain-language verdict on nonlinearity (tied to the
  manuscript's Fig. 3 / Table 3).
- The calibration adjudication per salt (model vs authors'-inline vs vial-ICP: relative
  slope, intercept sign/magnitude, which the physics supports) with the cell-constant
  reasoning and the multivalent-underprediction caveat.
- The Task-4 reproduction check (do you recover the manuscript's Table 3 MAPEs?).
- Sections/figures/refs updated, and a clean-compile confirmation (0 undefined, 0
  overfull, page count).
