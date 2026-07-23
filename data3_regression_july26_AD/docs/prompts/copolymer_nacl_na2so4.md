# Prompt for Sonnet — copolymer membranes: repeat the full workflow for NaCl (1:1) and Na2SO4 (1:2), a non-adsorbing-salt test

You are extending the DATA3 analysis to a **new membrane family** --- Laurianne
Estrada's COOH-functionalized copolymer membranes (her "spacer arm" manuscript) ---
and a **new salt, \ce{Na2SO4} (1:2, with the divalent \ce{SO4^{2-}} as the co-ion)**.
Repeat the **entire single-salt workflow** already applied to NF270 (data
processing/calibration, concentration-polarization prediction, Donnan regression,
AR(1) autocorrelation correction, uncertainty quantification) on four datasets, and
add a new report section. **Skip the exploratory $\chi(c)$ spline** for this dataset.
Work in `/Users/adowling/DowlingLab/Membranes/data3_regression`, env
`data3-regression`.

## Why this matters (the scientific framing --- keep it front and center)
Bill's ask: NaCl and \ce{Na2SO4} are **non-adsorbing salts** (unlike \ce{CaCl2}/
\ce{LaCl3}, whose multivalent \emph{counter}-ions adsorb and drove $|\chi|\to0$ with
systematic residual drift). So this is a **positive control**: the constant-charge
Donnan model is \emph{expected to fit well}, with well-identified, nonzero $|\chi|$
and no adsorption-style residual structure. Two things make it uniquely valuable:
1. **An independent check on $\chi$.** Estrada's SI reports the membrane fixed-charge
   density $\chi$ from \textbf{zeta-potential} measurements (Gouy--Chapman, SI
   Section~S2 / Table~S1) --- so for the first time in this program we can compare a
   \emph{regressed} $\chi$ against an \emph{independent physical} estimate. This is an
   accuracy test, not just precision.
2. **A new salt geometry.** \ce{Na2SO4} flips the valence: the divalent ion is now the
   \emph{co-ion}, so the Donnan polynomial differs from \ce{CaCl2}'s (see below).

## Data (already staged in the repo)
`data/copolymer/` contains:
- `Analysis of Select Data.xlsx` --- processed sheets, **standard 15-column schema
  identical to NF270** (B=$\cint$, D=$I$, F=$R$, H=$\cp$, J=$\Jw$, K=$B$, M=$t$,
  N=$\kappa_r$, O=$\kappa_p$; verified).
- Raw campaign workbooks `PEG0.1.xlsx`, `PEG6.5.xlsx`, `ALK1.2.xlsx` --- **same raw
  format as the NF270 campaign** (Time, Mass, Pressure, Temps, Conductivities, Vial
  Swap; Vial Data block with per-vial conductivity + ICP mg/L), so
  `analysis/preprocessing/preprocess_nacl.py::load_raw_sheet` parses them.
- `Spacer_arm_SI.docx` --- the manuscript SI. **Extract** the zeta-potential $\chi$
  (Table~S1) and the **water volume fraction $\phi_w$** (Section~S5 / Table~S3) for
  each membrane; you'll need $\phi_w$ for the Mackie--Meares membrane diffusivities
  and $\chi$ for the validation.

**The four focus datasets (all unadjusted pH~6):**
| Processed sheet | Membrane | Salt | Raw workbook / sheet |
|---|---|---|---|
| `PEG0.1 6.18.24NaCl` | PEG0.1 | NaCl (1:1) | `PEG0.1.xlsx` / `6.18.24NaCl` |
| `PEG6.5 6.27.24NaCl (2)` | PEG6.5 | NaCl (1:1) | `PEG6.5.xlsx` / `6.27.24NaCl (2)` |
| `PEG0.1 8.5.24Na2SO4` | PEG0.1 | \ce{Na2SO4} (1:2) | `PEG0.1.xlsx` / `8.5.24Na2SO4` |
| `ALK1.2 5.7.25Na2SO4` | ALK1.2 | \ce{Na2SO4} (1:2) | `ALK1.2.xlsx` / `5.7.25Na2SO4` |

Note the copolymer membranes run at **much lower water flux** than NF270
($\Jw\approx2\times10^{-6}$ vs.\ $\sim10^{-5}$ m/s) --- do not assume NF270 magnitudes.
Match SI $\chi$ to the specific membranes where the naming allows (e.g.\ PEG0 $\leftrightarrow$
PEG0.1); if a focus membrane has no SI $\chi$, say so.

## The \ce{Na2SO4} (1:2) Donnan model --- NEW (from Estrada's reworked equations)
Membrane negatively charged (COOH); co-ion $=\ce{SO4^{2-}}$ (divalent), counter-ion
$=\ce{Na+}$. Estrada's attached derivation gives the **co-ion cubic** (note the extra
$\tfrac14|\chi|^2$ term that \ce{CaCl2}'s cubic does \emph{not} have):
\begin{equation}
  (\cm{\co})^3 + |\chi|(\cm{\co})^2 + \tfrac14|\chi|^2\,\cm{\co} - K\,(\cs{\co})^3 = 0,
  \qquad K \coloneqq \dstar\,\delta^{\circ}_{\text{counter}},
\end{equation}
(lumping the non-ideality product into a single fitted $K\ge0$, exactly as \ce{CaCl2}/
\ce{LaCl3} did). Select the real positive root. Stoichiometry: **1 \ce{SO4^{2-}} per
\ce{Na2SO4}**, so $\cs{\co}=c$ (the reported salt concentration, i.e.\ $\cint$/$\cp$
directly --- \emph{not} $2c$); the counter-ion is $\cs{\text{ct}}=2c$. Ionic strength
$I=\tfrac12(1^2\cdot2c+2^2\cdot c)=3c$. Estrada also gives the counter-ion cubic
$(\cm{\text{ct}})^3-|\chi|(\cm{\text{ct}})^2-K(\cs{\text{ct}})^3=0$ (no $\tfrac14|\chi|^2$
term) --- keep it available for a rejection/counter-ion cross-check, but the primary
regression uses the **co-ion** root (matching the NF270 transformed-response
convention).

**Transformed response** (mirror NF270 §\ref{sec:reproMATLAB}/\ref{sec:cacl2}):
$\Delta c_{\mathrm m}=\Js\ell/\Dm$ with the co-ion solute flux $\Js=\cp\Jw$ (1
\ce{SO4} per salt, so no stoichiometric multiplier), $\ell$ the copolymer membrane
thickness (from the SI if given, else state the assumption), and the **ambipolar
(Nernst--Hartley) \ce{Na2SO4} membrane diffusivity**
$\Dm=\frac{(1+2)D_{\ce{Na}}D_{\ce{SO4}}}{D_{\ce{Na}}+2D_{\ce{SO4}}}$ evaluated at the
\emph{membrane}-phase $D$'s (Mackie--Meares from the copolymer $\phi_w$). Solution
diffusivities: $D_{\ce{Na}}=1.33\times10^{-9}$, $D_{\ce{SO4}}=1.065\times10^{-9}$
m$^2$/s. Prediction $=\cm{\co}(\cint)-\cm{\co}(\cp)$; fit $(|\chi|,K)$ **bounded
$\ge0$**. **NaCl (1:1)** reuses the existing quadratic co-ion model unchanged
(`analysis/step2_nacl`), just on the copolymer data and $\phi_w$.

## Workflow (repeat the NF270 pipeline; new code under `analysis/copolymer/`)
Reuse the existing machinery (`preprocessing/preprocess_nacl.py`,
`step2_nacl/regress_nacl.py`, `cacl2/regress_cacl2.py` as templates for the cubic
root-solve, `step8_autocorrelation/ar1_correction.py`, `step3_nacl_uncertainty`) ---
generalize, don't fork. Per dataset:
1. **Preprocessing/calibration from raw** --- $\Jw$ from the smoothed mass slope;
   bulk concentration from a per-experiment conductivity calibration fit from that
   run's vial (conductivity, ICP) pairs (ICP dilution fix + elemental molar masses:
   \ce{Na}=22.99; for \ce{Na2SO4} the ICP measures Na and/or S --- check which, and
   convert to salt mM accordingly). Extend the `Salt` config for a **divalent-co-ion
   salt** (\ce{Na2SO4}: $I=3c$, ambipolar $D$ above, co-ion $=$ salt conc); confirm
   NaCl behavior is unchanged.
2. **Concentration polarization** --- thin-film CP as the standing convention, on the
   authors' bulk $\cint$, with the \ce{Na2SO4}/NaCl solution ambipolar $D^{\mathrm s}$
   in the mass-transfer coefficient. Report the CP shift.
3. **Regression** --- bounded $(|\chi|,K)$ (or $(|\chi|,\dstar)$ for NaCl) with the
   correct co-ion polynomial; multi-start global %; scale-free fit metric.
4. **AR(1)** --- $\hat\phi$, $n_{\mathrm{eff}}$, GLS-whitened CI widening (Method A,
   $n-p$ df), per the corrected Step-8 convention (no double-count).
5. **Uncertainty** --- Wald + profile CIs; scipy, with a ParmEst cross-check if time
   permits. **Skip the $\chi(c)$ spline.**

**Adsorption diagnostics as the control test:** compute the residual-vs-$\cint$
sign-runs test and (optionally) the sliding-window effective-$|\chi|$ (Eq.~eq:effchi)
--- the \emph{expectation} here is the opposite of \ce{CaCl2}/\ce{LaCl3}: little/no
systematic residual structure and a well-identified, roughly constant $|\chi|$. State
clearly whether that expectation holds (it is the whole point of the "non-adsorbing"
control).

## Validation against the independent $\chi$ (the headline)
For each membrane, compare the **regressed $|\chi|$** to the **SI zeta-potential
$\chi$** (Table~S1). Agreement (same sign/order of magnitude) would be strong
independent support for the whole workflow; disagreement is itself an important
finding. Report the comparison explicitly per membrane.

## Laurianne's noise question (address it)
Several other copolymer runs (e.g.\ `ALK1.1 10.16.24NaCl`) have conductivity-probe
readings with large spikes/valleys that may interfere with regression. The four focus
datasets are the clean ones, so **do not regress the noisy ones here**, but add a
short, concrete methodological recommendation (a `\note` / discussion paragraph): e.g.
robust loss (Huber/soft-L1) in the fit; median/Hampel despiking or a rolling-median
prefilter on $\kappa$ before calibration; flagging air-bubble/vial-swap transients;
and that the AR(1) machinery already partially down-weights correlated excursions.
This connects to the raw-data-filtering step already on the project roadmap.

## Report (new section, after \ce{LaCl3})
Add a new `\section` **after §\ref{sec:lacl3}** (its own appendix too, matching the
per-salt standard), before the cross-salt synthesis. Cover: the copolymer membranes +
non-adsorbing-control framing; the \ce{Na2SO4} 1:2 co-ion cubic and how it differs
from \ce{CaCl2}; the four datasets; results table (regressed $|\chi|$, $K$/$\dstar$,
CIs, SSE, fit metric, AR(1) widening) with non-CP-vs-CP as elsewhere; the
**validation vs.\ zeta-potential $\chi$**; the adsorption-control verdict (does
constant-$\chi$ fit cleanly?); and the noise note. Publication-quality figures
(`analysis/plot_style.py`; math symbols, legends, **fonts $\ge$ body**, per the recent
figure conventions). Compile clean.

## Guardrails
- New code under `analysis/copolymer/`; new figures `docs/reports/figures/copolymer_*`.
  Reuse/generalize NF270 machinery; **do not change any NF270 (\ce{NaCl}/\ce{CaCl2}/
  \ce{LaCl3}) code, numbers, results sections, or the $\chi(c)$ spline work.**
- Edit `main.tex` only by adding the new section + appendix (and its entries); don't
  touch other sections, `refs.bib` (you may add one entry for Estrada's spacer-arm
  manuscript), or `data/` beyond reading `data/copolymer/`.
- Verify the processed-sheet CP status (the NF270 lesson: is col B already CP-corrected
  or bulk? use the $\cint/\kappa_r$ linearity check) before applying CP.
- Compile clean (`latexmk -pdf`: 0 undefined, 0 overfull; `latexmk -c`); rasterize +
  eyeball figures; report page count.

## Report back
Per dataset: regressed $(|\chi|,K/\dstar)\pm$CI, SSE, fit metric, multi-start %, CP
shift, AR(1) $\hat\phi$/$n_{\mathrm{eff}}$/widening; the **regressed-vs-zeta $\chi$
comparison** per membrane; the adsorption-control verdict (did constant-$\chi$ fit
cleanly, as expected for non-adsorbing salts?); the extracted SI $\phi_w$/$\chi$
values used; the \ce{Na2SO4} cubic implementation (and a root-solver sanity check);
the noise recommendation; sections/figures added; clean-compile confirmation; and
confirmation no NF270 results changed.
