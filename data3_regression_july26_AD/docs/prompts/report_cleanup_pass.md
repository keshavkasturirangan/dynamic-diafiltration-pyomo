# Prompt for Sonnet — full report cleanup + synthesis pass

You are doing a whole-document cleanup/consistency pass on the report
`docs/reports/main.tex` (`/Users/adowling/DowlingLab/Membranes/data3_regression`,
env `data3-regression`, build with `latexmk -pdf` in `docs/reports/`). The report
has grown to ~41 pp through many incremental edits by several sessions; it needs a
coherence pass plus two substantive additions. **This is a
documentation/presentation task — do NOT change any analysis results, numbers, or
the model mathematics.** Verify numbers against the code/tables and fix only
stale/inconsistent *presentation*. Result of an Opus audit — the specific items
below are real findings, not hypotheticals.

> **Run this AFTER the LaCl3 (Step 10) session has been committed** — it edits
> `main.tex` too, and the cross-salt synthesis (Task B) and agenda (Task A) must
> include LaCl3. Use the finalized `analysis/plot_style.py` for any figures.

## Guiding principle
The report is the team's **handoff + discussion document**: comprehensive,
self-contained, internally consistent, and honest about open questions. Favor
clarity and currency; keep the exploratory detail.

## Task A — consolidate the open decisions into a single discussion agenda
The ~19 `\todo`/`\note` flags and the §1.2 "Modeling questions" list are the raw
material, but §1.2 predates the newest results and there is no single actionable
agenda. **Build a consolidated "Open decisions for the collaborative meeting"** —
either by upgrading `sec:questions` (§1.2) into a table (columns: *item* |
*where flagged (\S ref)* | *the question* | *candidate resolution / options*) or a
crisp enumerated agenda that a reader can walk down in a meeting. It must **cover
every open item**, including ones §1.2 is currently missing:
- Sign/δ conventions; the factor-of-2 (bug); δ° for $z\ne1$; fixed-vs-estimated
  ($\ell,\phi_w,\Dm$); single→multi-salt transferability of $\chi$/$B$; transport
  mechanism reduction (M1–M6). *(already in §1.2)*
- Preprocessing: $A_{\mathrm m}$/$\rho$ inferred; $\Jw$ smoothing; calibration
  source **and the vial-vs-inline instrument/scale question**; ICP molar mass;
  CP-vs-non-CP (and that §2/§3 used non-CP); the "(2)"-sheet mostly-pasted column.
- **AR(1) (`sec:naclAR1`):** which correction convention is the headline — Method A
  (GLS/whitening, $n-p$) vs. Method B ($n_{\mathrm{eff}}$ heuristic); and AR(1)
  insufficiency (whitened ACF still structured).
- **Raw-measurement reformulation (`sec:altform`):** the estimates shift $\sim$2–3×,
  so is the transformed $|\chi|\approx\SI{50}{mM}$ robust? plus the non-stationary
  noise / need for CP+temperature+state-space model.
- **Concentration-dependent surface charge** for the multivalent salts (CaCl2 $\to$
  LaCl3): a data-driven $\chi(c)$/$\chi(t)$ (spline) is the proposed model extension.
- **Multi-salt identifiability:** conductivity is one scalar but multiple salts
  contribute (needs MSA multi-salt conductivity + ICP) — the barrier to the
  multicomponent goal.
Keep the in-place `\todo`/`\note` flags, but make sure none contradict the
consolidated agenda (fix any that do).

## Task B — add a document-level cross-salt synthesis
There is currently no place that ties the 1:1 / 2:1 / 3:1 story together (discussion
is NaCl-scoped in §2.5 plus per-salt verdicts). **Add a new top-level section**
(e.g. `\section{Cross-salt synthesis and discussion}` placed after §4 LaCl3 and
before §5 "Analysis plan with Claude Code") that synthesizes:
- The **valence trend**: NaCl (1:1) is a clean, well-identified Donnan fit
  ($|\chi|\approx\SI{50}{mM}$ corrected, $\dstar\approx0.33$); \ce{CaCl2} (2:1) and
  \ce{LaCl3} (3:1) drive $|\chi|\to0$ with strong structured residuals and
  effective-$\chi$ drift — i.e. a constant-charge Donnan model fails progressively
  for multivalent counter-ions, consistent with concentration-dependent
  adsorption. (Pull the actual LaCl3 numbers from `sec:lacl3` once Step 10 lands.)
- The **methodological cross-cuts**: identifiability depends on concentration range
  (MC2 vs MC5); i.i.d. confidence regions are $\sim$2.7× too narrow once
  autocorrelation is accounted for; and the estimate is sensitive to whether the
  regression is posed on transformed or raw measurements (`sec:altform`).
- What this implies for next steps (the $\chi(c)$ model; the multicomponent goal).
Keep it a synthesis (cross-reference the per-salt sections; don't restate every
number). Consider whether the NaCl-scoped §2.5 recommendations should stay in §2 or
fold selected cross-cutting points up into this new section.

## Task C — currency/consistency sweep (fix these specific stale spots)
- **§2.5 (`sec:naclDisc`) "Roadmap (i)–(v)" bullet** lists items that are now DONE
  (preprocessing, CaCl2, AR(1), raw reformulation, LaCl3) as if future work —
  update to the *current* roadmap (raw-data filtering; data-driven $\chi(c)$ spline;
  AR analysis for CaCl2/LaCl3 after the $\chi(c)$ model; multicomponent deferred
  pending feedback), or move the roadmap into Task B's new section.
- **§2.5 "Uncertainty" bullet** predates AR(1) ("Wald intervals adequate...
  precisions optimistic") — update to cite the `sec:naclAR1` result (i.i.d. CIs
  $\sim$2.7× too narrow) and the `sec:altform` formulation-sensitivity.
- **§5 (`sec:claude`)** "Implementation (Sonnet)" paragraph only lists Steps 2–3;
  update to reflect the full scope actually done (preprocessing pipeline, CaCl2,
  figure-compliance audit, AR(1)+both-ways correction, raw-measurement, LaCl3).
  Add to "Lessons learned": the Opus audit catching the AR(1) profile
  **double-counting** (a real statistical error caught before it reached
  collaborators), and the iterative report-refinement / figure-compliance rounds.
- **Abstract**: it still frames CaCl2/LaCl3 as future ("extend to…") and doesn't
  mention the AR(1) or raw-measurement findings — update to reflect that all three
  salts + those analyses are done.
- **Global number/reference sweep:** verify every quoted number matches its
  table/figure/code (spot-check |χ| values, n, φ, n_eff, fold-widenings,
  pseudo-R², calibration slopes); every `\ref`/`\eqref` resolves and points to the
  right target; terminology is consistent ($\dstar$ vs $\delta^\circ$, $|\chi|$,
  co-/counter-ion, "corrected model", $n_{\mathrm{eff}}$). Fix any residual stale
  claims (e.g. leftover "$8\times$" framing outside the deliberate "earlier-draft
  artifact" explanations).

## Task D — redundancy & flow (lighter touch)
- The adsorption hypothesis, the autocorrelation caveat, CP-vs-non-CP, and the
  narrow-range identifiability point each appear in several places. Self-contained
  sections are fine, but remove/fix any *contradictory or superseded* duplicates
  and cross-reference rather than restate where it reads repetitively.
- *(Optional, note don't force):* `sec:altform` (raw-measurement, a deep
  methodological aside on the primary dataset) currently sits between `sec:naclAR1`
  and the multi-dataset workflow `sec:naclAll`, interrupting the NaCl results
  progression. Consider whether the report reads better with the two "advanced
  single-dataset explorations" (AR(1), raw-measurement) grouped *after* the
  multi-dataset results. Only reorder if it's clearly cleaner and you can do it
  without breaking references.

## Task E — capture the nuanced analytical points (synthesized from the Opus audits of Steps 1–9)
These are the subtle scientific/statistical points that came up while auditing each
step. **Two jobs:** (1) verify each is stated in its home section (add a crisp
sentence where missing — most are present but confirm), and (2) **elevate the
recurring, cross-cutting ones into the new synthesis section (Task B)** as explicit
methodological takeaways, since they matter across all salts. Do not restate
numbers you can't verify against the code/tables.

### Per-section nuances (confirm each is present; add if missing)
- **§2.1 (`sec:reproMATLAB`):** the missing factor-of-2 is an *exact
  reparameterization* — $f(2\beta_1,4\beta_2)=2f(\beta_1,\beta_2)$, so it leaves the
  SSE surface (and fit quality) unchanged and only rescales $|\chi|\!\to\!2\times$,
  $\dstar\!\to\!4\times$; it is a bug in parameter *values*, not in the fit. The
  corrected $|\chi|\approx\SI{50}{mM}$ is closer to the independent BoE
  $\chi\approx\SI{44}{mM}$ than the uncorrected $\approx\SI{25}{mM}$ (weak physical
  support for the corrected form). $\ell$ and $\Dm$ are fixed but *confounded* with
  $\dstar$ via $B=\Dm H/\ell$.
- **§2.2 (`sec:naclUQ`):** the MATLAB non-convergence Bill saw is a *genuine local
  minimum at $|\chi|\to0$* that only reflective/trust-region-reflective (`dogbox`)
  methods fall into — **not** ill-conditioning (`lm`/`trf`/Ipopt reach the global
  optimum 100%); hence the fix is Levenberg–Marquardt, not rescaling. NaCl is
  well-identified and nearly linear near the optimum (nonlinear region $\approx$
  Wald ellipse). scipy$\leftrightarrow$ParmEst agree to 6–7 sig figs (mutual
  validation).
- **§2.4 (`sec:naclAll`):** MC2-vs-MC5 divergence is most plausibly a
  *concentration-range identifiability* effect (MC2's narrow $\lesssim\SI{44}{mM}$
  range → weak Donnan curvature → weak $|\chi|$), i.e. a less-informative
  experiment, not necessarily a different membrane. The 6 reconstructed runs are
  *reconstruction-dominated* (scatter, disagree among themselves) → the reliable
  conclusion rests on the 3 processed runs. **A tight CI $\neq$ accuracy** (e.g.
  reconstructed MC5 07.23.24 is precisely identified yet biased).
- **`sec:preproc`:** the ICP `mg/L` is the *diluted* sample (dilution-factor fix;
  elemental \ce{Na+} molar mass); the "(2)" sheet's $B$ column is *mostly pasted*
  (22/747 rows are the live formula); the two §2 datasets *lack per-vial
  conductivity* (calibration not independently checkable); and — flag if not
  already — the *vial-fit calibration differs from the inline-applicable embedded
  curve in slope AND intercept ($\sim$+\SI{175}{\micro S/cm})*, biasing the
  low-concentration end (a likely driver of reconstruction noise; open question
  whether vial and inline conductivity share an instrument/scale). $A_{\mathrm m}$,
  $\rho$ inferred; $\Jw$ smoothing $\sim$3.8\%; CP shift 12–22\%.
- **§3 (`sec:cacl2`):** $|\chi|\to0$ is *sharper* than "smaller than NaCl"; the
  effective-$\chi$ drift is *non-monotonic (inverted-U)*, complicating a simple
  monotonic-adsorption story; and that drift is a *windowed refit with $K$ fixed at
  the global $\chi\approx0$* — a **misfit diagnostic, not a direct $\chi(c)$
  measurement** (could be partly a fixed-$K$/CP artifact — both readings already
  flagged). MC5 06.27.26 fits worse ($R^2$ 0.84 vs 0.94) — possible fouling/noise
  state.
- **§4 (`sec:lacl3`):** confirm/add — LaCl3 drives $|\chi|$ **exactly onto its
  $\ge0$ bound** (a *boundary* solution, not an interior optimum; CI $(0,\SI{94}{mM})$
  is one-sided), so $|\chi|$ is essentially **unidentified**, not "measured as 0."
  **Critical caveat for the valence trend (elevate to Task B):** LaCl3's
  $|\chi|\to0$ / worst fit ($R^2\approx0.80$) is *confounded with its narrow 2–28 mM
  range and low La³⁺ rejection* (weak identifiability) — whereas CaCl2's
  $|\chi|\to0$ came from *wide*-range ($\to$\SI{108}{mM}) data. So the
  NaCl→CaCl2→LaCl3 adsorption trend rests firmly on **NaCl (wide, $|\chi|\approx50$)
  vs. CaCl2 (wide, $\to0$)**; the **LaCl3 point is consistent but range-confounded**
  (a single narrow-range dataset, no repeat — the raw-only MC4 `SLaCl3` is future
  work) → present as *suggestive, not confirmatory*. LaCl3's effective-$\chi$ drift
  is *monotonic* (rising 0.06→0.85 mM) vs. CaCl2's inverted-U — plausibly because
  the narrow range samples only the low-$c$ *rising leg* (consistent with CaCl2's
  low-$c$ behavior). Method A (scipy/Brent) vs. B (ParmEst nonlinear-constraint)
  agree on $K$ to 5 sig figs — clean cross-validation of the constraint approach
  that generalizes to multicomponent.
- **`sec:naclAR1`:** the profile "$\sim8\times$" was a *double-count* (whitening AND
  $n_{\mathrm{eff}}$ df — alternative, not complementary corrections); done
  consistently, both methods give $\sim$2.7×; the point estimate is essentially
  unchanged (the correction is to *uncertainty*, not the estimate); AR(1) is *not
  fully sufficient* (whitened ACF still structured, negative lag-1).
- **`sec:altform`:** the raw and transformed fits impose the *same* flux-balance
  equation but place the error differently, so the 2–3× shift is a genuine
  errors-in-variables consequence — **but** the near-unit-root residuals
  ($\hat\phi\to0.999$) mean the raw fit has strong *systematic misfit*, so **neither
  estimate is trustworthy yet**; it is a reduced stand-in for the full DATA~2.0 DAE.

### Cross-cutting methodological takeaways (elevate into the synthesis section, Task B)
1. **Precision $\neq$ accuracy** — a tight confidence region reflects internal
   precision, not correctness; low-information or reconstructed runs can be tightly
   identified yet biased.
2. **Identifiability requires concentration range / curvature** — the Donnan
   curvature that pins $|\chi|$ is a higher-concentration feature; narrow-range
   (MC2-NaCl, LaCl3) or low-rejection (MC3 07.09.24) experiments weakly identify or
   fail to converge. Design wide-range experiments for parameter estimation.
3. **Reconstruction-dominated $\neq$ physics** — parameters from raw$\to$processed
   reconstructed runs are limited by reconstruction uncertainty (calibration, $\Jw$
   smoothing, CP), not membrane physics alone.
4. **Working assumptions bias absolute values, not diagnostic *shapes*** — the
   lumped $K$, \ce{Cl-} stoichiometry, ambipolar $\Dm$, and the $\delta^\circ$
   treatment for $z\ne1$ shift the absolute $|\chi|/K$ but *not* the residual-vs-
   concentration / effective-$\chi$ adsorption diagnostics — so the qualitative
   conclusions (clean 1:1, adsorption-driven 2:1/3:1) are robust even while the
   absolute parameters await convention convergence.
5. **Autocorrelation is pervasive** — dense single-run time series make i.i.d.
   confidence regions systematically too narrow ($\sim$2.7× for NaCl); the same
   correction is needed for every salt.
6. **Formulation / errors-in-variables sensitivity** — transforming measurements
   into the response ($\Delta c_{\mathrm m}$) vs. regressing on raw signals
   ($\kappa_p$) materially changes the estimate; the trustworthy answer depends on
   where measurement error is placed and needs the full measurement model.
7. **Cross-validation of independent implementations** (scipy vs. ParmEst; the
   AR(1) both-ways check) is how several of these subtleties (the local minimum, the
   `measurement_error` trap, the double-count) were caught — worth stating as a
   practice.

## Guardrails
- Do **not** change analysis results/numbers, the model math (`sec:model`),
  notation, `\author`/preamble, `refs.bib`, `data/`, `prior_analysis/`, or any
  figure-generating code's numerics. This is text/structure/consistency only
  (you may add a synthesis section and an agenda table).
- Compile clean: `latexmk -pdf` → **0 undefined citations, 0 overfull hboxes**;
  **rasterize (`pdftoppm`) and eyeball** any new/changed table or the agenda; run
  `latexmk -c` after (`main.pdf` is gitignored).

## Report back
- A **section-by-section change log** (what was consolidated, the new synthesis
  section, the stale spots fixed, any reorder).
- The consolidated agenda's item count and where it lives.
- For **Task E**: which per-section nuances were already present vs. added, and
  which cross-cutting takeaways you elevated into the synthesis section.
- Confirmation of a clean compile, final page count, and that results/numbers,
  `refs.bib`, and `data/` are unchanged.
