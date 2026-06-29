# DATA3 NF270 — "Choosing the B(c) correlation" — session handoff

**Date:** 2026-06-16 · **Branch:** `refactored_unified_codebase` · **Repo root:** `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo`

This document captures everything from the working session so we can resume without losing context.

---

## 0. Goal

Determine the best equation correlating solute permeability **B** with interfacial concentration `c_in` (or ionic strength `I`) for the DATA3 NF270 **single-salt** campaign, so we can build good contours/fits (mass, permeate conc, retentate conc) and present it in a slide deck. Earlier finding (from the June 11 meeting): on most single-salt sheets, **B and σ were unidentifiable or railed to bounds** under a constant-B model.

**Salts / valence:** NaCl (1:1), CaCl₂ (2:1), LaCl₃ (3:1). **11 single-salt sheets.**
**Headline sheets (GOOD/OKAY/POOR):** `MC3.07.22.24_SNaCl` (NaCl), `MC2.05.07.24_NaCl`, `MC2.05.07.24_CaCl2`, `MC2.05.21.24_LaCl3`.

---

## 1. Library changes (all in `refactored_codes_v1/refactored_ucb_library.py`)

**ALL changes are DATA3-gated. DATA1/DATA2 stay byte-identical — verified 6/6 guards (`tests/regression_data1_data2_guards.py`) + DATA2 regression (`tests/test_data2_regression.py`, 2 passed/2 skipped) after every edit.**

1. **Two new B-forms** (mechanistic partition; `Js = Jw·B·(c_in−c_H)`, same Jw-scaled convention as the polynomial forms):
   - `B_form='sat'` — saturating exp: `B = B_inf·(1 − exp(−cc/c_star))`
   - `B_form='donnan'` — Donnan–dielectric: `B·2cc = P0·(−X + sqrt(X² + 4k²cc²))`
   - where `cc = m.cIn[n,t]` (or `k_I·cIn` if the ionic-strength toggle is on).
   - Wired everywhere `'convection'` is special-cased: `model_construct_inter` Var/Param decls (both `sim_opt` branches), `B_form_rule3`, `save_model` (param save + LOUD), `solve_model_B_fix` guard, `_multistart_theta_defaults`, `_multistart_theta_specs`, plus `param_in` seeds.
   - **CRITICAL fix:** `m.B[n,t]` for sat/donnan is **unbounded** (`Var(m.n_vial, m.tau)`) — a `bounds=(1e-6,50)` box made the DAE infeasible on transient IPOPT steps. Positivity comes from the defining constraint.

2. **DATA3-gated concentration-scaled β bounds** (in the polynomial `else` branch of the `sim_opt=False` decl). The legacy `β∈(−20,20)` let `β₂·c²` reach ~10⁵ at `c_in~95 mM`, making the **quadratic/cubic DAE fit infeasible** (that's why `bform_quad_test/` was always empty). Now `|β_k| ≤ 30/c_max^k` for DATA3 (DATA1/DATA2 keep ±20). **This is what made quad/cubic estimable for the first time.**

3. **`NF270_DONNAN_FIX_K`** module toggle — when set, `m.k_dd` becomes a fixed `Param` → donnan becomes a well-conditioned 2-param (P0,X) fit. Used by the donnan recovery.

4. **`fix_vars` kwarg on `solve_model_B_fix`** — when given a dict, fixes exactly those Vars at grid values and leaves everything else (incl. B(c) coefficients) **free to solve**. This is the inverse of the default (fix-B) behaviour and is the hook for **profile contours** (fix Lₚ,σ → solve B(c) per node). **Built but not yet run.**

---

## 2. Scripts created (in `refactored_codes_v1/`)

| script | what it does |
|---|---|
| `_run_bform_campaign.py` | Fits every B-form on every single-salt sheet (parallel `Pool` over sheets, thread-pinned). Multistart(N=10)+FIM on headline sheets; single-shot elsewhere. donnan/quad/cubic single-shot. Writes `bform_study/<sheet>/result_<form>.json` + `apparent_B.csv`. Resumable. Run: `python3 _run_bform_campaign.py all all 8` |
| `_run_bform_contours.py` | B-LEVEL contours (σ×Lₚ, B×Lₚ, σ×B) via `B_form='single'`, forward-sim WSSE, figure_s5 style + CSV. Cross-sheet parallel. Run: `python3 _run_bform_contours.py 30 6` |
| `_make_model_selection.py` | AICc table + analytic fitted **partition** B(c) curves (in-window + extrapolated) + B-vs-c/B-vs-I. **NOTE: the AICc averaging is biased by uneven n (forms converge on different sheets) — use per-sheet WSSE instead (see §4).** |
| `_make_pervial_overlay.py` | **The corrected method.** `B_form='pervial'` (free B per vial) → empirical apparent B = Jₛ/(c_in−c_H) vs c_in; overlays each correlation's B(c_in). Headline sheets. Output: `model_selection/pervial_overlay.png` |
| `_make_predictions.py` | Prediction-vs-measured (mass / perm / ret conc) per (sheet,form) via `run_data3_time_series_plots`. Use **NFE≤120** (NFE=150 hits a `'no attribute B'` quirk for numeric forms). |
| `_recover_donnan.py` + `_donnan_fit_worker.py` | Bounded fix-k donnan recovery — each (sheet,k) in an isolated subprocess with a hard timeout (kills the process group incl. orphaned ipopt). k∈{0.6,0.5,0.45,0.7,0.35}. Recovered LaCl₃ (k=0.6) + CaCl₂ SCaCl2 (k=0.5). |
| `_build_bform_deck.py` | Builds the deck from the template: swaps figure placeholders, rewrites the slide-10 table + §05 θ-boxes, QA-renders PDF via soffice. **Does NOT insert a β slide** (the template already has the user's). |
| `_make_aic_selection.py` | **TO BE WRITTEN** — DATA2-Table-3-style AIC ranking on DATA3 (see §6). |

---

## 3. Artifacts (`UnifiedFramework/DATA3/results/paper_artifacts/nf270/bform_study/`)

- `<sheet>/result_<form>.json` — fitted params + `WSSE3` (=obj_m+obj_cv+obj_cr, the deck metric, **NOT** the full `Obj` which carries the truncate term) + AICc + `n_Bparams` + `n_data` + FIM (headline). Forms: `single, poly1, poly2, poly3, sat, donnan`.
- `<sheet>/apparent_B.csv` — **flat / misleading** (apparent B from the *constant* fit is pinned constant). Superseded by the per-vial overlay. Do not use for B(c) trends.
- `contours/<sheet>/{objcontour,paper,contourdata}-*` — constant-B level contours, **9/11 sheets done** (MC3 SNaCl + MC5 S2NaCl didn't finish — slow integrator nodes; not needed).
- `model_selection/{fitted_Bc_curves.png, Bc_vs_ionic_strength.png, pervial_overlay.png, aicc_table.csv}`
- `predictions/<sheet>/<form>/{mass,concentration}-*.png` — single + sat for the 3 headline sheets.
- Legacy per-vial trend (rolling 3-vial window): `band_value_test/per_vial_trend-*.png` + `per_vial_trend_summary.json`.

---

## 4. Key scientific findings (per-sheet WSSE3 — apples-to-apples)

| salt (sheet) | const | linear | quad | cubic | sat | donnan |
|---|---|---|---|---|---|---|
| NaCl (MC3 SNaCl) | 42.2 | 40.0 | ERR | 30.7 | 47.6 | ERR |
| NaCl (MC2) | 114.5 | 108.4 | 102.1 | 97.5 | 112.1 | ERR |
| **CaCl₂ (MC2)** | 219.8 | 168.0 | 147.0 | 145.9 | 146.3 | (k-fix) |
| LaCl₃ (MC2) | 1244.1 | 1213.8 | 1191.8 | 1188.8 | 1240.5 | 1240.6 |

- **CaCl₂ is the clear B(c) win:** constant 220 → quad/cubic/sat ~146 (**~33% cut**). The **per-vial empirical B rises** ~2→~13 µm/s with c; the **saturating curve tracks it**, linear overshoots, quadratic turns over. (`pervial_overlay.png`, middle panel — the money figure.)
- **NaCl:** per-vial B is **flat ~12–15 µm/s** (lowest-c vial rails to the 50 bound → unidentified). B ≈ constant is adequate.
- **LaCl₃:** all forms ~1190–1244 (**model-limited** — B(c) gives only ~4%; points past B(c) alone, e.g. σ(c)/ion-pairing).
- **Campaign:** 50/66 fits. Failures: 9 donnan (fragile 3-param), 3 poly2, 3 sat, 1 poly1. donnan converged on 4 sheets (NaCl MC5 + MC5 S2NaCl; CaCl₂ MC3 SCaCl2; LaCl₃ MC2 via fix-k=0.6).

---

## 5. The deck

**File:** `refactored_codes_v1/DATA3_B_equation_for_contours_REAL.pptx` (+ `.pdf` QA). Built from the user's **updated 24-slide** template `/Users/kkasturi/Documents/Claude/Projects/Diafiltration/DATA3_B_equation_for_contours_2026-06-14.pptx` (the user added their own β slides 18 "WHY NOT A MAP PER β" + 19 "THE RIGHT β DIAGNOSTICS").

**Figure map (slide → figure):** 4→σ×Lₚ contour (MC2 NaCl); 9→`pervial_overlay`; 11→`fitted_Bc_curves`; 12→`pervial_overlay`; 13→`Bc_vs_ionic_strength`; 15→`fitted_Bc_curves`; 16→σ×B contour; 18→`fitted_Bc_curves`; 19→`fitted_Bc_curves`+σ×Lₚ; 21→NaCl **constant-B** prediction; 22→CaCl₂+LaCl₃ **constant-B** predictions.

**Deck corrections already made:**
- **Fixed a real error:** slide 10 previously claimed "LaCl₃ saturating halves WSSE (1240 vs 2618)" — that compared sat on ONE sheet to constant averaged over TWO (sat failed on the 2nd LaCl₃ sheet). On the same sheet sat≈constant. Replaced with the honest per-sheet finding (CaCl₂ ~33% cut; NaCl fine; LaCl₃ form-limited).
- §05 θ-boxes (slides 21/22) → real constant-B campaign fits + honest notes; figures switched to constant-B to match.
- slides 9/12 → per-vial empirical overlay; slide-10 demo subtitle replaced.

**Deck open items:** slide 11 still has DATA2-context numbers (β₁=0.248) — legitimate for the reconciliation but verify; slides 18/19 (user's β slides) need figure direction; `pervial_overlay.png` could be polished (CaCl₂ linear shoots off-panel — needs a y-cap).

---

## 6. IN PROGRESS — the DATA2-Taylor reframe (resume HERE)

**User's directive:** *"If cubic is good, stick to using a Taylor series expansion based on the derivation in DATA2_main.pdf. Don't simply choose cubic for the sake of it — have a good logical flow."*

**DATA2 paper** (`UnifiedFramework/DATA3/published_works/DATA2_main.pdf`) — the derivation:
- **Eq 21:** local solute flux `j_s = −D_m·dc/dz + J_w·c(z)` (diffusion + convection).
- **Eq 25 (EXACT):** `J_s = J_w·[c_in,f·H_f·e^Pe − H_p·c_in,p]/(e^Pe − 1)`, Pe = J_w·L/D_m. (The exponential / `convection` form.)
- **Eq 28:** `B = D_m·H/L` (constant) = the **diffusion-dominated limit** (Pe→0).
- **Eq 29–31:** **Taylor expansion when convection & diffusion both matter (Pe~1)** → `B(J_w,c_in,f) ≈ J_w·Σ_{i∈I} βᵢ·c_in,f^i`. Orders tried: `I ∈ {0},{0,−0.5},{0,0.5},{0,1},{0,1,2},{0,1,2,3}`.
- Paper attributes the c-dependence to **"electrostatic interactions [that] alter solute partitioning"** = Donnan/dielectric — i.e. **the saturating partition is the mechanism; the polynomial is its Taylor expansion.** (They even tried fractional `c^±0.5` — the `c^(1/z)` valence knee.)
- **Table 3 (p.11, KCl, concentrating):** linear `I={0,1}` **rank 1**, quadratic `I={0,1,2}` **rank 2**, **cubic `I={0,1,2,3}` rank 6 despite the LOWEST raw WLS (83.28)** — AIC rejects it for complexity. (The Δ values imply they used AICc with small effective n ≈ #concentration points ~12.)

**The logical flow to implement (this directly answers the user):**
1. The DATA2 transport model **derives** `B = J_w·Σβᵢc_in^i` (Taylor of the exact eq-25 solution, valid Pe~1) — physically grounded, not arbitrary.
2. **Select the order by AIC** (DATA2's own Table-3 method) — which **rejects cubic** (overfit) and selects linear/quadratic. The effective n is the number of **concentration points (~vials, ~10–12)**, NOT the ~2500 time points — that's why cubic is penalized.
3. The order **tracks valence** (partition curvature): **NaCl→linear** (B flat, matches DATA2 KCl which is also 1:1), **CaCl₂→quadratic** (B rises; const→quad cuts WSSE 33%, quad→cubic negligible ~0.7%), LaCl₃ model-limited.
4. **WHY the order tracks valence:** the βᵢ are the local Taylor coefficients of the underlying saturating Donnan/dielectric partition; higher valence → sharper knee (`B ∝ c^(1/z)`) → more curvature → higher order (and fractional `c^0.5` terms, which DATA2 anticipated).
5. **Use the polynomial in-window** (interpolation), where the Taylor expansion is valid. Don't extrapolate it (where it diverges to B<0); the saturating partition is the bounded form **if** transfer outside the window is ever needed.

**NEXT STEPS (do these to resume):**
- [ ] Write `_make_aic_selection.py` — compute the DATA2-Table-3-style **AICc ranking** per headline sheet over orders {0, 1, 2, 3} (and optionally fractional {0,0.5}), with **n = #concentration points (~vials)**, Akaike weights, ranks. Confirm: cubic rejected, linear/quad win, order tracks valence. Also print the raw WSSE3-by-order (showing diminishing returns past quadratic). Build a small table/figure for the deck.
- [ ] **Reframe deck §03** to this flow: lead with the DATA2-derived Taylor (eq 25→31) as the fitting form; AIC order selection (linear/quad, **cubic rejected**); order tracks valence; saturating partition = the *mechanism* (and extrapolation guard), not a competing fit. Update slide 10 title/table/takeaway accordingly. Possibly retitle slide 10 from "Inside the window: Taylor. Across it: saturating" to something like "Use the DATA2 Taylor in-window; the order follows valence."
- [ ] **Correct the earlier "cubic best in-window" framing** — it used raw WSSE; by DATA2's AIC it's over-parameterized. The honest message: parsimony selects linear (NaCl) / quadratic (CaCl₂).

**Also still open (lower priority):**
- [ ] Profile contours (fix Lₚ,σ → solve B(c) per node) via the new `fix_vars` hook — only if the (Lₚ,σ) identifiability landscape with B solved is still wanted (the per-vial method already answered "B vs c"). Thousands of fits; coarse grid + headline sheets.
- [ ] Polish `pervial_overlay.png` (y-cap so CaCl₂ linear doesn't shoot off-panel).
- [ ] Decide slides 18/19 figures with the user.

---

## 7. Gotchas learned this session

- **`apparent_B.csv` from the constant fit is flat** — the constant model pins B. Use **`B_form='pervial'`** (free B per vial) for the empirical B(c_in). (This is Xinhong Liu's method, `legacy/scripts/run_pre_B_dependence.py`.)
- **The casadi Simulator (idas) init is uncapped** — `solver_max_cpu_time` only caps IPOPT. Bad parameter regions (donnan, some contour nodes) can grind 10–18 min. Mitigate with process-level timeouts (see `_recover_donnan.py`).
- **NFE=150 forward-sim of numeric B-forms** raised `'ConcreteModel' object has no attribute 'B'`; **NFE≤120 works.**
- **macOS:** `timeout` cmd missing; `cat -A` missing (use python `repr`); Desktop is TCC-blocked from the shell (use the repo copy of the deck).
- **Use `WSSE3` (=obj_m+obj_cv+obj_cr)**, not the full `Obj` (truncate-dominated), as the comparison metric.
- IPOPT lives at `~/.idaes/bin/ipopt`; pandoc + LibreOffice (`/opt/homebrew/bin/soffice`) present for deck native-equations + QA.

---

## 8. One-liners to resume

```bash
cd /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1
python3 tests/regression_data1_data2_guards.py        # confirm DATA1/DATA2 still 6/6
python3 _make_pervial_overlay.py                       # empirical B(c_in) + correlations
python3 _build_bform_deck.py                           # rebuild deck + QA PDF
# NEXT: write & run _make_aic_selection.py, then reframe deck §03 (see §6)
```
