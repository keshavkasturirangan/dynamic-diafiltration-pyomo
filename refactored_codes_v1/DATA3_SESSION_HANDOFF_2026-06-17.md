# DATA3 session handoff — 2026-06-17 (resume in a new thread)

**Repo root:** `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo`  ·  work dir: `refactored_codes_v1/`
**Artifacts root (`$B`):** `UnifiedFramework/DATA3/results/paper_artifacts/nf270/bform_study/`

## TL;DR — to continue
```bash
cd refactored_codes_v1
bash RESUME_DATA3.sh            # resumes all campaigns (idempotent) + rebuilds both decks
bash RESUME_DATA3.sh monitor    # quick status snapshot
```
Memory files (auto-loaded in a new thread) hold the science: `data3-aic-and-contours-2026-06-16`,
`data3-peclet-pooling-deck-2026-06-17`, `data3-bform-deck-2026-06-16`.

## State at handoff
- per-B-form σ×Lp stacking: **slice 23/24 done, profile 6/24** (the profile pass is the long pole, ~1–2 h left).
- constant-B 3-plane stacking (σ×Lp, B×Lp, σ×B): **done** all 4 groups.
- per-experiment deck contours (`profile_contours/`): **20/27** (paused; resumes via the campaign cmd).
- AIC, pooling (stage1/xverify/sharesigma), Peclet, predictions: **done**.

## Deliverables
- **`DATA3_FITTING_4experiments.pptx`** (12 slides) — the meeting deck: 4 representative fits
  (NaCl diluting/concentrating=constant B, CaCl₂/LaCl₃=quadratic B(c)) + Peclet rationale + AIC rationale
  + 4 supporting σ×Lp contours + summary. Contour style = labeled line contours (matches `DATA3_contour_diagnostics_2026-06-04_v7.pptx`).
- **`DATA3_STORY_deck.pptx`** — full story (Peclet → AIC → identifiability → shared-σ + stacking → per-B-form matrix);
  auto-fills as figures complete (`python3 _build_story_deck.py`).

## Key scripts (all in `refactored_codes_v1/`)
| script | what |
|---|---|
| `_run_profile_contours.py` | profile/slice contours + stacking. Cmds: `stack <grp> <form> <density> <method>`, `bformsweep <grp> <density> <method>`, `campaign deck <density> <par>`, `plot <rid> <form> <plane>`. Contour render = `_contour_panel` (labeled line contours). |
| `_run_pooling.py` | `stage1` (FIM/optimality, σ-rail), `xverify` (transfer), `sharesigma` (Stage 3 shared-σ), `lacl3refit` |
| `_make_aic_selection.py` | DATA2 per-response AIC (eq 39) — the B(c) order selector |
| `_make_peclet.py` | Peclet regime (DATA2 Fig 7/8) |
| `_make_predictions.py <sheets> <forms>` | measured-vs-predicted time series |
| `_build_fitting_deck.py` / `_build_story_deck.py` | the two decks (idempotent, figure-aware) |

## Key scientific findings
- **Peclet regime explains B(c):** diluting-NaCl Pe→0 (diffusion, B const); CaCl₂/conc-NaCl/LaCl₃ Pe≫1 (convection+diffusion, B(c)).
- **AIC (per-response, eq 39):** CaCl₂ → quadratic/linear, cubic REJECTED (Δ=8–25); NaCl flat. NOT the pooled JSON AIC (over-picks cubic).
- **Identifiability:** σ rails 9/11 single experiments (cond 1e6–1e12). **Shared-σ pooling** (per salt+regime group) is ~free
  (cost 0.99–1.15×) and resolves the rail (diluting-NaCl σ*=0.75). **Stacked objective landscapes** (= summed Fisher info)
  localize where individuals are sloppy = experiments informing each other.
- **Same NF270 membrane** across DATA3 → Lₚ consistent; within-group scatter is identifiability noise.

## Gotchas
- casadi/idas Simulator init is UNCAPPED → grinds on stiff high-σ nodes. `skip_sim_init=True` (opt-in, added to
  `solve_model` + `solve_model_B_fix`; DATA1/DATA2 guards 6/6) skips it — good for contour nodes, HURTS full fits.
- contour campaigns: warm-start continuation + parent hard-kill on stall (process-group SIGKILL); resumable via nodes.jsonl.
- session-only crons die with the thread — re-arm a monitor in the new thread if you want auto-refresh (or just re-run the deck builders).
- 4 representative experiments: diluting NaCl=MC3.07.22.24_SNaCl, conc NaCl=MC2.05.07.24_NaCl, CaCl₂=MC2.05.07.24_CaCl2, LaCl₃=MC2.05.21.24_LaCl3.

## Open / next
- finish the per-B-form **profile** pass (RESUME_DATA3.sh does it).
- backfill NaCl-concentrating constant-B σ×Lp contour (deck currently uses its linear-B(c) one as nearest).
- restyle `_stackplot` to labeled line contours too (only `_plot` was restyled this session).
- Stage 4 (joint multi-experiment fit: shared σ + B(c) shape, warm-started from the stacked minima).
