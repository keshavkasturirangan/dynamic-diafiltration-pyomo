# DATA3 B(c) order selection — DATA2 per-response AIC (verified)

**Method.** Order chosen by the DATA2 eq-39 per-RESPONSE AIC `AIC = sum_r n_r*log(obj_r) + 2k`, k=n_Bparams+6 (=p+R, +6 cancels in dAIC), reconstructed from the library's own per-channel objects (obj_m/obj_cv/obj_cr = mean scaled squared residual) and per-channel counts. NOT raw WSSE3, and NOT the pooled AIC stored in result_*.json / aicc_table.csv (that pooled AIC ranks cubic #1 on CaCl2 — the per-channel eq-39 object is what inverts that).

**Headline (CaCl2 = the clean case).** AIC selects LOW ORDER and rejects cubic on all three CaCl2 sheets: quadratic (MC2.05.07.24 dAIC=8.2, MC3.07.12 dAIC=13.2), linear (MC3.07.11 dAIC=25.2). The cleanest evidence is the permeate channel (direct B observable): obj_cv collapses 55.9->9.7 with order — genuine low-order B(c).

**NaCl complication (disclosed, not hidden).** AIC ranks cubic #1 on 3 of 6 NaCl sheets, but this is a candidate-set + channel artifact, not evidence for cubic B(c): on the two dAIC>1000 sheets (MC3.07.22 diluting, MC5.07.23_S2 diluting) the QUADRATIC returned 'solver returned no feasible fit', so cubic was compared only against const/linear/sat; and the fitted cubic B(c) is non-monotonic in-window (unphysical for a saturating Donnan/dielectric partition). MC4.07.11 ranks saturating #1 (constant is a dAIC<2 tie-break).

**LaCl3.** Model-limited: AIC weakly selects quadratic/linear; all forms within a narrow band (cubic dAIC~2.5-2.9).

**Caveat (count-weighting).** The permeate channel — the only direct B observable — has n_cv=10 vs n_m,n_cr in the hundreds, so the count-weighted AIC is biased toward lower order. This is faithful to DATA2 (same weighting), but it means: (a) cubic-rejection is partly a count effect; (b) the cleanest B(c) evidence is the permeate-channel MAGNITUDE, not the scalar AIC. The mass-domination caveat applies symmetrically — the CaCl2 linear/quad picks also improve mass while worsening retentate.

**Bottom line for the deck.** Use the DATA2-derived Taylor B=Jw*sum beta_i c^i as the form; choose the order by the per-response AIC; prefer the simplest model within dAIC<2. Result: low-order (quadratic/linear) for CaCl2; constant/low-order for NaCl; cubic is flagged-and-excluded on physical (non-monotonic B(c)) and convergence (no feasible quadratic) grounds — NOT the claim 'AIC never selects cubic'.
