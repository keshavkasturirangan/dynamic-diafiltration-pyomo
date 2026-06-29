#!/usr/bin/env python3
"""DATA3 multi-experiment IDENTIFIABILITY / pooling study.

Staged plan (user-directed):
  stage1  per-experiment fit + full diagnostic suite (params, 3-channel objectives,
          FIM, A/D/E/modified-E optimality, AIC) -> identifiability baseline + the
          initial guesses every later stage warm-starts from.
  xverify cross-verification: fix the best-identified experiment's (B, sigma) and test
          it predicts the others in its group.                              [stage 2]
  sigma   share sigma across each salt-regime group (per-experiment Lp, B).  [stage 3]
  full    share sigma + B(c) shape (per-experiment Lp, B-scale).             [stage 4]

Grouping (salt AND regime; SAME NF270 membrane across the whole DATA3 campaign, so the
within-group Lp spread is identifiability noise, not coupon variability):
  NaCl_diluting       MC3.07.22.24_SNaCl, MC5.07.23.24_S2NaCl
  NaCl_concentrating  MC4.07.11.24_SNaCl, MC2.05.07.24_NaCl, MC5.07.23.24_NaCl, MC5.07.23.24_SNaCl
  CaCl2               MC2.05.07.24_CaCl2, MC3.07.11.24_SCaCl2, MC3.07.12.24_S2CaCl2
  LaCl3               MC2.05.21.24_LaCl3, MC4.07.11.24_SLaCl3

Output: bform_study/pooling/
Usage: python3 _run_pooling.py stage1 [form=single]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, csv, json, math, contextlib, io
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _run_bform_campaign as camp
import refactored_ucb_library as lib
lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0

STUDY = camp.OUT
OUT = STUDY / "pooling"
OUT.mkdir(parents=True, exist_ok=True)
NFE = 100

GROUPS = {
    "NaCl_diluting":      ["MC3.07.22.24_SNaCl", "MC5.07.23.24_S2NaCl"],
    "NaCl_concentrating": ["MC4.07.11.24_SNaCl", "MC2.05.07.24_NaCl", "MC5.07.23.24_NaCl", "MC5.07.23.24_SNaCl"],
    "CaCl2":              ["MC2.05.07.24_CaCl2", "MC3.07.11.24_SCaCl2", "MC3.07.12.24_S2CaCl2"],
    "LaCl3":              ["MC2.05.21.24_LaCl3", "MC4.07.11.24_SLaCl3"],
}
SHEET_GROUP = {s: g for g, ss in GROUPS.items() for s in ss}
PARAM_NAMES = ["Lp", "B", "sigma", "S0", "S"]   # constant-B parameter vector order


def load_single(sheet):
    p = STUDY / sheet / "result_single.json"
    if not p.exists():
        return None
    r = json.loads(p.read_text())
    return None if "error" in r else r


def optimality(eig):
    """A/D/E/modified-E criteria from the FIM eigenvalues (legacy doe_heatmap convention)."""
    eig = np.asarray([e for e in eig if e is not None], float)
    eig = eig[np.isfinite(eig)]
    if eig.size == 0:
        return {}
    pos = eig[eig > 0]
    return {
        "A_trace": float(np.sum(eig)),                       # total information (maximize)
        "D_logdet": float(np.sum(np.log(pos))) if pos.size else float("-inf"),  # log det FIM (maximize)
        "E_eigmin": float(np.min(eig)),                      # least-informative direction (maximize)
        "modE_cond": float(np.max(eig) / np.min(pos)) if pos.size else float("inf"),  # condition # (minimize)
    }


def fim_for(sheet, params, form="single"):
    """Reuse stored FIM if present, else compute calc_FIM (constant-B)."""
    r = load_single(sheet)
    if r and "FIM" in r and r["FIM"].get("std"):
        f = r["FIM"]
        eig = []
        if f.get("eig_min") is not None and f.get("eig_max") is not None:
            eig = [f["eig_min"], f["eig_max"]]   # only extremes stored; recompute full set below if needed
        return {"std": f.get("std"), "eig_min": f.get("eig_min"), "eig_max": f.get("eig_max"),
                "cond": f.get("cond"), "det": f.get("det"), "eig_val": None, "source": "stored"}
    # compute
    ds = camp.load(sheet)
    with contextlib.redirect_stdout(io.StringIO()):
        doe = lib.calc_FIM(ds, ds["mode"], theta=params, B_form=form, workflow_family="DATA3", nfe=NFE)
    ev = [float(x) for x in (doe.get("eig_val") or [])]
    return {"std": [float(x) for x in (doe.get("std") or [])],
            "eig_min": (min(ev) if ev else None), "eig_max": (max(ev) if ev else None),
            "cond": (max(ev) / min(ev) if ev and min(ev) > 0 else None),
            "det": float(doe["det"]) if doe.get("det") is not None else None,
            "eig_val": ev, "source": "computed"}


def aic_selected(sheet):
    p = STUDY / "model_selection" / "aic_selection.json"
    if not p.exists():
        return None, None
    d = json.loads(p.read_text()).get(sheet, {})
    return d.get("selected_label"), (d.get("forms", {}).get("single", {}) or {}).get("AIC")


def stage1(form="single"):
    rows = []
    print("=== STAGE 1: per-experiment identifiability diagnostics (constant-B baseline) ===")
    for g, sheets in GROUPS.items():
        print(f"\n-- {g} --")
        for sheet in sheets:
            r = load_single(sheet)
            if r is None:
                print(f"   {sheet:22s}  (no constant-B fit)"); continue
            p = r["parameters"]
            sg = p["sigma"]; railed = (sg <= 1e-3 or sg >= 0.999)
            fim = fim_for(sheet, p, form)
            eig = fim.get("eig_val")
            opt = optimality(eig) if eig else {
                "A_trace": None, "D_logdet": None,
                "E_eigmin": fim.get("eig_min"), "modE_cond": fim.get("cond")}
            sel, aic = aic_selected(sheet)
            # sigma std (index 2 of the param vector)
            std = fim.get("std") or []
            sg_std = std[2] if len(std) > 2 else None
            row = {"group": g, "salt": g.split("_")[0], "coupon": sheet.split(".")[0], "sheet": sheet,
                   "Lp": p["Lp"], "B": p.get("B"), "sigma": sg, "sigma_railed": railed,
                   "sigma_std": sg_std, "obj_m": r["obj_m"], "obj_cv": r["obj_cv"], "obj_cr": r["obj_cr"],
                   "WSSE3": r["WSSE3"], "cond_modE": fim.get("cond"), "eig_min_E": fim.get("eig_min"),
                   "A_trace": opt.get("A_trace"), "D_logdet": opt.get("D_logdet"),
                   "aic_selected_form": sel, "fim_source": fim.get("source")}
            rows.append(row)
            cond = fim.get("cond")
            print(f"   {sheet:22s} Lp={p['Lp']:5.2f} B={p.get('B'):6.2f} sig={sg:.2f}"
                  f"{'(RAIL)' if railed else '     '}  cond={cond:.2e}  "
                  f"Emin={fim.get('eig_min'):.1f}  selAIC={sel}  [{fim.get('source')}]")
    # write
    cols = ["group", "salt", "coupon", "sheet", "Lp", "B", "sigma", "sigma_railed", "sigma_std",
            "obj_m", "obj_cv", "obj_cr", "WSSE3", "cond_modE", "eig_min_E", "A_trace", "D_logdet",
            "aic_selected_form", "fim_source"]
    with open(OUT / "stage1_diagnostics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow(r)
    json.dump({"groups": GROUPS, "rows": rows}, open(OUT / "stage1.json", "w"), indent=2)

    # figure: condition number (identifiability) per experiment, grouped; sigma-rail marked
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    gc = {"NaCl_diluting": "#1f77b4", "NaCl_concentrating": "#4aa3df", "CaCl2": "#2ca02c", "LaCl3": "#d62728"}
    xs = list(range(len(rows)))
    conds = [r["cond_modE"] if r["cond_modE"] else np.nan for r in rows]
    cols_ = [gc[r["group"]] for r in rows]
    ax1.bar(xs, np.log10(conds), color=cols_)
    ax1.set_ylabel("log$_{10}$ condition number (modified-E)\nhigher = worse identifiability")
    ax1.set_xticks(xs); ax1.set_xticklabels([r["sheet"].split("_")[0] for r in rows], rotation=60, ha="right", fontsize=6.5)
    ax1.set_title("Per-experiment ill-conditioning (constant-B)")
    for i, r in enumerate(rows):
        if r["sigma_railed"]:
            ax1.text(i, np.log10(r["cond_modE"]) + 0.05, "σ↯", ha="center", fontsize=8, color="k")
    sgs = [r["sigma"] for r in rows]
    ax2.bar(xs, sgs, color=cols_)
    ax2.axhline(1.0, color="k", ls="--", lw=0.7); ax2.axhline(0.0, color="k", ls="--", lw=0.7)
    ax2.set_ylabel("fitted σ  (railed to 0/1 = unidentified)")
    ax2.set_xticks(xs); ax2.set_xticklabels([r["sheet"].split("_")[0] for r in rows], rotation=60, ha="right", fontsize=6.5)
    ax2.set_title("σ rails in most single experiments → the pooling target")
    from matplotlib.patches import Patch
    ax1.legend(handles=[Patch(color=c, label=g) for g, c in gc.items()], fontsize=7, loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "stage1_identifiability.png", dpi=150); plt.close(fig)

    # group summary
    print("\n=== group identifiability summary ===")
    for g in GROUPS:
        gr = [r for r in rows if r["group"] == g]
        nrail = sum(1 for r in gr if r["sigma_railed"])
        conds_g = [r["cond_modE"] for r in gr if r["cond_modE"]]
        print(f"  {g:20s}: {len(gr)} exp, σ railed in {nrail}/{len(gr)}, "
              f"cond range {min(conds_g):.1e}–{max(conds_g):.1e}")
    print("\nSTAGE1 DONE ->", OUT / "stage1_identifiability.png")


# ----------------------------------------------------------------- stage 2: cross-verify
# best-identified experiment per group (lowest cond / identified sigma) -> the anchor
ANCHOR = {"NaCl_diluting": "MC5.07.23.24_S2NaCl", "NaCl_concentrating": "MC4.07.11.24_SNaCl",
          "CaCl2": "MC2.05.07.24_CaCl2", "LaCl3": "MC2.05.21.24_LaCl3"}


def load_form(sheet, tag):
    p = STUDY / sheet / f"result_{tag}.json"
    if not p.exists():
        return None
    r = json.loads(p.read_text())
    return None if "error" in r else r


def xverify():
    """Fix each group anchor's linear B(c) (beta_0,beta_1) + sigma; refit only Lp on every
    group-mate; compare the transferred WSSE3 to that experiment's own best fit."""
    import contextlib, io
    rows = []
    print("=== STAGE 2: cross-verification (anchor linear B(c)+sigma transferred; Lp refit) ===")
    for g, sheets in GROUPS.items():
        anc = ANCHOR[g]
        ap = load_form(anc, "poly1")
        if ap is None:
            print(f"-- {g}: anchor {anc} has no linear fit -- skipped"); continue
        b0, b1, sga = ap["parameters"]["beta_0"], ap["parameters"]["beta_1"], ap["parameters"]["sigma"]
        print(f"\n-- {g}  anchor={anc}  (beta_0={b0:.3f}, beta_1={b1:.4f}, sigma={sga:.3f}) --")
        for sheet in sheets:
            own1 = load_form(sheet, "poly1")
            own_s = load_single(sheet)
            own_wsse = own1["WSSE3"] if own1 else (own_s["WSSE3"] if own_s else None)
            own_form = "linear" if own1 else "constant"
            own_lp = (own1 or own_s)["parameters"]["Lp"] if (own1 or own_s) else None
            ds = camp.load(sheet)
            seed = {"Lp": own_lp or 5.0, "beta_0": b0, "beta_1": b1, "sigma": sga}
            for kk in ("S0", "S"):
                if own_s and kk in own_s["parameters"]:
                    seed[kk] = own_s["parameters"][kk]
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    fit, _, _ = lib.solve_model_B_fix(ds, ds["mode"], theta=seed, sim_opt=False, B_form=1,
                                                      sigma_fixed=True, workflow_family="DATA3", nfe=NFE,
                                                      solver_max_cpu_time=90)
                if isinstance(fit, dict):
                    tw = sum(float(fit.get(k, 0.0)) for k in ("obj_m", "obj_cv", "obj_cr"))
                    tlp = fit["parameters"]["Lp"]
                    ratio = (tw / own_wsse) if own_wsse else None
                    status = "ok"
                else:
                    tw = tlp = ratio = None; status = "infeasible"
            except Exception as ex:
                tw = tlp = ratio = None; status = f"err:{str(ex)[:30]}"
            is_anchor = (sheet == anc)
            rows.append({"group": g, "anchor": anc, "experiment": sheet, "is_anchor": is_anchor,
                         "own_form": own_form, "own_WSSE3": own_wsse, "own_Lp": own_lp,
                         "transferred_WSSE3": tw, "refit_Lp": tlp, "ratio": ratio, "status": status})
            r_s = f"{ratio:.2f}x" if ratio else status
            lp_s = f"{own_lp:.2f}->{tlp:.2f}" if (own_lp and tlp) else "-"
            print(f"   {sheet:22s}{' (anchor)' if is_anchor else '         '}  own({own_form}) WSSE3="
                  f"{own_wsse:7.1f}  transferred={tw if tw is None else round(tw,1)}  ratio={r_s}  Lp {lp_s}")
    # write + summary
    cols = ["group", "anchor", "experiment", "is_anchor", "own_form", "own_WSSE3", "own_Lp",
            "transferred_WSSE3", "refit_Lp", "ratio", "status"]
    with open(OUT / "xverify_table.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow(r)
    json.dump(rows, open(OUT / "xverify.json", "w"), indent=2)

    # figure: transfer-degradation ratio per group-mate
    fig, ax = plt.subplots(figsize=(12, 5))
    gc = {"NaCl_diluting": "#1f77b4", "NaCl_concentrating": "#4aa3df", "CaCl2": "#2ca02c", "LaCl3": "#d62728"}
    plot = [r for r in rows if r["ratio"]]
    xs = list(range(len(plot)))
    ax.bar(xs, [r["ratio"] for r in plot], color=[gc[r["group"]] for r in plot],
           edgecolor=["k" if r["is_anchor"] else "none" for r in plot], linewidth=1.5)
    ax.axhline(1.0, color="k", ls="-", lw=0.8); ax.axhline(1.5, color="gray", ls="--", lw=0.7)
    ax.text(len(plot) - 0.5, 1.52, "1.5x (transfers OK below)", ha="right", fontsize=7)
    ax.set_ylabel("transferred WSSE3 / own-best WSSE3\n(1 = anchor correlation fits as well as own fit)")
    ax.set_xticks(xs); ax.set_xticklabels([r["experiment"].split("_")[0] for r in plot], rotation=60, ha="right", fontsize=6.5)
    ax.set_title("Stage 2 cross-verification: does the anchor's linear B(c)+σ predict its group? (black edge = anchor)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c, label=g) for g, c in gc.items()], fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / "xverify_transfer.png", dpi=150); plt.close(fig)

    print("\n=== transfer summary (ratio = transferred/own; ~1 transfers, >>1 does not) ===")
    for g in GROUPS:
        gr = [r for r in rows if r["group"] == g and r["ratio"]]
        if gr:
            rr = [r["ratio"] for r in gr if not r["is_anchor"]]
            print(f"  {g:20s}: group-mate ratios {[round(x,2) for x in rr]}  (anchor {ANCHOR[g].split('_')[0]})")
    print("\nXVERIFY DONE ->", OUT / "xverify_transfer.png")


# ----------------------------------------------------------------- stage 3: share sigma
def sharesigma(n_sigma=13):
    """Per group, profile a SHARED sigma: at each sigma fix it on every experiment and
    refit per-experiment (Lp, B constant); sum the group WSSE3.  sigma* = argmin.
    A clear interior minimum => pooling identifies sigma (rail resolved)."""
    import contextlib, io
    SIG = np.linspace(0.0, 1.0, n_sigma)
    summary = {}
    print("=== STAGE 3: shared-sigma pooling (per-exp Lp,B; warm-started over sigma) ===")
    for g, sheets in GROUPS.items():
        print(f"\n-- {g} --")
        dss = {s: camp.load(s) for s in sheets}
        seeds = {}
        indiv = {}   # individual best WSSE3 (per-exp sigma) for the baseline
        for s in sheets:
            r = load_single(s)
            seeds[s] = dict(r["parameters"]) if r else {"Lp": 5.0, "B": 5.0, "sigma": 1.0}
            indiv[s] = r["WSSE3"] if r else None
        per_sigma = {}
        for i in sorted(range(n_sigma), key=lambda i: -SIG[i]):   # march high->low sigma (warm-start)
            sg = float(SIG[i]); tot = 0.0; ok = True; rd = {}
            for s in sheets:
                seed = dict(seeds[s]); seed["sigma"] = sg
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        fit, _, _ = lib.solve_model_B_fix(dss[s], dss[s]["mode"], theta=seed, sim_opt=False,
                                                          B_form="single", workflow_family="DATA3", nfe=NFE,
                                                          solver_max_cpu_time=60, fix_vars={"sigma": sg})
                                                          # simulator init + warm-start over sigma -> good seeds, no grind
                    if isinstance(fit, dict):
                        w = sum(float(fit.get(k, 0.0)) for k in ("obj_m", "obj_cv", "obj_cr"))
                        seeds[s] = dict(fit["parameters"])
                        rd[s] = {"WSSE3": w, "Lp": fit["parameters"]["Lp"], "B": fit["parameters"].get("B")}
                        tot += w
                    else:
                        ok = False; rd[s] = None
                except Exception:
                    ok = False; rd[s] = None
            per_sigma[round(sg, 3)] = {"total": (tot if ok else None), "exp": rd}
            print(f"   sigma={sg:.3f}  group WSSE3={tot:8.1f}{'' if ok else '  (incomplete)'}")
        # sigma*
        valid = {s: v["total"] for s, v in per_sigma.items() if v["total"] is not None}
        sg_star = min(valid, key=valid.get) if valid else None
        indiv_sum = sum(v for v in indiv.values() if v is not None)
        pooled = valid.get(sg_star) if sg_star is not None else None
        summary[g] = {"sigma_star": sg_star, "pooled_WSSE3": pooled, "indiv_sum_WSSE3": indiv_sum,
                      "cost_ratio": (pooled / indiv_sum) if (pooled and indiv_sum) else None,
                      "profile": per_sigma, "indiv": indiv,
                      "exp_at_star": per_sigma.get(sg_star, {}).get("exp") if sg_star is not None else None}
        print(f"   => sigma* = {sg_star}  pooled WSSE3={pooled:.1f}  vs indiv-sum {indiv_sum:.1f}  "
              f"(cost {pooled/indiv_sum:.2f}x)" if pooled else "   => no feasible profile")
    json.dump(summary, open(OUT / "sharesigma.json", "w"), indent=2)

    # figure: group WSSE3(sigma) profile, sigma* marked, individual railed-sigma ticks
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.3))
    for ax, g in zip(axes, GROUPS):
        prof = summary[g]["profile"]
        xs = sorted(prof.keys()); ys = [prof[x]["total"] for x in xs]
        ax.plot(xs, ys, "o-", color="#333")
        ss = summary[g]["sigma_star"]
        if ss is not None:
            ax.axvline(ss, color="#2ca02c", lw=2, label=f"σ*={ss:.2f}")
        for s, p in summary[g]["indiv"].items():
            pass
        ax.set_title(f"{g}\npooled {summary[g]['cost_ratio']:.2f}x indiv" if summary[g]["cost_ratio"] else g, fontsize=9)
        ax.set_xlabel("shared σ"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
        if g == list(GROUPS)[0]:
            ax.set_ylabel("group Σ WSSE3")
    fig.suptitle("Stage 3 — shared-σ profile per group: a sharp interior minimum = pooling identifies σ "
                 "(resolves the single-experiment rail)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94]); fig.savefig(OUT / "sharesigma_profile.png", dpi=150); plt.close(fig)
    print("\nSHARESIGMA DONE ->", OUT / "sharesigma_profile.png")


# ----------------------------------------------------------------- LaCl3 stuck-optimum re-fit
def lacl3refit():
    """Re-fit the 2 LaCl3 sheets with multistart seeded from the better basin the
    cross-verification exposed (sigma~0.29, Lp~4.7).  Update result_<form>.json only if
    WSSE3 improves.  Corrects the stuck-optimum LaCl3 fits (affects AIC + pooling)."""
    import contextlib, io
    sheets = ["MC2.05.21.24_LaCl3", "MC4.07.11.24_SLaCl3"]
    forms = [("single", "single"), (1, "poly1"), (2, "poly2"), (3, "poly3"), ("sat", "sat")]
    NB = {"single": 1, "poly1": 2, "poly2": 3, "poly3": 4, "sat": 2}
    print("=== LaCl3 RE-FIT (multistart from the better basin sigma~0.29, Lp~4.7) ===")
    for sheet in sheets:
        ds = camp.load(sheet); mode = ds["mode"]
        print(f"\n-- {sheet} --")
        for bform, tag in forms:
            cur = load_form(sheet, tag)
            cur_w = cur["WSSE3"] if cur else float("inf")
            # seed set: current fit + better basin (low sigma) variants
            seeds = []
            if cur:
                seeds.append(dict(cur["parameters"]))
            for sg in (0.29, 0.5):
                base = {"Lp": 4.7, "sigma": sg, "B": 0.2, "beta_0": 0.12, "beta_1": -0.003,
                        "beta_2": 0.0, "beta_3": 0.0, "B_inf": 0.3, "c_star": 8.0, "S0": 0, "S": 0.35}
                seeds.append(base)
            best, best_w = None, cur_w
            for seed in seeds:
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        fit, _, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=False, B_form=bform,
                                                    workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=120,
                                                    LOUD=False)   # simulator init: good seed -> converges (skip_sim_init hurts full fits)
                    if isinstance(fit, dict):
                        w = sum(float(fit.get(k, 0.0) or 0.0) for k in ("obj_m", "obj_cv", "obj_cr"))
                        if w < best_w:
                            best, best_w = fit, w
                except Exception:
                    pass
            if best is not None and best_w < cur_w * 0.995:
                # update the result json (preserve schema fields used downstream)
                res = {"parameters": best["parameters"],
                       "WSSE3": best_w, "Obj": float(best.get("Obj", 0.0)),
                       "obj_m": float(best["obj_m"]), "obj_cv": float(best["obj_cv"]),
                       "obj_cr": float(best["obj_cr"]), "n_Bparams": NB[tag],
                       "n_data": (cur or {}).get("n_data"), "salt": "LaCl3", "form": str(bform),
                       "refit": "lacl3_basin"}
                (STUDY / sheet / f"result_{tag}.json").write_text(json.dumps(res, indent=2))
                print(f"   {tag:7s}: {cur_w:.1f} -> {best_w:.1f}  IMPROVED, updated")
            else:
                print(f"   {tag:7s}: {cur_w:.1f} (no improvement; kept)")
    print("\nLaCl3 REFIT DONE. Re-run: python3 _make_aic_selection.py all  (to refresh AIC with the better fits)")


def main():
    if len(sys.argv) < 2:
        print(__doc__); return
    if sys.argv[1] == "stage1":
        stage1(sys.argv[2] if len(sys.argv) > 2 else "single"); return
    if sys.argv[1] == "xverify":
        xverify(); return
    if sys.argv[1] == "sharesigma":
        sharesigma(int(sys.argv[2]) if len(sys.argv) > 2 else 13); return
    if sys.argv[1] == "lacl3refit":
        lacl3refit(); return
    print(__doc__)


if __name__ == "__main__":
    main()
