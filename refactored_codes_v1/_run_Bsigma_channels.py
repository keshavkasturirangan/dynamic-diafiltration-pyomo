#!/usr/bin/env python3
"""Task 5 — 7-channel σ×B WSSE identifiability contours at FIXED Lp (DATA3).

Extends the existing per-response WSSE slice contour (_run_Bsigma_fixedLp.py:
mass / ICP-permeate / ICP-retentate) to the FOUR new conductivity response
channels, each as its own contour panel.  SEVEN channels (ΔP excluded — it is a
control input, not a fitted response):

  1 mass                          fit.obj_m            (vial mass)            [g²]
  2 ICP retentate concentration   model cF(end) vs the end-of-run retentate
                                   grab samples (Final Tube + Retentate ICP)  [mM²]
  3 ICP permeate concentration     fit.obj_cv           (cV at vial close)    [mM²]
  4 conductivity-probe retentate   σ-space: Shedlovsky(model cF) vs the inline
                                   retentate probe trace cF_exp_conductivity  [µS²]
  5 conductivity-probe permeate    σ-space: Shedlovsky(model cH) vs the 5-Hz
                                   permeate probe trace cV_perm_cond           [µS²]
  6 cond-derived retentate conc    fit.obj_cr  (the model's continuous cF fit;
                                   for DATA3 the retentate is conductivity-
                                   derived in mM — the familiar "retentate")  [mM²]
  7 cond-derived permeate conc     c-space: model cH vs invert(cV_perm_cond)   [mM²]

METHOD = SLICE (square forward solve at every node — NOT profile likelihood).
Pin Lp = Lp* (well-identified), set (σ, B) at each grid node, forward-simulate the
transport DAE once (solve_model(sim_opt=True) → zero free parameters → square),
and EVALUATE all seven per-channel WSSE from that one solve.  No per-node
optimization.  The conductivity↔concentration maps reuse the already-coded
library converters (Shedlovsky single-salt): lib.concentration_to_conductivity /
lib.conductivity_to_concentration.

The grid (axis_sigma, axis_B, Lp_fixed) is reused from the existing
Bsigma_fixedLp_nofloor slice so channels 1/3/6 reproduce the familiar contour
and the four new channels overlay on the same axes.  Weighting matches the
canonical "nofloor" variant (cF 2 %, no absolute floor).

Output (bform_study/Bsigma_channels/<rid>/single/):
  nodes.jsonl  grid.json  contourdata-7ch.csv  contour7.png (v7 titled, 7 panels)

Usage:
  python3 _run_Bsigma_channels.py exp <rid> [density]
  python3 _run_Bsigma_channels.py all [density]      # CaCl2 + LaCl3 + NaCl headline
"""
import os, sys, json, time, io, contextlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _run_Bsigma_fixedLp as bs      # B-scaling, anchors, FORMARG, MAG
import _run_profile_contours as pc    # _study(), config

NFE = 120
CPU = 30                              # IPOPT cpu cap per square node (s)
REL_COND = 0.03                       # relative weight for conductivity channels
REL_ICP = 0.02                        # relative weight for the retentate ICP anchor
REL_PERM = 0.03                       # relative weight for the permeate channels
# Small ABSOLUTE floors on the measurement-error weight, so a near-zero startup
# measurement (e.g. vial-1 permeate ICP ≈ 0.08 mM) cannot dominate the WSSE via
# its tiny relative denominator. Same convention as NF270_PERMEATE_PROBE_FLOOR_MM.
FLOOR_MM = 0.5                        # mM, for concentration-space channels
FLOOR_US = 50.0                       # µS/cm, for conductivity-space channels
FORM = "single"

# Human labels + run_ids (the three the advisor asked for)
REPS = {
    "concentrating CaCl₂": "MC2.05.07.24_CaCl2",
    "concentrating LaCl₃": "MC2.05.21.24_LaCl3",
    "concentrating NaCl":       "MC2.05.07.24_NaCl",
}

# panel order, key, pretty title, units tag
CHANNELS = [
    ("mass",        "mass",                       "g$^2$"),
    ("icp_reten",   "ICP retentate conc.",        "mM$^2$"),
    ("icp_perm",    "ICP permeate conc.",         "mM$^2$"),
    ("cond_reten",  "conductivity-probe retentate", "µS$^2$"),
    ("cond_perm",   "conductivity-probe permeate",  "µS$^2$"),
    ("cderiv_reten","cond.-derived retentate conc.", "mM$^2$"),
    ("cderiv_perm", "cond.-derived permeate conc.",  "mM$^2$"),
]
CH_KEYS = [c[0] for c in CHANNELS]


def out_dir(rid):
    return pc._study() / "Bsigma_channels" / rid / FORM


def _set_weighting(lib):
    """Canonical nofloor variant — byte-consistent with Bsigma_fixedLp_nofloor."""
    lib.NF270_CF_RESIDUAL_SCALE_FRACTION = 0.02
    lib.NF270_CF_RESIDUAL_FLOOR_MM = None


def _existing_grid(rid):
    """Reuse the axes (and Lp_fixed / Bfeed_fit) of the existing nofloor slice."""
    g = (pc._study() / "Bsigma_fixedLp_nofloor" / rid / FORM / "grid.json")
    return json.loads(g.read_text()) if g.exists() else None


# ----------------------------------------------------------- precompute targets
def _precompute(ds, salt):
    """Per-vial measured probe traces + the (σ,B)-independent permeate
    conductivity-derived concentration (inverted ONCE via the library converter).
    Retentate c-space uses the model's own continuous fit (channel 6 = obj_cr)."""
    import numpy as np
    import refactored_ucb_library as lib
    n_v0 = int(ds["data_config"].get("n_v0", 1))
    pv = {}
    for i, row in enumerate(ds["data_raw"]):
        if (i + 1) < n_v0:
            continue
        t = np.asarray(row.get("time", []), dtype=float).reshape(-1)
        ret_sig = np.asarray(row.get("cF_exp_conductivity", []), dtype=float).reshape(-1)
        perm_sig = np.asarray(row.get("cV_perm_cond", []), dtype=float).reshape(-1)
        perm_deriv = (lib.conductivity_to_concentration(perm_sig, salt, method="shedlovsky")
                      if perm_sig.size else np.array([]))
        cv = np.asarray(row.get("cV_avg", []), dtype=float).reshape(-1)
        fin = np.where(np.isfinite(cv))[0]
        cv_icp = float(cv[fin[-1]]) if fin.size else float("nan")   # ICP at vial close
        pv[i] = {"t": t, "ret_sig": ret_sig, "perm_sig": perm_sig,
                 "perm_deriv": np.asarray(perm_deriv, dtype=float).reshape(-1),
                 "cv_icp": cv_icp}
    cfg = ds["data_config"]
    anchors = [v for v in (cfg.get("cF_final_meas"), cfg.get("cF_retentate_icp_mM"))
               if v is not None and np.isfinite(float(v)) and float(v) > 0]
    return {"pv": pv, "end_anchors": [float(a) for a in anchors], "n_v0": n_v0}


def _accum(pred, meas, rel, floor=0.0):
    """Σ[(pred-meas)/max(rel*|meas|, floor)]² and N over finite, positive meas."""
    import numpy as np
    pred = np.asarray(pred, float); meas = np.asarray(meas, float)
    m = np.isfinite(pred) & np.isfinite(meas) & (meas > 0)
    if not m.any():
        return 0.0, 0
    scale = np.maximum(rel * np.abs(meas[m]), floor)
    r = (pred[m] - meas[m]) / scale
    return float(np.sum(r * r)), int(m.sum())


def _seven_channels(fit, sim, precomp, salt):
    """Compute all 7 per-channel WSSE from ONE forward sim + precomputed traces.

    Mass and the conductivity-derived retentate use the model's own objective
    (mass weight is absolute; cF weight already floor-protected via the nofloor
    variant).  The permeate ICP and the conductivity channels are formed here at
    the measurement points with the same relative measurement-error weighting,
    plus a small absolute floor so a near-zero startup measurement cannot
    dominate the WSSE."""
    import numpy as np
    import refactored_ucb_library as lib
    ch = {"mass": float(fit.get("obj_m")),
          "cderiv_reten": float(fit.get("obj_cr"))}

    # 2 — ICP retentate end anchors vs model cF at the final vial close
    last = max(sim.keys())
    cF_end = float(np.asarray(sim[last]["cF"], float).reshape(-1)[-1])
    s2, n2 = _accum([cF_end] * len(precomp["end_anchors"]),
                    precomp["end_anchors"], REL_ICP, FLOOR_MM)
    ch["icp_reten"] = s2 / n2 if n2 else float("nan")

    # 3 — ICP permeate: model cV at vial close vs the per-vial ICP scalar (floored)
    # 4,5,7 — accumulate the conductivity channels over the fitted vials
    s3 = s4 = s5 = s7 = 0.0; n3 = n4 = n5 = n7 = 0
    for i, pv in precomp["pv"].items():
        if i not in sim:
            continue
        t_sim = np.asarray(sim[i]["time"], float).reshape(-1)
        if t_sim.size == 0:
            continue
        cF_sim = np.asarray(sim[i]["cF"], float).reshape(-1)
        cH_sim = np.asarray(sim[i]["cH"], float).reshape(-1)
        cV_close = float(np.asarray(sim[i]["cV"], float).reshape(-1)[-1])
        r = _accum([cV_close], [pv["cv_icp"]], REL_PERM, FLOOR_MM); s3 += r[0]; n3 += r[1]
        t = pv["t"]
        if pv["ret_sig"].size == t.size and t.size:
            cF_at = np.interp(t, t_sim, cF_sim)
            sig_pred = lib.concentration_to_conductivity(cF_at, salt, method="shedlovsky")
            r = _accum(sig_pred, pv["ret_sig"], REL_COND, FLOOR_US); s4 += r[0]; n4 += r[1]
        if pv["perm_sig"].size == t.size and t.size:
            cH_at = np.interp(t, t_sim, cH_sim)
            sig_pred = lib.concentration_to_conductivity(cH_at, salt, method="shedlovsky")
            r = _accum(sig_pred, pv["perm_sig"], REL_COND, FLOOR_US); s5 += r[0]; n5 += r[1]
            # 7 — c-space permeate: model cH vs inverted permeate probe
            r = _accum(cH_at, pv["perm_deriv"], REL_PERM, FLOOR_MM); s7 += r[0]; n7 += r[1]
    ch["icp_perm"] = s3 / n3 if n3 else float("nan")
    ch["cond_reten"] = s4 / n4 if n4 else float("nan")
    ch["cond_perm"] = s5 / n5 if n5 else float("nan")
    ch["cderiv_perm"] = s7 / n7 if n7 else float("nan")
    return ch


# --------------------------------------------------------------------- one node
def _node(ds, mode, anchor, sigma_val, B_val, precomp, salt):
    import numpy as np
    import refactored_ucb_library as lib
    seed = dict(anchor); seed["sigma"] = float(sigma_val); seed["B"] = float(B_val)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fit, sim, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=True,
                                          B_form=bs.FORMARG[FORM], workflow_family="DATA3",
                                          nfe=NFE, solver_max_cpu_time=CPU)
        if not isinstance(fit, dict) or sim is None:
            return None
        return _seven_channels(fit, sim, precomp, salt)
    except Exception:
        return None


# -------------------------------------------------------------------- per-exp
def run_exp(rid, density=11):
    import numpy as np
    import refactored_ucb_library as lib
    import _run_bform_campaign as camp
    _set_weighting(lib)
    anchor = bs._anchor(rid, FORM)
    if anchor is None:
        print(f"[exp] {rid}: no converged fit — skip", flush=True); return None
    ds = camp.load(rid); mode = ds["mode"]
    salt = str(ds["data_config"].get("namec", "")).strip()

    g = _existing_grid(rid)
    if g is not None:
        sig_axis = list(g["axis_sigma"]); B_axis = list(g["axis_B"])
        Lp_fixed = float(g["Lp_fixed"]); Bfit = float(g.get("Bfeed_fit", anchor.get("B", 1.0)))
    else:                                                  # fresh grid
        Bfit = float(anchor.get("B", 1.0)); Lp_fixed = float(anchor["Lp"])
        sig_axis = list(np.linspace(0.0, 1.0, density))
        B_axis = list(np.linspace(max(1e-4, 0.3 * Bfit), 1.7 * Bfit, density))
    anchor = dict(anchor); anchor["Lp"] = Lp_fixed         # pin Lp at optimum

    precomp = _precompute(ds, salt)
    od = out_dir(rid); od.mkdir(parents=True, exist_ok=True)
    grid = {"rid": rid, "form": FORM, "Lp_fixed": Lp_fixed, "Bfeed_fit": Bfit,
            "sigma_fit": float(anchor.get("sigma", 1.0)), "salt": salt,
            "axis_sigma": [float(s) for s in sig_axis], "axis_B": [float(b) for b in B_axis],
            "channels": CH_KEYS}
    (od / "grid.json").write_text(json.dumps(grid))
    jsonl = od / "nodes.jsonl"
    done = set()
    if jsonl.exists():
        for ln in jsonl.read_text().splitlines():
            try:
                r = json.loads(ln); done.add((r["i"], r["j"]))
            except Exception:
                pass
    fh = open(jsonl, "a"); nok = 0; t0 = time.time()
    for j, Bval in enumerate(B_axis):
        for i, sval in enumerate(sig_axis):
            if (i, j) in done:
                continue
            t = time.time()
            ch = _node(ds, mode, anchor, sval, Bval, precomp, salt)
            rec = {"i": i, "j": j, "sigma": float(sval), "B": float(Bval)}
            if ch is None:
                rec["status"] = "infeasible"
            else:
                rec["status"] = "ok"; rec.update({k: ch[k] for k in CH_KEYS}); nok += 1
            rec["dt"] = round(time.time() - t, 1)
            fh.write(json.dumps(rec) + "\n"); fh.flush()
    fh.close()
    (od / "WORKER_DONE").write_text("done")
    _export_csv(rid); plot_exp(rid)
    print(f"[exp] {rid}: {nok} feasible nodes  Lp*={Lp_fixed:.2f}  salt={salt}  "
          f"{time.time()-t0:.0f}s -> {od/'contour7.png'}", flush=True)
    return od


def _load_Z(od, sig_axis, B_axis):
    import numpy as np
    n = len(sig_axis); m = len(B_axis)
    Z = {k: np.full((m, n), np.nan) for k in CH_KEYS}
    jsonl = od / "nodes.jsonl"
    if jsonl.exists():
        for ln in jsonl.read_text().splitlines():
            try:
                r = json.loads(ln)
                if r.get("status") == "ok":
                    for k in CH_KEYS:
                        v = r.get(k)
                        Z[k][r["j"], r["i"]] = v if v is not None else np.nan
            except Exception:
                pass
    return Z


def _export_csv(rid):
    import numpy as np
    od = out_dir(rid)
    g = json.loads((od / "grid.json").read_text())
    rows = []
    for ln in (od / "nodes.jsonl").read_text().splitlines():
        try:
            r = json.loads(ln)
        except Exception:
            continue
        if r.get("status") != "ok":
            continue
        def lg(k):
            v = r.get(k)
            return float(np.log10(v)) if (v and v > 0 and np.isfinite(v)) else float("nan")
        rows.append([r["sigma"], r["B"]] + [lg(k) for k in CH_KEYS])
    rows.sort(key=lambda t: (t[1], t[0]))
    hdr = ["sigma", "B"] + [f"log10_{k}" for k in CH_KEYS]
    lines = [",".join(hdr)] + [",".join(f"{v:.6g}" for v in row) for row in rows]
    (od / "contourdata-7ch.csv").write_text("\n".join(lines) + "\n")


def plot_exp(rid):
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    od = out_dir(rid)
    g = json.loads((od / "grid.json").read_text())
    A = np.array(g["axis_sigma"]); Bv = np.array(g["axis_B"])
    Z = _load_Z(od, A.tolist(), Bv.tolist())
    S, B = np.meshgrid(A, Bv)
    label = next((k for k, v in REPS.items() if v == rid), rid)

    fig, axes = plt.subplots(2, 4, figsize=(20, 9.2))
    axflat = axes.flatten()
    for n, (ax, (key, title, unit)) in enumerate(zip(axflat, CHANNELS)):
        _panel(ax, S, B, Z[key], f"{n + 1}. {title}", unit)
    axflat[-1].axis("off")                       # 8th cell unused (7 channels)
    for r in range(2):
        axes[r, 0].set_ylabel("B  [µm s⁻¹]")
    # legend / reading guide in the unused 8th cell
    axflat[-1].text(0.02, 0.95,
                    "Square forward-solve SLICE\n(no per-node optimization):\n"
                    "L$_p$ pinned at its optimum; at each\n(σ, B) node the DAE is solved once\n"
                    "and every channel's WSSE evaluated.\n\n"
                    "▲  per-channel optimum (min WSSE)\n"
                    "lines = log₁₀ WSSE level curves\n\n"
                    "Channels 4–7 route concentration↔\nconductivity through the Shedlovsky\n"
                    "converter (single-salt).",
                    transform=axflat[-1].transAxes, va="top", ha="left", fontsize=9,
                    color="0.15")
    fig.suptitle(f"{label} ({rid.split('_')[0]}) — σ×B identifiability across 7 response "
                 f"channels · L$_p$ fixed at its optimum {g['Lp_fixed']:.2f} "
                 f"L m⁻² h⁻¹ bar⁻¹ · square forward-solve slice (log₁₀ WSSE)",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.subplots_adjust(hspace=0.42)
    fig.savefig(od / "contour7.png", dpi=155); plt.close(fig)


def _panel(ax, S, B, Zk, title, unit):
    """v7 titled labelled-line-contour panel (turbo, red ▲ at the minimum)."""
    import numpy as np
    ax.set_facecolor("white")
    ax.set_title(f"{title}\n[{unit}]", fontsize=10, fontweight="bold")
    ax.set_xlabel("σ")
    if Zk is None or not np.isfinite(Zk).any():
        ax.text(0.5, 0.5, "(no data)", transform=ax.transAxes, ha="center", color="0.5")
        return
    Zl = np.log10(np.clip(Zk, 1e-12, None))
    levels = np.linspace(np.nanmin(Zl), np.nanmax(Zl), 14)
    try:
        cs = ax.contour(S, B, Zl, levels=levels, cmap="turbo", linewidths=1.1)
        ax.clabel(cs, inline=True, fontsize=6, fmt="%.2f")
    except Exception:
        pass
    ax.grid(alpha=0.2, lw=0.4)
    j, i = np.unravel_index(np.nanargmin(Zk), Zk.shape)
    ax.plot(S[j, i], B[j, i], "^", color="red", ms=11, mec="k", mew=0.7, zorder=6)
    ax.text(0.03, 0.96, f"min: σ={S[j, i]:.2f}, B={B[j, i]:.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=7.5, color="red",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="red", lw=0.5, alpha=0.85))


def main():
    if len(sys.argv) < 2:
        print(__doc__); return
    cmd = sys.argv[1]
    if cmd == "exp":
        run_exp(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 11); return
    if cmd == "all":
        density = int(sys.argv[2]) if len(sys.argv) > 2 else 11
        for rid in REPS.values():
            run_exp(rid, density)
        return
    if cmd == "plot":
        plot_exp(sys.argv[2]); return
    print(__doc__)


if __name__ == "__main__":
    main()
