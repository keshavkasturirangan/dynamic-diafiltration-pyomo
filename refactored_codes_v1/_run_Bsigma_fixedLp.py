#!/usr/bin/env python3
"""DATA3 B–σ identifiability contour at FIXED Lp.

The user's directive: Lp is consistently WELL-IDENTIFIED across DATA3 (same NF270
membrane → Lp stable; see the pooling/identifiability finding).  So pin Lp at its fitted
optimum Lp* and map the objective over (σ, B) — exposing the residual σ–B correlation
that the model trades off, WITHOUT Lp absorbing the variation.  This is the classic
DATA1 identifiability triangle with the easy axis (Lp) removed.

Method = SLICE: fix Lp=Lp*, set σ and B at each grid node, forward-simulate the transport
model, read the per-channel objective (mass | permeate | retentate | combined WSSE3).
B is swept by scaling the form's magnitude coefficient(s) so the realized membrane B at
the feed concentration (read back from the library's own sim) equals the axis value.

Forms: the slice is exact + fully feasible for the UNBOUNDED-B forms
  single  : B is a scalar Param            → axis = B [µm/s] (σ-independent, exact)
  sat     : B = B_inf·(1−e^{−c/c*}), B Var unbounded → axis = B@feed [µm/s]
  donnan  : B from P0,X,k_dd, B Var unbounded       → axis = B@feed [µm/s]
The POLYNOMIAL forms (poly1/2/3) use a BOUNDED B Var(1e-6,50); off-optimum the curve hits
the bound and the square DAE goes locally-infeasible, so a pure slice cannot fill them
(verified empirically).  They are skipped here — their identifiability picture is the
existing σ×Lp PROFILE contour (B(c) re-solved), in _run_profile_contours.py.

Output (bform_study/Bsigma_fixedLp/):
  <rid>/<form>/nodes.jsonl,grid.json,contour.png   per-experiment
  _stack/<group>/<form>/stacked.png,stacked.json   summed group landscape (pooled)

Usage:
  python3 _run_Bsigma_fixedLp.py exp   <rid> <form> [density]     # one experiment
  python3 _run_Bsigma_fixedLp.py stack <group> <form> [density]   # one group (per-exp + sum)
  python3 _run_Bsigma_fixedLp.py campaign [density] [maxpar]      # all feasible forms, all exps+groups
"""
import os, sys, json, time, io, contextlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _run_profile_contours as pc   # reuse config + the v7 line-contour panel

# forms whose membrane B is unbounded → the slice is feasible; axis = B@feed [µm/s]
SLICE_FORMS = ["single", "sat", "donnan"]
# magnitude coefficient(s) to scale so B(c) ∝ k (uniform scale of the whole curve)
MAG = {"single": ["B"], "sat": ["B_inf"], "donnan": ["P0"]}
FORMARG = {f: pc.FORMS[f][0] for f in pc.FORMS}            # token → B_form arg
TAG = {f: pc.FORMS[f][1] for f in pc.FORMS}               # token → result json tag
CPU = 30                                                   # IPOPT cpu cap per node (s)

# ---- retentate-weighting variant (Lilonfe et al.: cF = 2% relative) ---------------
# "2pct_nofloor" (DEFAULT, canonical): pure 2%, NO absolute floor — byte-consistent
#   with the canonical fits, the 3D grids, and the DATA1 σ×Lp/B×σ panels.
# "2pct_floor1mm": 2% with a 1 mM absolute floor — the earlier (floored) variant,
#   kept SIDE-BY-SIDE for comparison.  Select via env BSIGMA_VARIANT.
WEIGHT_VARIANT = os.environ.get("BSIGMA_VARIANT", "2pct_nofloor")
CF_FRAC  = 0.02
CF_FLOOR = None if WEIGHT_VARIANT == "2pct_nofloor" else 1.0
OUT_SUBDIR = "Bsigma_fixedLp_nofloor" if WEIGHT_VARIANT == "2pct_nofloor" else "Bsigma_fixedLp"


def _set_weighting(lib):
    """Pin the retentate weighting for THIS variant (explicit; no module drift)."""
    lib.NF270_CF_RESIDUAL_SCALE_FRACTION = CF_FRAC
    lib.NF270_CF_RESIDUAL_FLOOR_MM = CF_FLOOR


def _study():
    return pc._study()


def out_dir(rid, form):
    return _study() / OUT_SUBDIR / rid / form


def _anchor(rid, form):
    rp = _study() / rid / f"result_{TAG[form]}.json"
    if not rp.exists():
        return None
    r = json.loads(rp.read_text())
    return r.get("parameters")


# --------------------------------------------------------------------------- solve
def _slice_node(ds, mode, anchor, form, sigma_val, k):
    """Forward-sim with Lp=Lp* (anchor), σ=sigma_val, B-curve scaled by k.  Returns
    (Bfeed, obj dict) or (None, None) if the node is infeasible."""
    import refactored_ucb_library as lib
    seed = dict(anchor)
    seed["sigma"] = sigma_val
    for mp in MAG[form]:
        if mp in seed:
            seed[mp] = seed[mp] * k
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fit, sim, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=True, B_form=FORMARG[form],
                                          workflow_family="DATA3", nfe=120, solver_max_cpu_time=CPU)
        if not isinstance(fit, dict) or sim is None:
            return None, None
        Bfeed = float(sim[0]["B"][0])
        obj = {"obj_m": float(fit.get("obj_m")), "obj_cv": float(fit.get("obj_cv")),
               "obj_cr": float(fit.get("obj_cr"))}
        obj["WSSE3"] = obj["obj_m"] + obj["obj_cv"] + obj["obj_cr"]
        return Bfeed, obj
    except Exception:
        return None, None


def _bfeed_fit(ds, mode, anchor, form):
    Bfeed, obj = _slice_node(ds, mode, anchor, form, anchor.get("sigma", 1.0), 1.0)
    return Bfeed


# --------------------------------------------------------------------------- per-experiment
def run_exp(rid, form, density=9, B_axis=None):
    """Compute the (σ, B@feed) slice landscape for one experiment with Lp pinned.
    If B_axis is given (common grid for stacking) use it; else build it around B@feed_fit."""
    import numpy as np
    import _run_bform_campaign as camp
    import refactored_ucb_library as lib
    _set_weighting(lib)
    anchor = _anchor(rid, form)
    if anchor is None:
        print(f"[exp] {rid} {form}: no converged fit — skip", flush=True)
        return None
    ds = camp.load(rid); mode = ds["mode"]
    Bfit = _bfeed_fit(ds, mode, anchor, form)
    if Bfit is None or not np.isfinite(Bfit) or Bfit <= 0:
        print(f"[exp] {rid} {form}: fit sim infeasible — skip", flush=True)
        return None
    sig_axis = list(np.linspace(0.0, 1.0, density))
    if B_axis is None:
        B_axis = list(np.linspace(max(1e-4, 0.3 * Bfit), 1.7 * Bfit, density))
    B_axis = [float(b) for b in B_axis]
    od = out_dir(rid, form); od.mkdir(parents=True, exist_ok=True)
    jsonl = od / "nodes.jsonl"
    # idempotency: if a prior run used a DIFFERENT grid, its (i,j) nodes are stale -> wipe.
    gridf = od / "grid.json"
    if gridf.exists():
        old = json.loads(gridf.read_text())
        same = (old.get("axis_sigma") == sig_axis and old.get("axis_B") == B_axis)
        if same and (od / "WORKER_DONE").exists():
            print(f"[exp] {rid} {form}: already done (grid match) — skip", flush=True)
            return od
        if not same:
            jsonl.unlink(missing_ok=True); (od / "WORKER_DONE").unlink(missing_ok=True)
    grid = {"rid": rid, "form": form, "Lp_fixed": float(anchor["Lp"]), "Bfeed_fit": float(Bfit),
            "sigma_fit": float(anchor.get("sigma", 1.0)), "axis_sigma": sig_axis, "axis_B": list(B_axis)}
    gridf.write_text(json.dumps(grid))
    done = set()
    if jsonl.exists():
        for ln in jsonl.read_text().splitlines():
            try:
                r = json.loads(ln); done.add((r["i"], r["j"]))
            except Exception:
                pass
    fh = open(jsonl, "a")
    nok = 0
    for j, Bval in enumerate(B_axis):           # y = B@feed
        for i, sval in enumerate(sig_axis):     # x = sigma
            if (i, j) in done:
                continue
            t = time.time()
            Bfeed, obj = _slice_node(ds, mode, anchor, form, float(sval), float(Bval) / Bfit)
            rec = {"i": i, "j": j, "sigma": float(sval), "B": float(Bval)}
            if obj is None:
                rec["status"] = "infeasible"
            else:
                rec["status"] = "ok"; rec["Bfeed"] = Bfeed; rec.update(obj); nok += 1
            rec["dt"] = round(time.time() - t, 1)
            fh.write(json.dumps(rec) + "\n"); fh.flush()
    fh.close()
    (od / "WORKER_DONE").write_text("done")
    plot_exp(rid, form)
    _export_csv(rid, form)
    print(f"[exp] {rid} {form}: {nok} feasible nodes  Lp*={anchor['Lp']:.2f}  B@feed_fit={Bfit:.3f}  [{WEIGHT_VARIANT}] -> {od/'contour.png'}", flush=True)
    return od


def _load_Z(od, density):
    import numpy as np
    Z = {k: np.full((density, density), np.nan) for k in ("obj_m", "obj_cv", "obj_cr", "WSSE3")}
    jsonl = od / "nodes.jsonl"
    if jsonl.exists():
        for ln in jsonl.read_text().splitlines():
            try:
                r = json.loads(ln)
                if r.get("status") == "ok":
                    for k in Z:
                        Z[k][r["j"], r["i"]] = r[k]
            except Exception:
                pass
    return Z


def _export_csv(rid, form):
    """Write the slice as a DATA1-schema contour CSV (log10 objectives) for quick retrieval:
    columns [sigma, B, Obj_mass, Obj_concentration, Obj_retentate_concentration]."""
    import numpy as np
    od = out_dir(rid, form)
    jsonl = od / "nodes.jsonl"
    if not jsonl.exists():
        return
    rows = []
    for ln in jsonl.read_text().splitlines():
        try:
            r = json.loads(ln)
        except Exception:
            continue
        if r.get("status") != "ok":
            continue
        def lg(k):
            v = r.get(k)
            return float(np.log10(v)) if (v and v > 0) else float("nan")
        rows.append((r["sigma"], r["B"], lg("obj_m"), lg("obj_cv"), lg("obj_cr")))
    rows.sort(key=lambda t: (t[1], t[0]))
    lines = ["sigma,B,Obj_mass,Obj_concentration,Obj_retentate_concentration"]
    lines += [",".join(f"{v:.6g}" for v in row) for row in rows]
    (od / "contourdata-x_sigma-y_B.csv").write_text("\n".join(lines) + "\n")


def plot_exp(rid, form):
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    od = out_dir(rid, form)
    grid = json.loads((od / "grid.json").read_text())
    A = np.array(grid["axis_sigma"]); Bv = np.array(grid["axis_B"]); n = len(A)
    Z = _load_Z(od, n)
    X, Y = np.meshgrid(A, Bv)
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))
    panels = [("obj_m", "Mass Residual$^2$ [g$^2$]"), ("obj_cv", "Permeate-conc Residual$^2$ [mM$^2$]"),
              ("obj_cr", "Retentate-conc Residual$^2$ [mM$^2$]"), ("WSSE3", "Combined WSSE$_3$")]
    for ax, (key, title) in zip(axes, panels):
        pc._contour_panel(ax, X, Y, Z[key], A, Bv, title, "$\\sigma$ [-]", "B@feed [$\\mu$m s$^{-1}$]")
    fig.suptitle(f"{rid}  ·  B-form: {form}  ·  B–σ slice at FIXED L$_p$={grid['Lp_fixed']:.2f} "
                 f"(well-identified)  ·  B@feed(fit)={grid['Bfeed_fit']:.3f}", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(od / "contour.png", dpi=150); plt.close(fig)


# --------------------------------------------------------------------------- stacked (pooled)
def run_stack(group, form, density=9):
    """Per-experiment slices on a COMMON (σ, B@feed) grid, then SUM (pooled objective =
    summed Fisher information).  Sharper pooled min vs scattered per-exp minima = the
    experiments informing each other on (σ, B) once Lp is pinned."""
    import numpy as np
    import _run_bform_campaign as camp
    import refactored_ucb_library as lib
    _set_weighting(lib)
    sheets = pc.POOL_GROUPS[group]
    fits = {}
    for s in sheets:
        a = _anchor(s, form)
        if a is None:
            continue
        ds = camp.load(s); mode = ds["mode"]
        bf = _bfeed_fit(ds, mode, a, form)
        if bf and np.isfinite(bf) and bf > 0:
            fits[s] = bf
    if not fits:
        print(f"[stack] {group} {form}: no usable fits — skip", flush=True)
        return
    lo = min(fits.values()) * 0.3; hi = max(fits.values()) * 1.7
    B_axis = list(np.linspace(max(1e-4, lo), hi, density))
    for s in fits:
        run_exp(s, form, density, B_axis=B_axis)
    _stackplot(group, form, density, B_axis, list(fits))


def _stackplot(group, form, density, B_axis, sheets):
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    A = np.linspace(0.0, 1.0, density); Bv = np.array(B_axis); X, Y = np.meshgrid(A, Bv)
    keys = ["obj_m", "obj_cv", "obj_cr", "WSSE3"]
    stackZ = {k: np.zeros((density, density)) for k in keys}
    cov = np.zeros((density, density)); per_exp_min = []
    for s in sheets:
        Z = _load_Z(out_dir(s, form), density)
        w = Z["WSSE3"]
        if np.all(np.isnan(w)):
            continue
        jm, im = np.unravel_index(np.nanargmin(w), w.shape)
        per_exp_min.append((float(A[im]), float(Bv[jm]), s.split("_")[0]))
        for k in keys:
            m = ~np.isnan(Z[k]); stackZ[k][m] += Z[k][m]
        cov += (~np.isnan(w)).astype(float)
    full = cov >= len(sheets)
    for k in keys:
        stackZ[k][~full] = np.nan
    comb = stackZ["WSSE3"]; star = None
    if np.any(~np.isnan(comb)):
        jm, im = np.unravel_index(np.nanargmin(comb), comb.shape)
        star = (float(A[im]), float(Bv[jm]), float(comb[jm, im]))
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))
    titles = {"obj_m": "Σ mass", "obj_cv": "Σ permeate-conc", "obj_cr": "Σ retentate-conc", "WSSE3": "Σ combined WSSE$_3$"}
    for ax, k in zip(axes, keys):
        pc._contour_panel(ax, X, Y, stackZ[k], A, Bv, titles[k], "$\\sigma$ [-]", "B@feed [$\\mu$m s$^{-1}$]")
        for (xx, yy, nm) in per_exp_min:
            ax.plot(xx, yy, "o", ms=6, mfc="white", mec="0.2", zorder=5)
    sub = f"pooled min ★ σ={star[0]:.2f}, B@feed={star[1]:.2f}" if star else "no full-overlap region"
    fig.suptitle(f"STACKED σ×B at FIXED L$_p$ — {group} · B-form={form}: {len(sheets)} experiments summed.   "
                 f"○ = each experiment's own min;  ▲ = pooled min / joint-fit seed.   {sub}", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    root = _study() / OUT_SUBDIR / "_stack" / group / form
    root.mkdir(parents=True, exist_ok=True)
    fig.savefig(root / "stacked.png", dpi=150); plt.close(fig)
    (root / "stacked.json").write_text(json.dumps(
        {"group": group, "form": form, "Lp": "fixed at per-exp optimum",
         "pooled_min": ({"sigma": star[0], "B": star[1], "WSSE3": star[2]} if star else None),
         "per_exp_min": [{"sigma": x, "B": y, "exp": n} for (x, y, n) in per_exp_min]}, indent=2))
    print(f"[stack] {group} {form}: {'min σ=%.2f B=%.2f' % (star[0], star[1]) if star else 'no overlap'} -> {root/'stacked.png'}", flush=True)


def main():
    if len(sys.argv) < 2:
        print(__doc__); return
    cmd = sys.argv[1]
    if cmd == "exp":
        run_exp(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 9); return
    if cmd == "stack":
        run_stack(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 9); return
    if cmd == "campaign":
        density = int(sys.argv[2]) if len(sys.argv) > 2 else 9
        for form in SLICE_FORMS:
            for group in pc.POOL_GROUPS:
                run_stack(group, form, density)
        return
    if cmd == "single":
        # DATA3 single-salt scope: constant-B form only, all 11 sheets (per-exp via the
        # 4 salt×regime groups) + the 4 pooled stacks.  Variant via env BSIGMA_VARIANT.
        density = int(sys.argv[2]) if len(sys.argv) > 2 else 11
        print(f"[single] variant={WEIGHT_VARIANT}  cf_frac={CF_FRAC} floor={CF_FLOOR}  -> {OUT_SUBDIR}/", flush=True)
        for group in pc.POOL_GROUPS:
            run_stack(group, "single", density)
        return
    print(__doc__)


if __name__ == "__main__":
    main()
