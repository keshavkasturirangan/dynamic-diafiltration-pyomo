#!/usr/bin/env python3
"""Task 8 (B(c) forms) — square coefficient × σ WSSE identifiability contours, fixed Lp.

The constant-B contours (`_run_Bsigma_channels.py`) plotted B × σ.  For a B(c)
correlation the single "B" axis is replaced by ONE axis per shape coefficient:

  linear      B(c)=β₀+β₁c                 → β₀×σ, β₁×σ
  quadratic   B(c)=β₀+β₁c+β₂c²            → β₀×σ, β₁×σ, β₂×σ
  cubic       B(c)=β₀+β₁c+β₂c²+β₃c³       → β₀×σ, β₁×σ, β₂×σ, β₃×σ
  saturating  B(c)=B∞(1−e^(−c/c*))        → B∞×σ, c*×σ
                                       (c = INTERFACIAL concentration cIn)

METHOD = SQUARE SLICE (no per-node optimization).  At each (coeff, σ) node:
  • Lp is pinned at its fitted optimum,
  • the plotted coefficient and σ are set to the node values,
  • every OTHER coefficient (and S0, S) is held at the form's WARM-START fit θ*,
so the model has zero free parameters → the DAE is forward-solved ONCE (square)
and all seven response-channel WSSE are evaluated (same 7 channels as
`_run_Bsigma_channels`).  This is the slice / conditional landscape through θ*.

Robustness (these square solves are stiff off-optimum):
  • lib.NF270_BFORM_SLICE_UNBOUNDED_B = True — the realized B Var is unbounded
    in sim_opt (matches sat/donnan); the [1e-6,50] box otherwise makes the DAE
    locally infeasible when an off-optimum coefficient pin drives B(cIn) out.
  • a cheap PRE-FILTER skips nodes whose implied B(cIn) (on the warm cIn grid)
    leaves the physical window — avoids most stiff grinds.
  • a subprocess WORKER marches the grid (warm-start order, near the optimum
    first) writing each node to JSONL; the PARENT watchdog kills+restarts the
    worker if a node stalls past PER_NODE_TIMEOUT (records it as a hole) and the
    worker resumes from the JSONL.  Holes are expected far from the optimum and
    are themselves informative (the coefficient is feasibility-constrained).

Output (bform_study/Bsigma_coeffs/<rid>/<form>/<coeff>/): nodes.jsonl, grid.json
       bform_study/Bsigma_coeffs/<rid>/<form>/contour_grid.png  (coeff rows × 7 channels)

Usage:
  python3 _run_Bsigma_coeffs.py sweep <rid> <form> <coeff> [density]   # one coeff axis
  python3 _run_Bsigma_coeffs.py expform <rid> <form> [density]         # all coeffs of a form
  python3 _run_Bsigma_coeffs.py all [density] [maxpar]                 # 3 exp × 4 forms
  python3 _run_Bsigma_coeffs.py plot <rid> <form>                      # render composite
  python3 _run_Bsigma_coeffs.py --worker <gridpath>                    # internal
"""
import os, sys, json, time, signal, subprocess
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _run_profile_contours as pc
import _run_Bsigma_channels as ch   # 7-channel machinery (precompute, _seven_channels, panel, CHANNELS)

NFE = 120
IPOPT_CPU = 18
PER_NODE_TIMEOUT = 25
B_LO, B_HI = 1e-4, 250.0            # physical pre-filter window on B(cIn) [µm s⁻¹ partition]

FORMS = {
    "poly1": {"Bform": 1,     "tag": "poly1", "coeffs": ["beta_0", "beta_1"]},
    "poly2": {"Bform": 2,     "tag": "poly2", "coeffs": ["beta_0", "beta_1", "beta_2"]},
    "poly3": {"Bform": 3,     "tag": "poly3", "coeffs": ["beta_0", "beta_1", "beta_2", "beta_3"]},
    "sat":   {"Bform": "sat", "tag": "sat",   "coeffs": ["B_inf", "c_star"]},
}
FORM_EQ = {
    "poly1": "B(c) = β₀ + β₁·c",
    "poly2": "B(c) = β₀ + β₁·c + β₂·c²",
    "poly3": "B(c) = β₀ + β₁·c + β₂·c² + β₃·c³",
    "sat":   "B(c) = B∞·(1 − e^(−c/c*))",
}
COEFF_LABEL = {"beta_0": "β₀", "beta_1": "β₁", "beta_2": "β₂", "beta_3": "β₃",
               "B_inf": "B∞", "c_star": "c*"}
# positivity-constrained coefficients (axis kept ≥ small positive)
POS_COEFF = {"beta_0", "B_inf", "c_star"}
# per-coefficient absolute half-width floor (so near-zero coefficients still span)
FLOOR = {"beta_0": 0.4, "beta_1": 0.04, "beta_2": 0.004, "beta_3": 0.0004,
         "B_inf": 0.4, "c_star": 0.5}
REPS = ["MC2.05.07.24_CaCl2", "MC2.05.21.24_LaCl3", "MC2.05.07.24_NaCl"]


def out_dir(rid, form, coeff=None):
    d = pc._study() / "Bsigma_coeffs" / rid / form
    return d / coeff if coeff else d


def _warm(rid, form):
    rp = pc._study() / rid / f"result_{FORMS[form]['tag']}.json"
    if not rp.exists():
        return None
    r = json.loads(rp.read_text())
    return r.get("parameters")


def _coeff_axis(name, warm, density):
    import numpy as np
    w = float(warm[name])
    if name == "c_star":
        h = max(0.6 * abs(w), FLOOR[name])
    else:
        h = max(abs(w), FLOOR[name])
    lo, hi = w - h, w + h
    if name in POS_COEFF:
        lo = max(lo, 1e-4)
    return [float(v) for v in np.linspace(lo, hi, density)]


def _Bpoly(form, warm, coeff_name, coeff_val, c):
    """Implied B(c) for a node's coefficient set (others at warm), on conc grid c."""
    import numpy as np
    p = dict(warm); p[coeff_name] = coeff_val
    c = np.asarray(c, float)
    if form == "sat":
        return p["B_inf"] * (1.0 - np.exp(-c / max(p["c_star"], 1e-9)))
    B = np.zeros_like(c) + p.get("beta_0", 0.0)
    for k, nm in ((1, "beta_1"), (2, "beta_2"), (3, "beta_3")):
        if nm in p:
            B = B + p[nm] * c ** k
    return B


# ----------------------------------------------------------------------- worker
def _worker(gridpath):
    import contextlib, io
    import numpy as np
    import refactored_ucb_library as lib
    import _run_bform_campaign as camp
    ch._set_weighting(lib)
    lib.NF270_BFORM_SLICE_UNBOUNDED_B = True

    grid = json.loads(Path(gridpath).read_text())
    rid, form, coeff = grid["rid"], grid["form"], grid["coeff"]
    Bform = FORMS[form]["Bform"]
    warm = grid["warm"]; warm = dict(warm); warm["Lp"] = grid["Lp_fixed"]
    A = np.array(grid["axis_sigma"]); Bv = np.array(grid["axis_coeff"])
    ai, aj = grid["anchor_ij"]
    od = Path(grid["out_dir"]); od.mkdir(parents=True, exist_ok=True)
    jsonl = od / "nodes.jsonl"; curf = od / "nodes.jsonl.cur"; skipf = od / "nodes.jsonl.skip"

    ds = camp.load(rid); mode = ds["mode"]
    salt = str(ds["data_config"].get("namec", "")).strip()
    precomp = ch._precompute(ds, salt)
    # warm cIn grid for the cheap physical pre-filter
    cIn_warm = np.asarray(grid.get("cIn_warm", []), float)

    done, skip = set(), set()
    if jsonl.exists():
        for ln in jsonl.read_text().splitlines():
            try:
                r = json.loads(ln); done.add((r["i"], r["j"]))
            except Exception:
                pass
    if skipf.exists():
        for ln in skipf.read_text().splitlines():
            try:
                a, b = ln.split(","); skip.add((int(a), int(b)))
            except Exception:
                pass

    nA, nB = len(A), len(Bv)
    order = sorted([(i, j) for i in range(nA) for j in range(nB)],
                   key=lambda ij: (ij[0] - ai) ** 2 + (ij[1] - aj) ** 2)
    fh = open(jsonl, "a")
    for (i, j) in order:
        if (i, j) in done or (i, j) in skip:
            continue
        curf.write_text(f"{i},{j}")
        sval, cval = float(A[i]), float(Bv[j])
        rec = {"i": i, "j": j, "sigma": sval, "coeff": cval}
        # pre-filter: implied B(cIn_warm) must stay physical
        if cIn_warm.size:
            Bc = _Bpoly(form, warm, coeff, cval, cIn_warm)
            if not (np.isfinite(Bc).all() and Bc.min() > B_LO and Bc.max() < B_HI):
                rec["status"] = "unphysical"
                fh.write(json.dumps(rec) + "\n"); fh.flush(); os.fsync(fh.fileno())
                continue
        seed = dict(warm); seed["sigma"] = sval; seed[coeff] = cval
        t = time.time()
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                fit, sim, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=True, B_form=Bform,
                                              workflow_family="DATA3", nfe=NFE,
                                              solver_max_cpu_time=IPOPT_CPU)
            if isinstance(fit, dict) and sim is not None:
                cc = ch._seven_channels(fit, sim, precomp, salt)
                rec["status"] = "ok"; rec.update({k: cc[k] for k in ch.CH_KEYS})
            else:
                rec["status"] = "infeasible"
        except Exception as e:
            rec["status"] = "error"; rec["err"] = str(e)[:60]
        rec["dt"] = round(time.time() - t, 1)
        fh.write(json.dumps(rec) + "\n"); fh.flush(); os.fsync(fh.fileno())
    fh.close()
    if curf.exists():
        curf.unlink()
    (od / "WORKER_DONE").write_text("done")


# ----------------------------------------------------------------------- parent
def _launch_worker(gridpath):
    return subprocess.Popen([sys.executable, str(HERE / "_run_Bsigma_coeffs.py"),
                             "--worker", str(gridpath)], start_new_session=True, cwd=str(HERE))


def sweep(rid, form, coeff, density=9):
    import numpy as np
    import contextlib, io
    import refactored_ucb_library as lib
    import _run_bform_campaign as camp
    warm = _warm(rid, form)
    if warm is None or coeff not in warm:
        print(f"[sweep] {rid} {form} {coeff}: no warm-start — skip", flush=True); return None
    od = out_dir(rid, form, coeff); od.mkdir(parents=True, exist_ok=True)
    sig_axis = list(np.linspace(0.0, 1.0, density))
    coeff_axis = _coeff_axis(coeff, warm, density)
    ai = int(np.argmin(np.abs(np.array(sig_axis) - warm.get("sigma", 1.0))))
    aj = int(np.argmin(np.abs(np.array(coeff_axis) - warm[coeff])))
    # warm cIn trajectory for the pre-filter
    ch._set_weighting(lib); lib.NF270_BFORM_SLICE_UNBOUNDED_B = True
    ds = camp.load(rid); mode = ds["mode"]
    cIn_warm = []
    try:
        seed = dict(warm); seed["Lp"] = warm["Lp"]
        with contextlib.redirect_stdout(io.StringIO()):
            fw, sw, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=True,
                                        B_form=FORMS[form]["Bform"], workflow_family="DATA3",
                                        nfe=NFE, solver_max_cpu_time=IPOPT_CPU)
        if sw is not None:
            cv = np.concatenate([np.asarray(sw[k]["cIn"], float).reshape(-1) for k in sw])
            cIn_warm = [float(x) for x in cv[np.isfinite(cv)]]
    except Exception:
        cIn_warm = []
    grid = {"rid": rid, "form": form, "coeff": coeff, "Lp_fixed": float(warm["Lp"]),
            "warm": {k: float(v) for k, v in warm.items()},
            "axis_sigma": sig_axis, "axis_coeff": coeff_axis, "anchor_ij": [ai, aj],
            "cIn_warm": cIn_warm, "out_dir": str(od)}
    (od / "grid.json").write_text(json.dumps(grid))

    jsonl = od / "nodes.jsonl"; curf = od / "nodes.jsonl.cur"; skipf = od / "nodes.jsonl.skip"
    donef = od / "WORKER_DONE"; donef.unlink(missing_ok=True)
    ntot = density * density
    print(f"[sweep] {rid} {form} {coeff}: {ntot} nodes (timeout {PER_NODE_TIMEOUT}s/node)", flush=True)
    t0 = time.time()
    while not donef.exists():
        nd = sum(1 for _ in jsonl.open()) if jsonl.exists() else 0
        ns = sum(1 for _ in skipf.open()) if skipf.exists() else 0
        if nd + ns >= ntot:
            break
        proc = _launch_worker(od / "grid.json")
        last, lastp = nd, time.time()
        while proc.poll() is None:
            time.sleep(2)
            cur = sum(1 for _ in jsonl.open()) if jsonl.exists() else 0
            if cur > last:
                last, lastp = cur, time.time()
            elif time.time() - lastp > PER_NODE_TIMEOUT:
                stuck = curf.read_text().strip() if curf.exists() else None
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except Exception:
                    pass
                proc.wait()
                if stuck:
                    with open(skipf, "a") as sf:
                        sf.write(stuck + "\n")
                break
        else:
            proc.wait()
    nd = sum(1 for _ in jsonl.open()) if jsonl.exists() else 0
    nok = sum(1 for ln in jsonl.open() if '"status": "ok"' in ln) if jsonl.exists() else 0
    print(f"[sweep] {rid} {form} {coeff}: {nok} feasible / {nd} nodes in {time.time()-t0:.0f}s", flush=True)
    return od


def expform(rid, form, density=9):
    for coeff in FORMS[form]["coeffs"]:
        sweep(rid, form, coeff, density)
    plot_expform(rid, form)


def all_runs(density=9, maxpar=6, rids=None, forms=None):
    rids = rids or REPS
    forms = forms or list(FORMS)
    jobs = [(rid, form, coeff) for rid in rids for form in forms for coeff in FORMS[form]["coeffs"]]
    todo = [j for j in jobs if not (out_dir(*j) / "WORKER_DONE").exists()]
    print(f"[all] {len(jobs)} coeff-sweeps, {len(todo)} to run, {maxpar}-way parallel", flush=True)
    running = {}
    q = list(todo)
    while q or running:
        while q and len(running) < maxpar:
            rid, form, coeff = q.pop(0)
            log = out_dir(rid, form, coeff); log.mkdir(parents=True, exist_ok=True)
            fh = open(log / "run.log", "w")
            p = subprocess.Popen([sys.executable, str(HERE / "_run_Bsigma_coeffs.py"),
                                  "sweep", rid, form, coeff, str(density)],
                                 stdout=fh, stderr=subprocess.STDOUT, cwd=str(HERE))
            running[p] = (rid, form, coeff, fh)
            print(f"  launch {rid} {form} {coeff} ({len(running)} running, {len(q)} queued)", flush=True)
        done = [p for p in running if p.poll() is not None]
        for p in done:
            rid, form, coeff, fh = running.pop(p); fh.close()
            print(f"  done   {rid} {form} {coeff}", flush=True)
        if not done:
            time.sleep(4)
    for rid in rids:
        for form in forms:
            try:
                plot_expform(rid, form)
            except Exception as e:
                print(f"[plot] {rid} {form} FAILED: {str(e)[:80]}", flush=True)
    print("[all] COMPLETE", flush=True)


# -------------------------------------------------------------------- plotting
def _load_Z(od, n_sig, n_coeff):
    import numpy as np
    Z = {k: np.full((n_coeff, n_sig), np.nan) for k in ch.CH_KEYS}
    jsonl = od / "nodes.jsonl"
    if jsonl.exists():
        for ln in jsonl.read_text().splitlines():
            try:
                r = json.loads(ln)
                if r.get("status") == "ok":
                    for k in ch.CH_KEYS:
                        v = r.get(k)
                        Z[k][r["j"], r["i"]] = v if (v is not None) else np.nan
            except Exception:
                pass
    return Z


def plot_expform(rid, form):
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    coeffs = FORMS[form]["coeffs"]
    grids = {}
    for cf in coeffs:
        gp = out_dir(rid, form, cf) / "grid.json"
        if gp.exists():
            grids[cf] = json.loads(gp.read_text())
    if not grids:
        print(f"[plot] {rid} {form}: no grids", flush=True); return
    nrows, ncols = len(coeffs), len(ch.CHANNELS)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.7 * nrows), squeeze=False)
    Lp = next(iter(grids.values()))["Lp_fixed"]
    warm = next(iter(grids.values()))["warm"]
    for r, cf in enumerate(coeffs):
        if cf not in grids:
            for c in range(ncols):
                axes[r, c].axis("off")
            continue
        g = grids[cf]
        A = np.array(g["axis_sigma"]); Cv = np.array(g["axis_coeff"])
        Z = _load_Z(out_dir(rid, form, cf), len(A), len(Cv))
        S, C = np.meshgrid(A, Cv)
        for c, (key, title, unit) in enumerate(ch.CHANNELS):
            ax = axes[r, c]
            _panel(ax, S, C, Z[key], (title if r == 0 else None), unit,
                   ylabel=(COEFF_LABEL[cf] if c == 0 else None), warm_y=warm.get(cf))
    fig.suptitle(f"{_label(rid)} ({rid.split('_')[0]}) · {form}: {FORM_EQ[form]} · "
                 f"square coefficient×σ slice, L$_p$={Lp:.2f} fixed, other coeffs held at fit "
                 f"(log₁₀ WSSE; rows = coefficient, cols = response channel)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.subplots_adjust(hspace=0.34, wspace=0.28)
    out = out_dir(rid, form) / "contour_grid.png"
    fig.savefig(out, dpi=140); plt.close(fig)
    print(f"[plot] {rid} {form} -> {out}", flush=True)


def _panel(ax, S, C, Zk, title, unit, *, ylabel=None, warm_y=None, fs=8.5, clab=5.5):
    import numpy as np
    ax.set_facecolor("white")
    if title:
        ax.set_title(f"{title}\n[{unit}]", fontsize=fs, fontweight="bold")
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=fs + 5, fontweight="bold")
    ax.set_xlabel("σ", fontsize=fs)
    ax.tick_params(labelsize=fs - 1.5)
    ax.set_xlim(float(S.min()), float(S.max()))      # always span the full grid
    ax.set_ylim(float(C.min()), float(C.max()))
    nfin = int(np.isfinite(Zk).sum()) if Zk is not None else 0
    if warm_y is not None:
        ax.axhline(warm_y, color="0.3", lw=0.7, ls="--", alpha=0.6)
    # need a few feasible nodes spread over ≥2 rows/cols to draw a contour
    ok2d = False
    if nfin >= 6:
        rows = np.where(np.isfinite(Zk).any(axis=1))[0]
        cols = np.where(np.isfinite(Zk).any(axis=0))[0]
        ok2d = rows.size >= 2 and cols.size >= 2 and np.nanmin(Zk) < np.nanmax(Zk)
    if not ok2d:
        msg = ("infeasible\n(B(c) leaves the\nphysical region —\nno square solve)"
               if nfin == 0 else f"too few feasible\nnodes ({nfin}) to\nform a contour")
        ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center",
                color="0.5", fontsize=fs - 0.5)
        if nfin:
            j, i = np.unravel_index(np.nanargmin(Zk), Zk.shape)
            ax.plot(S[j, i], C[j, i], "^", color="red", ms=9, mec="k", mew=0.6, zorder=6)
        return
    Zl = np.log10(np.clip(Zk, 1e-12, None))
    levels = np.linspace(np.nanmin(Zl), np.nanmax(Zl), 12)
    try:
        cs = ax.contour(S, C, Zl, levels=levels, cmap="turbo", linewidths=1.0)
        ax.clabel(cs, inline=True, fontsize=clab, fmt="%.1f")
    except Exception:
        pass
    j, i = np.unravel_index(np.nanargmin(Zk), Zk.shape)
    ax.plot(S[j, i], C[j, i], "^", color="red", ms=10, mec="k", mew=0.6, zorder=6)


def _label(rid):
    return {"MC2.05.07.24_CaCl2": "concentrating CaCl₂", "MC2.05.21.24_LaCl3": "concentrating LaCl₃",
            "MC2.05.07.24_NaCl": "concentrating NaCl"}.get(rid, rid)


def deck_composite(rid, form, channels=("mass", "icp_reten", "icp_perm"), out=None):
    """Clean coeff(rows) × selected-channels(cols) composite for slide use
    (no long suptitle; larger panels). Returns the saved PNG path."""
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    title_of = {k: t for k, t, u in ch.CHANNELS}
    unit_of = {k: u for k, t, u in ch.CHANNELS}
    coeffs = [c for c in FORMS[form]["coeffs"]
              if (out_dir(rid, form, c) / "grid.json").exists()]
    nrows, ncols = len(coeffs), len(channels)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.4 * ncols, 3.0 * nrows), squeeze=False)
    warm = None
    for r, cf in enumerate(coeffs):
        g = json.loads((out_dir(rid, form, cf) / "grid.json").read_text())
        warm = g["warm"]
        A = np.array(g["axis_sigma"]); Cv = np.array(g["axis_coeff"])
        Z = _load_Z(out_dir(rid, form, cf), len(A), len(Cv))
        S, C = np.meshgrid(A, Cv)
        for c, key in enumerate(channels):
            _panel(axes[r, c], S, C, Z[key], (title_of[key] if r == 0 else None),
                   unit_of[key], ylabel=(COEFF_LABEL[cf] if c == 0 else None),
                   warm_y=warm.get(cf), fs=12, clab=8)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.34, wspace=0.26)
    if out is None:
        out = pc._study() / "Bsigma_coeffs" / rid / f"deck_{form}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[deck] {rid} {form} -> {out}", flush=True)
    return out


A4_LANDSCAPE = (11.69, 8.27)   # inches


def render_a4(rid, forms=None, min_feasible=0, suffix=""):
    """A4-landscape multi-page PDF for one experiment: ONE PAGE per (form,
    coefficient), the coefficient × σ slice across all 7 response channels in a
    readable 2×4 layout. Re-renders from saved nodes (no solving).

    ``min_feasible`` (>0) drops any coefficient page whose grid has fewer than
    that many feasible nodes — a "clean" PDF with only the forms/coefficients
    that have a real feasible region (no degenerate all-blank pages)."""
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    forms = forms or list(FORMS)
    base = pc._study() / "Bsigma_coeffs" / rid
    out = base / f"{rid}_coeff_contours_A4{suffix}.pdf"
    pages = 0
    with PdfPages(out) as pdf:
        for form in forms:
            for cf in FORMS[form]["coeffs"]:
                gp = out_dir(rid, form, cf) / "grid.json"
                if not gp.exists():
                    continue
                g = json.loads(gp.read_text())
                A = np.array(g["axis_sigma"]); Cv = np.array(g["axis_coeff"])
                Z = _load_Z(out_dir(rid, form, cf), len(A), len(Cv))
                S, C = np.meshgrid(A, Cv)
                nfeas = int(np.isfinite(Z[ch.CH_KEYS[0]]).sum())
                if min_feasible and nfeas < min_feasible:
                    continue   # degenerate — skip the page
                fig, axes = plt.subplots(2, 4, figsize=A4_LANDSCAPE)
                axf = axes.flatten()
                for n, (key, title, unit) in enumerate(ch.CHANNELS):
                    _panel(axf[n], S, C, Z[key], f"{n + 1}. {title}", unit,
                           ylabel=(COEFF_LABEL[cf] if n % 4 == 0 else None),
                           warm_y=g["warm"].get(cf), fs=9.5, clab=6.5)
                axf[-1].axis("off")
                axf[-1].text(0.0, 0.96,
                             f"{_label(rid)}  ·  {form}\n{FORM_EQ[form]}\n\n"
                             f"Axis: {COEFF_LABEL[cf]} (y) × σ (x)\n"
                             f"L$_p$ = {g['Lp_fixed']:.2f} fixed\n"
                             f"other coeffs held at the\n{form} fit (warm start).\n\n"
                             f"▲ = per-channel optimum\n"
                             f"– – = fitted {COEFF_LABEL[cf]} value\n"
                             f"lines = log₁₀ WSSE\n\n"
                             f"Square forward-solve slice\n(no per-node optimization).\n"
                             f"Blank panels = the pinned\ncoefficients drive B(c) out of\n"
                             f"the physical region (no\nfeasible square solve there).",
                             transform=axf[-1].transAxes, va="top", ha="left", fontsize=8.5,
                             color="0.15")
                fig.suptitle(f"{_label(rid)} ({rid.split('_')[0]}) · {form}: {FORM_EQ[form]} · "
                             f"{COEFF_LABEL[cf]} × σ square slice (B as a function of interfacial conc.)",
                             fontsize=12)
                fig.tight_layout(rect=(0, 0, 1, 0.955))
                fig.subplots_adjust(hspace=0.42, wspace=0.30)
                pdf.savefig(fig); plt.close(fig); pages += 1
    print(f"[a4] {rid}: {pages} pages -> {out}", flush=True)
    return out


def main():
    if len(sys.argv) < 2:
        print(__doc__); return
    cmd = sys.argv[1]
    if cmd == "--worker":
        _worker(sys.argv[2]); return
    if cmd == "sweep":
        sweep(sys.argv[2], sys.argv[3], sys.argv[4],
              int(sys.argv[5]) if len(sys.argv) > 5 else 9); return
    if cmd == "expform":
        expform(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 9); return
    if cmd == "all":
        all_runs(int(sys.argv[2]) if len(sys.argv) > 2 else 9,
                 int(sys.argv[3]) if len(sys.argv) > 3 else 6); return
    if cmd == "plot":
        plot_expform(sys.argv[2], sys.argv[3]); return
    if cmd == "a4":
        render_a4(sys.argv[2]); return
    if cmd == "a4all":
        for rid in REPS:
            render_a4(rid)
        return
    if cmd == "allrids":
        rids = sys.argv[2].split(",")
        forms = sys.argv[5].split(",") if len(sys.argv) > 5 else None
        all_runs(int(sys.argv[3]) if len(sys.argv) > 3 else 9,
                 int(sys.argv[4]) if len(sys.argv) > 4 else 6, rids=rids, forms=forms)
        for rid in rids:
            render_a4(rid, forms=forms)
        return
    print(__doc__)


if __name__ == "__main__":
    main()
