#!/usr/bin/env python3
"""DATA3 B(c) PROFILE contours in (Lp, B, sigma) space.

The deck/DATA1 question: under a constant B, sigma-B-Lp are coupled and sigma rails
(unidentifiable).  Does letting B vary with concentration restore identifiability?
A profile contour answers it: FIX the two axis parameters at each grid node and
RE-SOLVE everything else (incl. the B(c) coefficients) -> the WSSE landscape with B
genuinely free.  (cf. the new solve_model_B_fix(fix_vars=...) hook.)

Two hard problems, both handled here:
  1. Cold DAE solves are slow/divergent off the optimum.  FIX: warm-start
     continuation -- march the grid outward from the campaign optimum, seeding each
     node's optimizer from an already-converged neighbour (per the user's directive:
     hold the axes, sweep, solve the third for a good initial guess).  Feasible nodes
     then converge in ~3-8 s instead of grinding.
  2. The casadi/idas Simulator init is UNCAPPED (solver_max_cpu_time only caps IPOPT);
     intrinsically stiff (high-sigma) nodes grind for minutes regardless of seed.
     FIX: the worker marches in-process (library imported once) writing each node to a
     JSONL immediately; the PARENT monitors progress and, if a node stalls past
     PER_NODE_TIMEOUT, kills the worker process group (orphaned ipopt/idas included,
     cf. _recover_donnan), records the stuck node as a hole, and restarts the worker,
     which resumes from the JSONL.  Resumable + bounded.

Planes (axes -> profiled-out free params):
  sigLp : (Lp, sigma) -> B(c) coeffs (+ S)        [ALL B-forms; the key panel]
  BLp   : (Lp, B)     -> sigma (+ S)              [constant-B only -- B is scalar]
  sigB  : (sigma, B)  -> Lp (+ S)                 [constant-B only]

Output (bform_study/profile_contours/<sheet>/<form>/<plane>/):
  nodes.jsonl       one line per grid node (resumable)
  grid.json         axis vectors + anchor + meta
  contour.png       legacy plot_contour style: per-channel (mass|perm|reten) + combined

Usage:
  python3 _run_profile_contours.py plan                      # list the campaign + node counts
  python3 _run_profile_contours.py run <rid> <form> <plane> [density] [pad]
  python3 _run_profile_contours.py --worker <gridpath>       # internal
  python3 _run_profile_contours.py plot <rid> <form> <plane> # (re)render from nodes.jsonl
"""
import os, sys, json, time, signal, subprocess, math
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

PER_NODE_TIMEOUT = 45      # wall seconds before the parent declares a node truly hung and kills the worker
IPOPT_CPU = 10             # IPOPT max_cpu_time per node; with skip_sim_init a hard node fails fast under this
NFE = 100
HEADLINE_4 = ["MC3.07.22.24_SNaCl", "MC2.05.07.24_NaCl", "MC2.05.07.24_CaCl2", "MC2.05.21.24_LaCl3"]
# salt-regime groups for STACKED objective-landscape pooling (experiments informing each other)
POOL_GROUPS = {
    "NaCl_diluting":      ["MC3.07.22.24_SNaCl", "MC5.07.23.24_S2NaCl"],
    "NaCl_concentrating": ["MC4.07.11.24_SNaCl", "MC2.05.07.24_NaCl", "MC5.07.23.24_NaCl", "MC5.07.23.24_SNaCl"],
    "CaCl2":              ["MC2.05.07.24_CaCl2", "MC3.07.11.24_SCaCl2", "MC3.07.12.24_S2CaCl2"],
    "LaCl3":              ["MC2.05.21.24_LaCl3", "MC4.07.11.24_SLaCl3"],
}

# form token -> (B_form arg for solve_model_B_fix, result json tag, n free B-params)
FORMS = {"single": ("single", "single"), "poly1": (1, "poly1"), "poly2": (2, "poly2"),
         "poly3": (3, "poly3"), "sat": ("sat", "sat"), "donnan": ("donnan", "donnan")}
# which planes make sense per form (B-axis planes need a scalar B -> constant-B only)
PLANES_FOR = {"single": ["sigLp", "BLp", "sigB"], "poly1": ["sigLp"], "poly2": ["sigLp"],
              "poly3": ["sigLp"], "sat": ["sigLp"], "donnan": ["sigLp"]}

SALT_SHEETS = {
    "NaCl": ["MC3.07.22.24_SNaCl", "MC5.07.23.24_S2NaCl", "MC4.07.11.24_SNaCl",
             "MC2.05.07.24_NaCl", "MC5.07.23.24_NaCl", "MC5.07.23.24_SNaCl"],
    "CaCl2": ["MC2.05.07.24_CaCl2", "MC3.07.11.24_SCaCl2", "MC3.07.12.24_S2CaCl2"],
    "LaCl3": ["MC2.05.21.24_LaCl3", "MC4.07.11.24_SLaCl3"],
}


def _study():
    import _run_bform_campaign as camp
    return camp.OUT


def out_dir(rid, form, plane):
    return _study() / "profile_contours" / rid / form / plane


# ----------------------------------------------------------------------------- worker
def _worker(gridpath):
    """March the grid in warm-start order, writing each node to nodes.jsonl.
    Resumes from whatever is already written; skips nodes listed in nodes.jsonl.skip."""
    import contextlib, io
    import numpy as np
    import _run_bform_campaign as camp
    import refactored_ucb_library as lib
    lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0

    grid = json.loads(Path(gridpath).read_text())
    rid, form, plane = grid["rid"], grid["form"], grid["plane"]
    A, Bv = np.array(grid["axisA"]), np.array(grid["axisB"])   # axisA=x, axisB=y vectors
    anchor = grid["anchor"]; ai, aj = grid["anchor_ij"]
    od = Path(grid["out_dir"]); od.mkdir(parents=True, exist_ok=True)
    jsonl = od / "nodes.jsonl"; skipf = od / "nodes.jsonl.skip"; curf = od / "nodes.jsonl.cur"
    Bform = FORMS[form][0]
    ds = camp.load(rid); mode = ds["mode"]

    # axis -> fixed variable names
    xname, yname = grid["xname"], grid["yname"]

    done, solved = set(), {}
    if jsonl.exists():
        for line in jsonl.read_text().splitlines():
            try:
                r = json.loads(line)
                done.add((r["i"], r["j"]))
                if r["status"] == "ok":
                    solved[(r["i"], r["j"])] = r["params"]
            except Exception:
                pass
    skip = set()
    if skipf.exists():
        for line in skipf.read_text().splitlines():
            try:
                a, b = line.split(","); skip.add((int(a), int(b)))
            except Exception:
                pass
    if (ai, aj) not in solved:
        solved[(ai, aj)] = dict(anchor)

    nA, nB = len(A), len(B := Bv)
    # march order: by grid distance from the anchor (every node has a nearer solved neighbour)
    order = sorted([(i, j) for i in range(nA) for j in range(nB)],
                   key=lambda ij: (ij[0] - ai) ** 2 + (ij[1] - aj) ** 2)

    def nearest_seed(i, j):
        best, bd = None, 1e18
        for (pi, pj), th in solved.items():
            d = (pi - i) ** 2 + (pj - j) ** 2
            if d < bd:
                bd, best = d, th
        return best or anchor

    method = grid.get("method", "profile")
    fh = open(jsonl, "a")
    for (i, j) in order:
        if (i, j) in done or (i, j) in skip:
            continue
        curf.write_text(f"{i},{j}")          # tell the parent which node is in flight
        xval, yval = float(A[i]), float(B[j])
        rec = {"i": i, "j": j, "x": xval, "y": yval}
        t = time.time()
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                if method == "slice":
                    # legacy calc_contour_2d style: HOLD B(c) at the fitted curve, set the two
                    # axis params, forward-simulate the transport model, read the objective.
                    seed = dict(anchor); seed[xname] = xval; seed[yname] = yval
                    fit, _, _ = lib.solve_model(ds, mode, theta=seed, sim_opt=True, B_form=Bform,
                                                workflow_family="DATA3", nfe=120, solver_max_cpu_time=30)
                else:
                    # profile: fix the two axis params, RE-SOLVE the others (B(c) coeffs + S)
                    seed = dict(nearest_seed(i, j)); seed[xname] = xval; seed[yname] = yval
                    fit, _, _ = lib.solve_model_B_fix(ds, mode, theta=seed, sim_opt=False, B_form=Bform,
                                                      workflow_family="DATA3", nfe=NFE,
                                                      solver_max_cpu_time=IPOPT_CPU,
                                                      fix_vars={xname: xval, yname: yval},
                                                      skip_sim_init=grid.get("skip_sim_init", True))
            if isinstance(fit, dict):
                rec.update(status="ok", obj_m=float(fit.get("obj_m")), obj_cv=float(fit.get("obj_cv")),
                           obj_cr=float(fit.get("obj_cr")),
                           WSSE3=float(sum(fit.get(k, 0.0) for k in ("obj_m", "obj_cv", "obj_cr"))),
                           params=fit["parameters"], dt=round(time.time() - t, 1))
                if method != "slice":
                    solved[(i, j)] = fit["parameters"]
            else:
                rec.update(status="infeasible", dt=round(time.time() - t, 1))
        except Exception as e:
            rec.update(status="error", err=str(e)[:80], dt=round(time.time() - t, 1))
        fh.write(json.dumps(rec) + "\n"); fh.flush(); os.fsync(fh.fileno())
    fh.close()
    if curf.exists():
        curf.unlink()
    (od / "WORKER_DONE").write_text("done")


# ----------------------------------------------------------------------------- parent
def _launch_worker(gridpath):
    return subprocess.Popen([sys.executable, str(HERE / "_run_profile_contours.py"), "--worker", str(gridpath)],
                            start_new_session=True, cwd=str(HERE))


def _build_grid(rid, form, plane, density, pad):
    import numpy as np
    study = _study()
    tag = FORMS[form][1]
    rp = study / rid / f"result_{tag}.json"
    base = json.loads(rp.read_text())
    if "parameters" not in base:
        raise SystemExit(f"no converged {form} fit for {rid} (need it as the anchor)")
    p = base["parameters"]
    Lp0 = float(p["Lp"]); sg0 = float(p.get("sigma", 1.0))
    B0 = float(p.get("B", p.get("beta_0", 1.0)))
    # axis ranges
    Lp_axis = np.linspace(max(0.2, Lp0 * (1 - pad)), Lp0 * (1 + pad), density)
    sg_axis = np.linspace(0.0, 1.0, density)
    B_axis = np.linspace(max(1e-3, B0 * (1 - pad)), B0 * (1 + pad), density)
    if plane == "sigLp":
        xname, yname, A, Bv = "Lp", "sigma", Lp_axis, sg_axis
    elif plane == "BLp":
        xname, yname, A, Bv = "Lp", "B", Lp_axis, B_axis
    elif plane == "sigB":
        xname, yname, A, Bv = "sigma", "B", sg_axis, B_axis
    else:
        raise SystemExit(f"unknown plane {plane}")
    ai = int(np.argmin(abs(A - p[xname]))) if xname in p else int(np.argmin(abs(A - (Lp0 if xname == "Lp" else sg0))))
    aj = int(np.argmin(abs(Bv - p[yname]))) if yname in p else int(np.argmin(abs(Bv - (sg0 if yname == "sigma" else B0))))
    od = out_dir(rid, form, plane)
    od.mkdir(parents=True, exist_ok=True)
    grid = {"rid": rid, "form": form, "plane": plane, "xname": xname, "yname": yname,
            "axisA": A.tolist(), "axisB": Bv.tolist(), "anchor": p, "anchor_ij": [ai, aj],
            "out_dir": str(od), "skip_sim_init": True}
    (od / "grid.json").write_text(json.dumps(grid))
    return od, len(A) * len(Bv)


def _run(rid, form, plane, density, pad):
    od, ntot = _build_grid(rid, form, plane, density, pad)
    gridpath = od / "grid.json"
    jsonl = od / "nodes.jsonl"; skipf = od / "nodes.jsonl.skip"; curf = od / "nodes.jsonl.cur"
    donef = od / "WORKER_DONE"
    if donef.exists():
        donef.unlink()
    print(f"[contour] {rid} {form} {plane}  {ntot} nodes  (timeout {PER_NODE_TIMEOUT}s/node)", flush=True)
    t0 = time.time()
    while True:
        if donef.exists():
            break
        n_done = sum(1 for _ in jsonl.open()) if jsonl.exists() else 0
        n_skip = sum(1 for _ in skipf.open()) if skipf.exists() else 0
        if n_done + n_skip >= ntot:
            break
        proc = _launch_worker(gridpath)
        # monitor progress: kill+restart if no new node within PER_NODE_TIMEOUT
        last_count = n_done
        last_progress = time.time()
        while proc.poll() is None:
            time.sleep(2)
            cur = sum(1 for _ in jsonl.open()) if jsonl.exists() else 0
            if cur > last_count:
                last_count = cur; last_progress = time.time()
            elif time.time() - last_progress > PER_NODE_TIMEOUT:
                # stuck node -> record as skip, kill worker, restart
                stuck = curf.read_text().strip() if curf.exists() else None
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except Exception:
                    pass
                proc.wait()
                if stuck:
                    with open(skipf, "a") as sf:
                        sf.write(stuck + "\n")
                    print(f"  [kill] stuck node {stuck} -> skip ({last_count}/{ntot} done)", flush=True)
                break
        else:
            proc.wait()
        if donef.exists():
            break
    n_done = sum(1 for _ in jsonl.open()) if jsonl.exists() else 0
    n_skip = sum(1 for _ in skipf.open()) if skipf.exists() else 0
    print(f"[contour] {rid} {form} {plane} complete: {n_done} solved, {n_skip} holes, "
          f"{ntot} total in {time.time()-t0:.0f}s", flush=True)
    _plot(rid, form, plane)


# ----------------------------------------------------------------------------- plot
def _plot(rid, form, plane):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    od = out_dir(rid, form, plane)
    grid = json.loads((od / "grid.json").read_text())
    A, Bv = np.array(grid["axisA"]), np.array(grid["axisB"])
    nA, nB = len(A), len(Bv)
    Z = {k: np.full((nB, nA), np.nan) for k in ("obj_m", "obj_cv", "obj_cr", "WSSE3")}
    jsonl = od / "nodes.jsonl"
    if jsonl.exists():
        for line in jsonl.read_text().splitlines():
            try:
                r = json.loads(line)
                if r["status"] == "ok":
                    for k in Z:
                        Z[k][r["j"], r["i"]] = r[k]
            except Exception:
                pass
    X, Y = np.meshgrid(A, Bv)
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))
    panels = [("obj_m", "Mass Residual$^2$ [g$^2$]"), ("obj_cv", "Permeate-conc Residual$^2$ [mM$^2$]"),
              ("obj_cr", "Retentate-conc Residual$^2$ [mM$^2$]"), ("WSSE3", "Combined WSSE$_3$")]
    for ax, (key, title) in zip(axes, panels):
        _contour_panel(ax, X, Y, Z[key], A, Bv, title, _LAB.get(grid["xname"], grid["xname"]),
                       _LAB.get(grid["yname"], grid["yname"]))
    fig.suptitle(f"{rid}  ·  B-form: {form}  ·  profile contour ({grid['xname']}×{grid['yname']}, "
                 f"B(c) re-solved at each node)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(od / "contour.png", dpi=150)
    plt.close(fig)
    print(f"  plotted -> {od/'contour.png'}", flush=True)


# ----------------------------------------------------------------------------- cli
def _plan():
    study = _study()
    total = 0
    for salt, sheets in SALT_SHEETS.items():
        for rid in sheets:
            for form in FORMS:
                tag = FORMS[form][1]
                if not (study / rid / f"result_{tag}.json").exists():
                    continue
                r = json.loads((study / rid / f"result_{tag}.json").read_text())
                if "parameters" not in r:
                    continue
                for plane in PLANES_FOR[form]:
                    total += 1
                    print(f"  {salt:6s} {rid:22s} {form:7s} {plane}")
    print(f"\n{total} contours to compute (each {{21x21=441}} nodes at deck density).")


PLANE_AXES = {"sigLp": ("Lp", "sigma"), "BLp": ("Lp", "B"), "sigB": ("sigma", "B")}
_LAB = {"Lp": "L$_p$ [L m$^{-2}$ h$^{-1}$ bar$^{-1}$]", "sigma": "$\\sigma$ [-]", "B": "B [$\\mu$m s$^{-1}$]"}


def _contour_panel(ax, X, Y, Zk, A, Bv, title, xlab, ylab):
    """Labeled line-contour panel (DATA1 plot_contour.m style): colored level lines with inline
    labels on a white background, red ▲ at the minimum + a min annotation.  log10 of the residual."""
    import numpy as np
    ax.set_facecolor("white"); ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlab, fontsize=9); ax.set_ylabel(ylab, fontsize=9)
    if Zk is None or np.all(np.isnan(Zk)):
        ax.text(0.5, 0.5, "(no data)", transform=ax.transAxes, ha="center", color="0.5"); return
    Zl = np.log10(np.clip(Zk, 1e-6, None))
    try:
        cs = ax.contour(X, Y, Zl, levels=12, cmap="turbo", linewidths=1.4)
        ax.clabel(cs, inline=True, fontsize=6, fmt="%.2f")
    except Exception:
        pass
    ax.grid(alpha=0.2, lw=0.4)
    jm, im = np.unravel_index(np.nanargmin(Zk), Zk.shape)
    ax.plot(A[im], Bv[jm], "^", ms=11, color="red", mec="white", mew=1.0, zorder=6)
    ax.text(0.5, -0.24, f"min ({A[im]:.2f}, {Bv[jm]:.2f})  log₁₀={Zl[jm, im]:.2f}",
            transform=ax.transAxes, ha="center", fontsize=7, color="red")


def _common_axis(name, ranges, density):
    import numpy as np
    if name == "sigma":
        return list(np.linspace(0.0, 1.0, density))
    lo, hi = min(ranges[name]), max(ranges[name])
    return list(np.linspace(max(0.05, lo * 0.55), hi * 1.45, density))


def _grid_common(rid, form, plane, xname, yname, anchor, Xax, Yax, od, method="profile"):
    """Write a grid.json for one experiment on a COMMON xname x yname grid (for stacking)."""
    import numpy as np
    od.mkdir(parents=True, exist_ok=True)
    ai = int(np.argmin(abs(np.array(Xax) - anchor.get(xname, 5.0))))
    aj = int(np.argmin(abs(np.array(Yax) - anchor.get(yname, 1.0))))
    grid = {"rid": rid, "form": form, "plane": plane, "xname": xname, "yname": yname,
            "axisA": list(Xax), "axisB": list(Yax), "anchor": anchor, "anchor_ij": [ai, aj],
            "out_dir": str(od), "skip_sim_init": True, "method": method}
    (od / "grid.json").write_text(json.dumps(grid))
    return len(Xax) * len(Yax)


def _done(od, ntot):
    j = od / "nodes.jsonl"; s = od / "nodes.jsonl.skip"
    nd = sum(1 for _ in j.open()) if j.exists() else 0
    ns = sum(1 for _ in s.open()) if s.exists() else 0
    return (od / "WORKER_DONE").exists() or nd + ns >= ntot


def _drive_parallel(jobs, maxpar=4):
    """Run a set of common-grid workers concurrently; kill+restart any that stall a node."""
    pending = [j for j in jobs if not _done(j[0], j[1])]
    active = {}
    while pending or active:
        while pending and len(active) < maxpar:
            od, ntot, name = pending.pop(0)
            if _done(od, ntot):
                continue
            if (od / "WORKER_DONE").exists():
                (od / "WORKER_DONE").unlink()
            p = _launch_worker(od / "grid.json")
            active[od] = {"p": p, "ntot": ntot, "name": name, "last": 0, "lastp": time.time()}
            print(f"  [stack] launch {name} ({len(active)} running, {len(pending)} queued)", flush=True)
        time.sleep(2)
        for od, st in list(active.items()):
            p = st["p"]; j = od / "nodes.jsonl"
            cur = sum(1 for _ in j.open()) if j.exists() else 0
            if p.poll() is not None:
                if _done(od, st["ntot"]):
                    print(f"  [stack] done {st['name']}", flush=True)
                else:
                    pending.append((od, st["ntot"], st["name"]))
                del active[od]
            elif cur > st["last"]:
                st["last"] = cur; st["lastp"] = time.time()
            elif time.time() - st["lastp"] > PER_NODE_TIMEOUT:
                stuck = (od / "nodes.jsonl.cur").read_text().strip() if (od / "nodes.jsonl.cur").exists() else None
                try:
                    os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                except Exception:
                    pass
                p.wait()
                if stuck:
                    with open(od / "nodes.jsonl.skip", "a") as sf:
                        sf.write(stuck + "\n")
                pending.append((od, st["ntot"], st["name"]))
                del active[od]


def _load_grid(od, key, nA, nB):
    import numpy as np
    Z = np.full((nB, nA), np.nan)
    j = od / "nodes.jsonl"
    if j.exists():
        for line in j.read_text().splitlines():
            try:
                r = json.loads(line)
                if r["status"] == "ok":
                    Z[r["j"], r["i"]] = r[key]
            except Exception:
                pass
    return Z


def stack(group, form="single", density=11, planes=None, method="profile"):
    """STACK per-experiment objective landscapes of a group on COMMON grids, for EACH
    parameter plane (default all three: σ×Lp, B×Lp, σ×B).  The summed surface is the pooled
    objective (= sum of Fisher information); its minimum is the joint-fit initial guess, and
    its sharpness vs the individuals shows the experiments informing each other.  B-axis
    planes need a scalar B -> form must be 'single'.  method='profile' re-solves the other
    params at each node (identifiability); method='slice' holds B(c) at the fit and forward-
    simulates (legacy calc_contour_2d, cheaper sensitivity)."""
    import json as _json
    study = _study(); tag = FORMS[form][1]
    planes = planes or (["sigLp", "BLp", "sigB"] if form == "single" else ["sigLp"])
    anchors, ranges = {}, {"Lp": [], "sigma": [], "B": []}
    for s in POOL_GROUPS[group]:
        rp = study / s / f"result_{tag}.json"
        if not rp.exists():
            continue
        r = _json.loads(rp.read_text())
        if "parameters" in r:
            p = r["parameters"]; anchors[s] = p
            ranges["Lp"].append(p.get("Lp", 5.0)); ranges["sigma"].append(p.get("sigma", 1.0))
            ranges["B"].append(p.get("B", p.get("beta_0", 1.0)))
    if not anchors:
        print(f"[stack] {group}: no {form} fits"); return
    for plane in planes:
        xname, yname = PLANE_AXES[plane]
        Xax = _common_axis(xname, ranges, density); Yax = _common_axis(yname, ranges, density)
        root = study / "stack_contours" / group / form / plane / method
        jobs = []
        for s in anchors:
            od = root / s
            ntot = _grid_common(s, form, plane, xname, yname, anchors[s], Xax, Yax, od, method)
            jobs.append((od, ntot, s.split("_")[0]))
        print(f"[stack] {group} ({form}) {plane} [{method}]: {len(jobs)} exp on common {density}x{density} {xname}×{yname}", flush=True)
        _drive_parallel(jobs, maxpar=min(4, len(jobs)))
        _stackplot(group, form, plane, method, xname, yname, Xax, Yax, list(anchors))


def _stackplot(group, form, plane, method, xname, yname, Xax, Yax, sheets):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    study = _study(); root = study / "stack_contours" / group / form / plane / method
    A, Bv = np.array(Xax), np.array(Yax); nA, nB = len(A), len(Bv)
    X, Y = np.meshgrid(A, Bv)
    keys = ["obj_m", "obj_cv", "obj_cr", "WSSE3"]
    stackZ = {k: np.zeros((nB, nA)) for k in keys}
    cov = np.zeros((nB, nA)); per_exp_min = []
    for s in sheets:
        od = root / s
        w = _load_grid(od, "WSSE3", nA, nB)
        if np.all(np.isnan(w)):
            continue
        jm, im = np.unravel_index(np.nanargmin(w), w.shape)
        per_exp_min.append((float(A[im]), float(Bv[jm]), s.split("_")[0]))
        for k in keys:
            Zk = _load_grid(od, k, nA, nB); m = ~np.isnan(Zk); stackZ[k][m] += Zk[m]
        cov += (~np.isnan(w)).astype(float)
    full = cov >= len(sheets)
    for k in keys:
        stackZ[k][~full] = np.nan
    comb = stackZ["WSSE3"]; star = None
    if np.any(~np.isnan(comb)):
        jm, im = np.unravel_index(np.nanargmin(comb), comb.shape)
        star = (float(A[im]), float(Bv[jm]), float(comb[jm, im]))
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))
    titles = {"obj_m": "Σ mass", "obj_cv": "Σ permeate-conc", "obj_cr": "Σ retentate-conc", "WSSE3": "Σ combined"}
    for ax, k in zip(axes, keys):
        Zk = stackZ[k]
        ax.set_xlabel(_LAB[xname]); ax.set_ylabel(_LAB[yname])
        if np.all(np.isnan(Zk)):
            ax.set_title(f"{titles[k]}\n(no overlap)"); continue
        Zl = np.log10(np.clip(Zk, 1e-6, None))
        cs = ax.contourf(X, Y, Zl, levels=18, cmap="viridis")
        ax.contour(X, Y, Zl, levels=18, colors="k", linewidths=0.3, alpha=0.4)
        for (xx, yy, nm) in per_exp_min:
            ax.plot(xx, yy, "o", ms=6, mfc="white", mec="0.2")
        if star:
            ax.plot(star[0], star[1], "r*", ms=18, mec="white")
        fig.colorbar(cs, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(titles[k], fontsize=10)
    sub = f"stacked min ★ {xname}={star[0]:.2f}, {yname}={star[1]:.2f}" if star else "no full-overlap region"
    how = "other params re-solved" if method == "profile" else "B(c) held at fit, forward-sim"
    fig.suptitle(f"STACKED {xname}×{yname} — {group} · B-form={form} · {method} ({how}): {len(sheets)} experiments summed.  "
                 f"○ = each experiment's own min (scattered = sloppy); ★ = pooled min / joint-fit seed.   {sub}", fontsize=9.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(root / "stacked.png", dpi=150); plt.close(fig)
    (root / "stacked.json").write_text(json.dumps(
        {"group": group, "form": form, "plane": plane, "method": method,
         "stacked_min": ({xname: star[0], yname: star[1], "WSSE3": star[2]} if star else None),
         "per_exp_min": [{xname: x, yname: y, "exp": n} for (x, y, n) in per_exp_min]}, indent=2))
    print(f"[stack] {group} {form} {plane} [{method}]: min {xname}={star[0]:.2f} {yname}={star[1]:.2f} -> {root/'stacked.png'}"
          if star else f"[stack] {group} {form} {plane} [{method}]: no overlap", flush=True)


def _contours_for(scope):
    """Enumerate (rid, form, plane) contours for a scope ('deck' or 'full')."""
    study = _study()
    sheets = HEADLINE_4 if scope == "deck" else [s for ss in SALT_SHEETS.values() for s in ss]
    out = []
    for rid in sheets:
        for form in FORMS:
            tag = FORMS[form][1]
            rp = study / rid / f"result_{tag}.json"
            if not rp.exists():
                continue
            try:
                if "parameters" not in json.loads(rp.read_text()):
                    continue
            except Exception:
                continue
            for plane in PLANES_FOR[form]:
                out.append((rid, form, plane))
    return out


def _campaign(scope, density, maxpar):
    contours = _contours_for(scope)
    # skip contours already complete
    todo = [c for c in contours if not (out_dir(*c) / "WORKER_DONE").exists()]
    print(f"[campaign] scope={scope} density={density}: {len(contours)} contours, "
          f"{len(todo)} to run, {maxpar}-way parallel", flush=True)
    running = {}   # proc -> (rid,form,plane)
    queue = list(todo)
    while queue or running:
        while queue and len(running) < maxpar:
            rid, form, plane = queue.pop(0)
            log = out_dir(rid, form, plane) / "run.log"
            log.parent.mkdir(parents=True, exist_ok=True)
            fh = open(log, "w")
            proc = subprocess.Popen(
                [sys.executable, str(HERE / "_run_profile_contours.py"), "run", rid, form, plane,
                 str(density), "0.45"], stdout=fh, stderr=subprocess.STDOUT, cwd=str(HERE))
            running[proc] = (rid, form, plane, fh)
            print(f"  launch {rid} {form} {plane}  ({len(running)} running, {len(queue)} queued)", flush=True)
        done = [p for p in running if p.poll() is not None]
        for p in done:
            rid, form, plane, fh = running.pop(p)
            fh.close()
            print(f"  done   {rid} {form} {plane}  rc={p.returncode}", flush=True)
        if not done:
            time.sleep(5)
    print(f"[campaign] scope={scope} COMPLETE", flush=True)


def main():
    if len(sys.argv) < 2:
        print(__doc__); return
    cmd = sys.argv[1]
    if cmd == "--worker":
        _worker(sys.argv[2]); return
    if cmd == "plan":
        _plan(); return
    if cmd == "campaign":
        scope = sys.argv[2] if len(sys.argv) > 2 else "deck"
        density = int(sys.argv[3]) if len(sys.argv) > 3 else 15
        maxpar = int(sys.argv[4]) if len(sys.argv) > 4 else 6
        _campaign(scope, density, maxpar); return
    if cmd == "plot":
        _plot(sys.argv[2], sys.argv[3], sys.argv[4]); return
    if cmd == "stack":
        group = sys.argv[2]; form = sys.argv[3] if len(sys.argv) > 3 else "single"
        density = int(sys.argv[4]) if len(sys.argv) > 4 else 11
        method = sys.argv[5] if len(sys.argv) > 5 else "profile"
        planes = None if form == "single" else ["sigLp"]
        stack(group, form, density, planes, method); return
    if cmd == "bformsweep":
        # per-B-form sigLp stacking for a group: all 6 forms for one method
        grp = sys.argv[2]; density = int(sys.argv[3]) if len(sys.argv) > 3 else 11
        method = sys.argv[4] if len(sys.argv) > 4 else "slice"
        for form in ("single", "poly1", "poly2", "poly3", "sat", "donnan"):
            try:
                stack(grp, form, density, ["sigLp"], method)
            except Exception as e:
                print(f"[bformsweep] {grp} {form} {method} FAILED: {str(e)[:80]}", flush=True)
        return
    if cmd == "stackplot":
        grp = sys.argv[2]; form = sys.argv[3] if len(sys.argv) > 3 else "single"
        plane = sys.argv[4] if len(sys.argv) > 4 else "sigLp"
        method = sys.argv[5] if len(sys.argv) > 5 else "profile"
        xname, yname = PLANE_AXES[plane]
        root = _study() / "stack_contours" / grp / form / plane / method
        sheets = [s for s in POOL_GROUPS[grp] if (root / s / "grid.json").exists()]
        g0 = json.loads((root / sheets[0] / "grid.json").read_text())
        _stackplot(grp, form, plane, method, xname, yname, g0["axisA"], g0["axisB"], sheets); return
    if cmd == "run":
        rid, form, plane = sys.argv[2], sys.argv[3], sys.argv[4]
        density = int(sys.argv[5]) if len(sys.argv) > 5 else 21
        pad = float(sys.argv[6]) if len(sys.argv) > 6 else 0.4
        _run(rid, form, plane, density, pad); return
    print(__doc__)


if __name__ == "__main__":
    main()
