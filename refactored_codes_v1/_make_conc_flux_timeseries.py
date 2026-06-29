#!/usr/bin/env python3
"""DATA3 §6 — concentration & flux time-series (the DATA2-Fig-7 analogue), per sheet.

Forward-simulate a constant-B fit (default result_single.json) at the canonical 2% no-floor
weighting and plot, concatenated over the diafiltration vials vs time:
   A: interface feed conc c_in,f (= cIn, conc-polarized wall) and permeate-side c_h
   B: water flux J_w(t)
   C: solute flux J_s(t)

Also supports a SECOND theta overlay (for §7b/§8: single-sheet fit vs pooled estimate).
Output: bform_study/conc_flux/<run_id>/conc_flux.png  (+ .csv)

Usage:
  python3 _make_conc_flux_timeseries.py <run_id>            # single fit
  python3 _make_conc_flux_timeseries.py all                 # all 11 single-salt reps/sheets
"""
import sys, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

REPS = {"MC3.07.22.24_SNaCl": "NaCl — diluting", "MC2.05.07.24_NaCl": "NaCl — concentrating",
        "MC2.05.07.24_CaCl2": "CaCl₂ — concentrating", "MC2.05.21.24_LaCl3": "LaCl₃ — concentrating"}


def _study():
    import _run_bform_campaign as camp
    return camp.OUT


def _sim_series(rid, theta, form="single"):
    """Forward-sim at theta; return concatenated (t_min, cIn, cH, Jw, Js) over vials."""
    import _run_bform_campaign as camp
    import refactored_ucb_library as lib
    lib.NF270_CF_RESIDUAL_SCALE_FRACTION = 0.02      # canonical 2% no-floor
    lib.NF270_CF_RESIDUAL_FLOOR_MM = None
    ds = camp.load(rid); mode = ds["mode"]
    Bform = "single"
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        fit, sim, _ = lib.solve_model(ds, mode, theta=dict(theta), sim_opt=True, B_form=Bform,
                                      workflow_family="DATA3", nfe=160, solver_max_cpu_time=60)
    if sim is None:
        return None
    t, cIn, cH, Jw, Js = [], [], [], [], []
    for i in sorted(sim.keys()):
        s = sim[i]
        t.extend(list(s["time"])); cIn.extend(list(s["cIn"])); cH.extend(list(s["cH"]))
        Jw.extend(list(s["Jw"])); Js.extend(list(s.get("Js", [np.nan] * len(s["time"]))))
    t = np.asarray(t, float); order = np.argsort(t)
    t = t[order]; t = (t - t[0]) / 60.0   # minutes from start
    # Jw is stored in cm/s (mass balance: Am·ρ·Jw, lib line ~1955) → ×1e4 for µm/s.
    # Js is already µmol·cm⁻²·s⁻¹ (dcF eqn: cF·Jw − Js, lib line ~2017) → no scaling.
    return {"t": t, "cIn": np.asarray(cIn, float)[order], "cH": np.asarray(cH, float)[order],
            "Jw": np.asarray(Jw, float)[order] * 1e4, "Js": np.asarray(Js, float)[order]}


def _meas_jw(rid):
    """Measured water flux per vial from the permeate mass accumulation:
    Jw = (dm/dt)/(ρ·Am)·1e4  [µm s⁻¹], with dm/dt the per-vial slope of collected
    mass (g) vs time (s).  Returns (t_min, Jw_meas) — one point per vial, on the
    same 'minutes from start' axis as the simulated trace."""
    import _run_bform_campaign as camp
    ds = camp.load(rid)
    cfg = ds.get("data_config", {})
    Am = float(cfg.get("Am", 4.1)); rho = float(cfg.get("rho", 1.0))
    # Drop the startup vial collected BEFORE steady permeate flow begins
    # (t_perm_start_s): it accumulates ~no mass, so its slope reads ~0 and would
    # drag the measured-Jw line to zero at the first point.
    t_perm_start = cfg.get("t_perm_start_s")
    rows = ds.get("data_raw", [])
    t0 = None
    for r in rows:
        tt = np.asarray(r.get("time", []), float).reshape(-1)
        if tt.size:
            t0 = tt[0]; break
    tm, jw = [], []
    for r in rows:
        t = np.asarray(r.get("time", []), float).reshape(-1)
        m = np.asarray(r.get("mass", []), float).reshape(-1)
        n = min(t.size, m.size)
        if n < 3 or t0 is None:
            continue
        if t_perm_start is not None and float(t[:n].mean()) < float(t_perm_start):
            continue                              # skip the pre-permeate startup vial
        slope = np.polyfit(t[:n], m[:n], 1)[0]   # g/s
        if not np.isfinite(slope) or slope <= 0:  # guard NaN (nan<=0 is False)
            continue
        tm.append((float(t[:n].mean()) - t0) / 60.0)
        jw.append(slope / (rho * Am) * 1e4)
    return np.asarray(tm, float), np.asarray(jw, float)


def make(rid, theta=None, theta2=None, label1="fit", label2=None, title=None):
    study = _study()
    if theta is None:
        theta = json.loads((study / rid / "result_single.json").read_text())["parameters"]
    s1 = _sim_series(rid, theta)
    if s1 is None:
        print(f"[conc_flux] {rid}: sim failed"); return None
    s2 = _sim_series(rid, theta2) if theta2 is not None else None
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    # A: cIn,f and cH
    ax[0].plot(s1["t"], s1["cIn"], color="tab:green", lw=2, label=f"c$_{{in,f}}$ ({label1})")
    ax[0].plot(s1["t"], s1["cH"], color="tab:red", lw=2, label=f"c$_h$ ({label1})")
    if s2 is not None:
        ax[0].plot(s2["t"], s2["cIn"], color="tab:green", lw=2, ls="--", label=f"c$_{{in,f}}$ ({label2})")
        ax[0].plot(s2["t"], s2["cH"], color="tab:red", lw=2, ls="--", label=f"c$_h$ ({label2})")
    ax[0].set_xlabel("Time [min]"); ax[0].set_ylabel("Concentration [mM]"); ax[0].set_title("A · interface concentrations"); ax[0].legend(fontsize=8)
    # B: Jw — model line + measured flux (from permeate-mass slope) as a dotted overlay
    ax[1].plot(s1["t"], s1["Jw"], color="tab:blue", lw=2, label=f"J$_w$ ({label1})")
    if s2 is not None:
        ax[1].plot(s2["t"], s2["Jw"], color="tab:blue", lw=2, ls="--", label=f"J$_w$ ({label2})")
    _tm, _jw = _meas_jw(rid)
    if _tm.size:
        ax[1].plot(_tm, _jw, color="k", ls=":", lw=1.6, marker="o", ms=4,
                   label="J$_w$ (measured)")
    ax[1].set_xlabel("Time [min]"); ax[1].set_ylabel("J$_w$ [µm s$^{-1}$]"); ax[1].set_title("B · water flux"); ax[1].legend(fontsize=8)
    # C: Js
    ax[2].plot(s1["t"], s1["Js"], color="black", lw=2, label=f"J$_s$ ({label1})")
    if s2 is not None:
        ax[2].plot(s2["t"], s2["Js"], color="tab:orange", lw=2, ls="--", label=f"J$_s$ ({label2})")
    ax[2].set_xlabel("Time [min]"); ax[2].set_ylabel("J$_s$ [µmol cm$^{-2}$ s$^{-1}$]"); ax[2].set_title("C · solute flux"); ax[2].legend(fontsize=8)
    fig.suptitle(title or f"{REPS.get(rid, rid)}  ·  {rid}  ·  concentration & flux time-series", fontsize=12)
    for a in ax: a.grid(alpha=0.25, lw=0.4)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    od = study / "conc_flux" / rid; od.mkdir(parents=True, exist_ok=True)
    out = od / ("conc_flux_compare.png" if s2 is not None else "conc_flux.png")
    fig.savefig(out, dpi=150); plt.close(fig)
    # CSV
    hdr = "t_min,cIn_mM,cH_mM,Jw_ums,Js"
    rows = np.column_stack([s1["t"], s1["cIn"], s1["cH"], s1["Jw"], s1["Js"]])
    np.savetxt(od / "conc_flux.csv", rows, delimiter=",", header=hdr, comments="")
    print(f"[conc_flux] {rid}: -> {out}")
    return out


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "all"
    rids = list(REPS) if arg == "all" else [arg]
    for rid in rids:
        make(rid)
