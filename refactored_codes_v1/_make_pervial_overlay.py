#!/usr/bin/env python3
"""DATA3 — empirical B(c_in) from a per-vial fit, with each correlation overlaid.

This is the *correct* "which B–concentration correlation" diagnostic (cf. Xinhong
Liu's run_pre_B_dependence.py):
  * B_form='pervial' floats a FREE B per vial -> the model-agnostic empirical
    apparent permeability B = Js/(c_in-c_H) vs interfacial concentration (the
    reference cloud).  Vials whose B rails to the 50 µm/s bound are flagged
    (open markers) — there B is unidentified, not measured.
  * each candidate correlation (constant / linear / quadratic / cubic / saturating
    / Donnan) is forward-simulated at its campaign fit and its apparent B(c_in)
    drawn through the same cloud.

Output: bform_study/model_selection/pervial_overlay.png  (one panel per salt)
        + per-sheet pervial_overlay_<rid>.png
Usage: python3 _make_pervial_overlay.py
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json
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
OUT = STUDY / "model_selection"
OUT.mkdir(parents=True, exist_ok=True)
NFE = 100

SALT_SHEET = {"NaCl": "MC3.07.22.24_SNaCl", "CaCl2": "MC2.05.07.24_CaCl2", "LaCl3": "MC2.05.21.24_LaCl3"}
FORMS = [("single", "constant", "#6c757d"), (1, "linear", "#1f77b4"), (2, "quadratic", "#2ca02c"),
         (3, "cubic", "#9467bd"), ("sat", "saturating", "#ff7f0e"), ("donnan", "Donnan", "#d62728")]


def apparent_from_sim(sim, ds, per_vial=False):
    """apparent B = Js*1e4/(c_in-c_H) [µm/s].  per_vial -> one (c,B) per vial
    (median over time); else a sorted (c,B) curve over all interior points."""
    nv0 = int(ds["data_config"].get("n_v0", 1))
    cs, bs, rail = [], [], []
    for k, v in sim.items():
        if (k + 1) < nv0:
            continue
        cIn = np.asarray(v["cIn"], float); cH = np.asarray(v["cH"], float); Js = np.asarray(v["Js"], float)
        Bv = v.get("B")
        dc = cIn - cH; m = dc > 1e-6
        if not m.any():
            continue
        Bapp = Js[m] * 1e4 / dc[m]
        if per_vial:
            cs.append(float(np.median(cIn[m]))); bs.append(float(np.median(Bapp)))
            # railed if the fitted per-vial B sits at the 50 bound
            bb = np.median(np.asarray(Bv)[m]) if Bv is not None else np.median(Bapp)
            rail.append(bool(bb >= 49.9))
        else:
            cs += list(cIn[m]); bs += list(Bapp)
    o = np.argsort(cs)
    return np.array(cs)[o], np.array(bs)[o], (np.array(rail)[o] if per_vial else None)


def pervial_empirical(ds):
    fit, sim, _ = lib.solve_model(ds, ds["mode"], sim_opt=False, B_form="pervial",
                                  workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=300, LOUD=False)
    if not isinstance(fit, dict):
        return None
    return apparent_from_sim(sim, ds, per_vial=True)


def form_curve(ds, form, params):
    _, sim, _ = lib.solve_model(ds, ds["mode"], theta=params, sim_opt=True, B_form=form,
                                workflow_family="DATA3", nfe=NFE, solver_max_cpu_time=120, LOUD=False)
    c, b, _ = apparent_from_sim(sim, ds, per_vial=False)
    return c, b


def main():
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    for ax, (salt, rid) in zip(axes, SALT_SHEET.items()):
        ds = camp.load(rid)
        # empirical per-vial cloud
        emp = pervial_empirical(ds)
        if emp is not None:
            c, b, rail = emp
            interior = ~rail
            ax.scatter(c[interior], b[interior], s=55, facecolor="white", edgecolor="black", zorder=6,
                       label="per-vial B (free)")
            if rail.any():
                ax.scatter(c[rail], b[rail], s=55, marker="^", facecolor="none", edgecolor="0.6",
                           zorder=6, label="per-vial (railed → unident.)")
        # correlation curves
        for form, lbl, col in FORMS:
            tag = form if isinstance(form, str) else f"poly{form}"
            rp = STUDY / rid / f"result_{tag}.json"
            if not rp.exists():
                continue
            r = json.loads(rp.read_text())
            if "parameters" not in r or "error" in r:
                continue
            try:
                c2, b2 = form_curve(ds, form, r["parameters"])
                ax.plot(c2, b2, "-", lw=2, color=col, label=lbl, alpha=0.9)
            except Exception as e:
                print(f"  [curve fail] {rid} {tag}: {e}", flush=True)
        ax.set_title(f"{salt}  ({rid.split('_')[0]})", fontsize=10)
        ax.set_xlabel("interfacial concentration  c_in  [mM]")
        ax.set_ylabel("apparent B  =  Jₛ/(c_in−c_H)   [µm s⁻¹]")
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=6.5); ax.grid(alpha=0.25)
        print(f"[overlay] {salt} {rid} done", flush=True)
    fig.suptitle("Empirical per-vial B(c_in) (markers) vs each fitted correlation (lines)", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "pervial_overlay.png", dpi=160)
    plt.close(fig)
    print("PERVIAL OVERLAY DONE", flush=True)


if __name__ == "__main__":
    main()
