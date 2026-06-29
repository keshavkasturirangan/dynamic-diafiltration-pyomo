#!/usr/bin/env python3
"""DATA3 B(c) curves + model-selection table from the B-form campaign.

For each sheet:
  * forward-sim every fitted form at its theta, compute apparent B(c)=Js*1e4/(cIn-cH),
  * plot apparent-B points (from the constant fit) with each form's B(c) overlaid
    (figure_s5 palette), -> bform_study/curves/<rid>_Bc.png
For the campaign:
  * aggregate WSSE / AICc / BIC per (sheet,form) -> model_selection.csv + figure,
  * fit each functional form to the pooled apparent-B-vs-c points per salt -> the
    slide-10 'which form to fit' table (R^2 in-window + extrapolation).

Usage: python3 _make_bc_analysis.py [sheets=all|<rid,..>]
Reads/writes only under paper_artifacts/nf270/bform_study/.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib
try:
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
STUDY = ART / "bform_study"
CURVES = STUDY / "curves"
CURVES.mkdir(parents=True, exist_ok=True)
lib.NF270_CF_RESIDUAL_FLOOR_MM = 1.0

FORMS = ["single", 1, 2, 3, "sat", "donnan"]
LABEL = {"single": "constant", "1": "linear", "2": "quadratic", "3": "cubic",
         "sat": "saturating exp", "donnan": "Donnan–dielectric"}
COLOR = {"single": "#6c757d", "1": "#1f77b4", "2": "#2ca02c", "3": "#9467bd",
         "sat": "#ff7f0e", "donnan": "#d62728"}
SHEETS = [
    ("MC3.07.22.24_SNaCl", "NaCl"), ("MC5.07.23.24_S2NaCl", "NaCl"),
    ("MC4.07.11.24_SNaCl", "NaCl"), ("MC2.05.07.24_NaCl", "NaCl"),
    ("MC5.07.23.24_NaCl", "NaCl"), ("MC5.07.23.24_SNaCl", "NaCl"),
    ("MC2.05.07.24_CaCl2", "CaCl2"), ("MC3.07.11.24_SCaCl2", "CaCl2"),
    ("MC3.07.12.24_S2CaCl2", "CaCl2"), ("MC2.05.21.24_LaCl3", "LaCl3"),
    ("MC4.07.11.24_SLaCl3", "LaCl3"),
]


def load_ds(rid):
    fam = lib.NF270_RUN_REGISTRY[rid]
    return lib.loadxlsx(lib.NF270_DEFAULT_ROOT / fam["workbook"], sheet=fam["sheet"])["data_stru"]


def apparent_curve(ds, form, params):
    """Forward-sim form at params; return (cIn[], Bapp[µm/s]) over interior vials."""
    _, sim, _ = lib.solve_model(ds, ds["mode"], theta=params, sim_opt=True, B_form=form,
                                workflow_family="DATA3", nfe=120, solver_max_cpu_time=120, LOUD=False)
    nv0 = int(ds["data_config"].get("n_v0", 1))
    cc, bb = [], []
    for k, v in sim.items():
        if (k + 1) < nv0:
            continue
        cIn = np.asarray(v["cIn"], float); cH = np.asarray(v["cH"], float); Js = np.asarray(v["Js"], float)
        dc = cIn - cH
        m = dc > 1e-6
        if m.any():
            cc.append(float(np.median(cIn[m]))); bb.append(float(np.median(Js[m] * 1e4 / dc[m])))
    o = np.argsort(cc)
    return np.array(cc)[o], np.array(bb)[o]


def sheet_curve(rid, salt):
    sdir = STUDY / rid
    if not sdir.exists():
        return None
    results = {}
    for form in FORMS:
        tag = form if isinstance(form, str) else f"poly{form}"
        rp = sdir / f"result_{tag}.json"
        if rp.exists():
            r = json.loads(rp.read_text())
            if "parameters" in r and "error" not in r:
                results[form] = r
    if not results:
        return None
    ds = load_ds(rid)
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    # observed apparent-B points from the constant fit (the data's realized partition)
    apc = sdir / "apparent_B.csv"
    if apc.exists():
        d = np.genfromtxt(apc, delimiter=",", names=True)
        ax.scatter(np.atleast_1d(d["c_mM"]), np.atleast_1d(d["B_apparent_ums"]),
                   s=42, facecolor="white", edgecolor="black", zorder=5,
                   label="apparent B = Jₛ/(c_in−c_H)")
    rows = []
    for form in FORMS:
        if form not in results:
            continue
        tag = str(form)
        try:
            cc, bb = apparent_curve(ds, form, results[form]["parameters"])
            ax.plot(cc, bb, "-", lw=2.2, color=COLOR[tag], label=LABEL[tag])
        except Exception as e:
            print(f"  [curve fail] {rid} {tag}: {e}", flush=True)
        rr = results[form]
        rows.append({"run_id": rid, "salt": salt, "form": tag,
                     "WSSE": rr.get("WSSE3", rr.get("Obj")), "AICc": rr.get("AICc"), "BIC": rr.get("BIC"),
                     "n_Bparams": rr.get("n_Bparams"),
                     "FIM_cond": (rr.get("FIM") or {}).get("cond")})
    ax.set_xlabel("interfacial concentration  c_in  [mM]")
    ax.set_ylabel("apparent solute permeability  B  [µm s⁻¹]")
    ax.set_title(f"{rid}  ·  {salt}", fontsize=10)
    ax.legend(fontsize=7, framealpha=0.9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(CURVES / f"{rid}_Bc.png", dpi=160)
    plt.close(fig)
    print(f"[curve] {rid}", flush=True)
    return rows


def main():
    sel = sys.argv[1] if len(sys.argv) > 1 else "all"
    sheets = SHEETS if sel == "all" else [s for s in SHEETS if s[0] in sel.split(",")]
    allrows = []
    for rid, salt in sheets:
        try:
            r = sheet_curve(rid, salt)
            if r:
                allrows += r
        except Exception as e:
            print(f"[sheet fail] {rid}: {e}", flush=True)
    # model-selection summary: best (lowest) AICc per sheet
    with open(STUDY / "model_selection.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["run_id", "salt", "form", "WSSE", "AICc", "BIC", "n_Bparams", "FIM_cond"])
        w.writeheader()
        for r in allrows:
            w.writerow(r)
    # per-sheet AICc winner
    by_sheet = {}
    for r in allrows:
        by_sheet.setdefault(r["run_id"], []).append(r)
    print("\n=== AICc winner per sheet ===")
    for rid, rs in by_sheet.items():
        rs2 = [r for r in rs if r.get("AICc") is not None]
        if rs2:
            best = min(rs2, key=lambda r: r["AICc"])
            print(f"  {rid:22s} -> {LABEL.get(best['form'], best['form'])} (AICc={best['AICc']:.1f})")
    print("ANALYSIS DONE", flush=True)


if __name__ == "__main__":
    main()
