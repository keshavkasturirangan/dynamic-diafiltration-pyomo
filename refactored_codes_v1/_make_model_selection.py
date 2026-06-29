#!/usr/bin/env python3
"""DATA3 'which B(c) form to fit' — model-selection evidence for deck §03, built
from the REAL campaign DAE fits (not the flat constant-B apparent scatter).

In-window quality   : campaign AICc per form (which form best fits the full
                      mass + permeate + retentate data, penalised for parameters).
Shape / extrapolation: the fitted partition curve B(c) of each form, analytic from
                      the fitted coefficients, drawn in-window and extended to
                      1.5*c_max — polynomials diverge, saturating/Donnan plateau.

Outputs (bform_study/model_selection/):
  aicc_table.csv             per salt x form: mean WSSE3, mean AICc, rank, n_ok
  fitted_Bc_curves.png       fitted partition B(c) per salt (in-window + extrapolated)
  Bc_vs_ionic_strength.png   fitted B(c) on a c-axis vs an I-axis (per salt)
Usage: python3 _make_model_selection.py
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270" / "bform_study"
OUT = STUDY / "model_selection"
OUT.mkdir(parents=True, exist_ok=True)

SALT_SHEETS = {
    "NaCl": ["MC3.07.22.24_SNaCl", "MC5.07.23.24_S2NaCl", "MC4.07.11.24_SNaCl",
             "MC2.05.07.24_NaCl", "MC5.07.23.24_NaCl", "MC5.07.23.24_SNaCl"],
    "CaCl2": ["MC2.05.07.24_CaCl2", "MC3.07.11.24_SCaCl2", "MC3.07.12.24_S2CaCl2"],
    "LaCl3": ["MC2.05.21.24_LaCl3", "MC4.07.11.24_SLaCl3"],
}
HEAD = {"NaCl": "MC3.07.22.24_SNaCl", "CaCl2": "MC2.05.07.24_CaCl2", "LaCl3": "MC2.05.21.24_LaCl3"}
KI = {"NaCl": 1.0, "CaCl2": 3.0, "LaCl3": 6.0}
FORMS = ["single", "poly1", "poly2", "poly3", "sat", "donnan"]
LABEL = {"single": "constant", "poly1": "linear", "poly2": "quadratic", "poly3": "cubic",
         "sat": "saturating", "donnan": "Donnan"}
COLOR = {"single": "#6c757d", "poly1": "#1f77b4", "poly2": "#2ca02c", "poly3": "#9467bd",
         "sat": "#ff7f0e", "donnan": "#d62728"}


def res(sheet, form):
    p = STUDY / sheet / f"result_{form}.json"
    if not p.exists():
        return None
    try:
        r = json.loads(p.read_text())
        return None if "error" in r else r
    except Exception:
        return None


def cmax_of(sheet):
    p = STUDY / sheet / "apparent_B.csv"
    if p.exists():
        d = np.genfromtxt(p, delimiter=",", names=True)
        return float(np.nanmax(np.atleast_1d(d["c_mM"])))
    return 100.0


def partition_curve(form, params, c):
    """Analytic partition B(c) (= the Jw bracket) from fitted coefficients."""
    p = params
    if form == "poly1":
        return p["beta_0"] + p["beta_1"] * c
    if form == "poly2":
        return p["beta_0"] + p["beta_1"] * c + p["beta_2"] * c**2
    if form == "poly3":
        return p["beta_0"] + p["beta_1"] * c + p["beta_2"] * c**2 + p["beta_3"] * c**3
    if form == "sat":
        return p["B_inf"] * (1 - np.exp(-c / p["c_star"]))
    if form == "donnan":
        cc = np.maximum(c, 1e-6)
        return p["P0"] * (-p["X"] + np.sqrt(p["X"]**2 + 4 * p["k_dd"]**2 * cc**2)) / (2 * cc)
    return None


def main():
    # ---- AICc table per salt ----
    rows = []
    for salt, sheets in SALT_SHEETS.items():
        for form in FORMS:
            ws, ac, n = [], [], 0
            for s in sheets:
                r = res(s, form)
                if r and r.get("AICc") is not None:
                    ws.append(r["WSSE3"]); ac.append(r["AICc"]); n += 1
            if n:
                rows.append({"salt": salt, "form": LABEL[form],
                             "mean_WSSE3": float(np.mean(ws)), "mean_AICc": float(np.mean(ac)), "n_ok": n})
    # rank within salt by AICc (lower better)
    for salt in SALT_SHEETS:
        sub = [r for r in rows if r["salt"] == salt]
        for rank, r in enumerate(sorted(sub, key=lambda x: x["mean_AICc"]), 1):
            r["rank"] = rank
    with open(OUT / "aicc_table.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["salt", "form", "mean_WSSE3", "mean_AICc", "rank", "n_ok"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("=== in-window AICc ranking (lower AICc = better) ===")
    for salt in SALT_SHEETS:
        print(f"-- {salt} --")
        for r in sorted([x for x in rows if x["salt"] == salt], key=lambda x: x["rank"]):
            print(f"   #{r['rank']} {r['form']:11s} WSSE3={r['mean_WSSE3']:.1f}  AICc={r['mean_AICc']:.0f}  (n={r['n_ok']})")

    # ---- fitted partition B(c) curves (in-window + extrapolated) ----
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3))
    for ax, salt in zip(axes, ["NaCl", "CaCl2", "LaCl3"]):
        sheet = HEAD[salt]
        cmax = cmax_of(sheet)
        c_in = np.linspace(0.5, cmax, 120)
        c_ex = np.linspace(cmax, 1.5 * cmax, 60)
        ax.axvspan(0.5, cmax, color="0.92", zorder=0, label="measured window")
        for form in ["poly1", "poly2", "poly3", "sat", "donnan"]:
            r = res(sheet, form)
            if not r:
                continue
            try:
                y_in = partition_curve(form, r["parameters"], c_in)
                y_ex = partition_curve(form, r["parameters"], c_ex)
                ax.plot(c_in, y_in, "-", lw=2, color=COLOR[form], label=LABEL[form])
                ax.plot(c_ex, y_ex, "--", lw=1.6, color=COLOR[form], alpha=0.8)
            except Exception:
                pass
        ax.set_title(f"{salt}  ({sheet.split('_')[0]})", fontsize=10)
        ax.set_xlabel("interfacial concentration  c_in  [mM]")
        ax.set_ylabel("partition  B(c) / Jw   [—]")
        ax.axhline(0, color="k", lw=0.6)
        ax.legend(fontsize=6); ax.grid(alpha=0.25)
    fig.suptitle("Fitted B(c) partition — solid = measured window, dashed = extrapolation", fontsize=10)
    fig.tight_layout(); fig.savefig(OUT / "fitted_Bc_curves.png", dpi=160); plt.close(fig)

    # ---- B(c) on concentration vs ionic-strength axis (saturating form across salts) ----
    fig2, (axc, axi) = plt.subplots(1, 2, figsize=(11, 4.3))
    sc = {"NaCl": "#1f77b4", "CaCl2": "#2ca02c", "LaCl3": "#d62728"}
    for salt in ["NaCl", "CaCl2", "LaCl3"]:
        sheet = HEAD[salt]; cmax = cmax_of(sheet)
        c = np.linspace(0.5, cmax, 120)
        # prefer the saturating curve; fall back to cubic
        for form in ("sat", "poly3", "poly2"):
            r = res(sheet, form)
            if r:
                y = partition_curve(form, r["parameters"], c)
                axc.plot(c, y, color=sc[salt], lw=2, label=f"{salt} ({LABEL[form]})")
                axi.plot(KI[salt] * c, y, color=sc[salt], lw=2, label=f"{salt}")
                break
    axc.set_xlabel("concentration c [mM]"); axc.set_ylabel("partition B(c) [—]"); axc.set_title("vs concentration"); axc.legend(fontsize=7); axc.grid(alpha=0.25)
    axi.set_xlabel("ionic strength I = ½Σcᵢzᵢ² [mM]"); axi.set_ylabel("partition B(c) [—]"); axi.set_title("vs ionic strength"); axi.legend(fontsize=7); axi.grid(alpha=0.25)
    fig2.tight_layout(); fig2.savefig(OUT / "Bc_vs_ionic_strength.png", dpi=160); plt.close(fig2)
    print("MODELSEL DONE", flush=True)


if __name__ == "__main__":
    main()
