#!/usr/bin/env python3
"""Taylor-series vs Donnan B(c) overlay for ANY concentrating NF270 sheet, using
the REAL fitted parameters from the DATA3 NF270 bform_study campaign.  This is the
generalized version of _make_taylor_vs_donnan_lacl3.py — it accepts one or more
run_ids and produces the SAME clean two-panel figure for each.

Truth  : the mechanistic Donnan+dielectric B(c) (solid black).
Approx : poly1 (linear) / poly2 (quadratic) / poly3 (cubic) Taylor forms, plus
         the saturating exp (the bounded practical form).

ROBUSTNESS: any form whose result_<form>.json lacks a "parameters" block (e.g. the
solver returned no feasible fit, as for NaCl poly1 + sat) is silently SKIPPED — it
is never fabricated and never crashes the figure.  So a sheet that only converged
const+quad+cubic shows exactly those vs Donnan.

B(c) math is REUSED verbatim from _make_model_selection.partition_curve (the
analytic Jw-bracket from the fitted coefficients).  The partition curve is
converted to apparent permeability in µm/s by the SAME factor the campaign uses,
B_apparent_ums = B_partition * (Jw*1e4)  [_run_bform_campaign.apparent_partition,
lines 86-87: Bmu = Js*1e4/dc, Bp = Js/(Jw*dc)].  That scalar Jw*1e4 is read straight
out of the sheet's apparent_B.csv (constant across vials because σ pins the flux to
Jw≈Lp·ΔP), so no re-solving is needed and the scaling matches the deck.

Two panels: within measured window | extrapolated to 1.5*c_max.

Usage:
    python3 _make_taylor_vs_donnan.py                  # default sheets (CaCl2, NaCl)
    python3 _make_taylor_vs_donnan.py MC3.07.11.24_SCaCl2 MC5.07.23.24_NaCl
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STUDY = (HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts"
         / "nf270" / "bform_study")
OUTDIR = STUDY / "taylor_vs_donnan"

DEFAULT_SHEETS = ["MC3.07.11.24_SCaCl2", "MC5.07.23.24_NaCl"]

# Approximation forms (Taylor) + the bounded saturating form; Donnan is the truth.
APPROX = ["poly1", "poly2", "poly3", "sat"]
LABEL = {"poly1": "linear (Taylor n=1)", "poly2": "quadratic (Taylor n=2)",
         "poly3": "cubic (Taylor n=3)", "sat": "saturating exp",
         "donnan": "Donnan+dielectric (mechanistic)"}
COLOR = {"poly1": "#1f77b4", "poly2": "#2ca02c", "poly3": "#9467bd",
         "sat": "#ff7f0e", "donnan": "#000000"}

# Pretty salt names for the title.
SALT_TEX = {"CaCl2": "CaCl₂", "NaCl": "NaCl", "LaCl3": "LaCl₃"}


def partition_curve(form, params, c):
    """Analytic partition B(c) (= the Jw bracket) from fitted coefficients.
    Verbatim from _make_model_selection.partition_curve."""
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


def load_form(sdir, form):
    """Return the fitted result dict for `form`, or None if it does not exist or
    has no converged "parameters" block (e.g. {"error": ...})."""
    p = sdir / f"result_{form}.json"
    if not p.exists():
        return None
    r = json.loads(p.read_text())
    if "parameters" not in r or not r.get("parameters"):
        return None        # skip non-converged forms — never fabricate
    return r


def jw_factor_and_window(sdir):
    """Read the (constant) Jw*1e4 partition->µm/s factor and the measured window
    [c_min, c_max] straight from the sheet's apparent_B.csv (the campaign output),
    so the scaling matches the deck exactly."""
    d = np.genfromtxt(sdir / "apparent_B.csv", delimiter=",", names=True)
    c = np.atleast_1d(d["c_mM"])
    bp = np.atleast_1d(d["B_partition"])
    bs = np.atleast_1d(d["B_apparent_ums"])
    jw_fac = float(np.median(bs / bp))   # = Jw * 1e4  (µm/s per partition unit)
    return jw_fac, float(np.nanmin(c)), float(np.nanmax(c))


def make_figure(sheet):
    sdir = STUDY / sheet
    if not sdir.exists():
        raise FileNotFoundError(f"sheet dir not found: {sdir}")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"{sheet}_taylor_vs_donnan.png"

    jw_fac, c_min, c_max = jw_factor_and_window(sdir)

    # concentration grids: within window and modest extrapolation to 1.5*c_max
    c_lo = max(0.1, 0.5 * c_min)
    c_in = np.linspace(c_lo, c_max, 200)
    c_ex = np.linspace(c_lo, 1.5 * c_max, 260)

    donnan = load_form(sdir, "donnan")

    # which approx forms actually converged (skip the rest entirely)
    present = [f for f in APPROX if load_form(sdir, f) is not None]
    skipped = [f for f in APPROX if load_form(sdir, f) is None]

    # fitted-parameter report (donnan + converged approx + single)
    used = {}
    for form in ["donnan"] + APPROX + ["single"]:
        r = load_form(sdir, form)
        if r:
            used[form] = r["parameters"]

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6), sharey=False)

    for ax, (cg, title, is_ext) in zip(
            axes,
            [(c_in, "within measured window", False),
             (c_ex, "extrapolated  (to 1.5 · c_max)", True)]):
        ax.axvspan(c_min, c_max, color="0.92", zorder=0, label="measured window")

        # Donnan truth — solid black
        if donnan:
            yb = partition_curve("donnan", donnan["parameters"], cg) * jw_fac
            ax.plot(cg, yb, "-", lw=3.0, color=COLOR["donnan"],
                    label=LABEL["donnan"], zorder=5)

        # Taylor approximations + saturating (only the converged ones). The
        # saturating-exp curve collapses onto the (flat) Donnan plateau within the
        # window, so plot it DASHED and on top of the black reference — otherwise it
        # hides exactly under the Donnan line and its legend entry looks phantom.
        for form in present:
            r = load_form(sdir, form)
            y = partition_curve(form, r["parameters"], cg) * jw_fac
            if form == "sat":
                ax.plot(cg, y, "--", lw=2.0, color=COLOR[form], alpha=0.95,
                        dashes=(4, 3), zorder=6, label=LABEL[form])
            else:
                ax.plot(cg, y, "-", lw=1.9, color=COLOR[form], alpha=0.95,
                        label=LABEL[form])

        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel("interfacial concentration  c$_{in}$  [mM]")
        ax.set_ylabel("solute permeability  B(c)  [µm/s]")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.25)
        ax.set_xlim(cg.min(), cg.max())

    # keep the within-window panel from being dominated by extrapolation blow-up
    if donnan:
        yb_win = partition_curve("donnan", donnan["parameters"], c_in) * jw_fac
        ymax_win = max(0.6, float(np.nanmax(yb_win)) * 2.2)
        axes[0].set_ylim(-0.05, ymax_win)

    axes[1].legend(fontsize=7, loc="best", framealpha=0.9)

    salt = (donnan or {}).get("salt") or sheet.split("_")[-1]
    salt_tex = SALT_TEX.get(salt, salt)
    rid_short = sheet.split("_")[0]
    fig.suptitle(
        f"{salt_tex} ({rid_short})  —  Taylor-series B(c) vs mechanistic Donnan, "
        f"REAL fitted parameters\n"
        f"Donnan is the mechanistically-grounded, physically-bounded form; the "
        f"polynomials are local Taylor approximations that diverge on extrapolation",
        fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(out, dpi=160)
    plt.close(fig)

    # ---- report ----
    print(f"\n================  {sheet}  ================")
    print(f"saved: {out}")
    print(f"Jw conversion factor (Jw*1e4, µm/s per partition unit) = {jw_fac:.6f}")
    print(f"   => implied Jw = {jw_fac/1e4:.6e} cm/s   (constant across vials)")
    print(f"measured window c_in = [{c_min:.3f}, {c_max:.3f}] mM; "
          f"extrapolated to {1.5*c_max:.3f} mM")
    print("B_apparent[µm/s] = partition_curve(form) * (Jw*1e4)")
    print(f"converged approx forms plotted : {present}")
    print(f"SKIPPED (no feasible fit)      : {skipped if skipped else 'none'}")
    if "single" in used:
        print(f"single-form B@feed (µm/s)      = {used['single'].get('B'):.4f}")
    print(f"=== {sheet} fitted parameters used (per form) ===")
    for form, p in used.items():
        print(f"-- {form} --")
        for k, v in p.items():
            print(f"     {k:8s} = {v}")
    return out


def main():
    sheets = sys.argv[1:] or DEFAULT_SHEETS
    outs = []
    for s in sheets:
        outs.append(make_figure(s))
    print("\n==== ALL FIGURES ====")
    for o in outs:
        print(o)


if __name__ == "__main__":
    main()
