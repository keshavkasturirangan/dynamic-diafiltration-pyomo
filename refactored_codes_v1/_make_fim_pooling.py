#!/usr/bin/env python3
"""DATA3 §7 headline — FIM-based parameter estimability: individual vs pooled std errors.

For the POORLY-IDENTIFIED groups (conc NaCl, CaCl2, LaCl3 — diluting NaCl excluded: it is
well-identified from a single sheet), compute the Fisher Information Matrix at each sheet's
constant-B fit (2% no-floor covariance) and the marginal parameter std errors σ_θ = sqrt
(diag(FIM^-1)) for (Lp, B, σ).  Pooling SHARES the membrane properties Lp and σ across the
group (same NF270 coupon, same salt-regime): the pooled shared-parameter information is the
SUM of the per-sheet (profiled) informations, so

    std_pooled(θ_shared) = 1 / sqrt( Σ_i 1/std_i(θ)^2 )   ≤  min_i std_i(θ)

a rigorous ~1/sqrt(N) estimability gain.  σ rails at the σ=1 bound on several sheets, so its
local FIM std is a LOWER BOUND — the shared-σ PROFILE (pooling/sharesigma) is the robust σ view.

Output: bform_study/pooling/fim_estimability.png  (+ fim_estimability.json)
"""
import sys, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import _run_profile_contours as pc

GROUPS = ["NaCl_concentrating", "CaCl2", "LaCl3"]   # diluting NaCl excluded (well-identified)
PNAMES = ["Lp", "B", "sigma"]                        # first 3 FIM params


def _study():
    return pc._study()


def _fim_std(rid):
    """calc_FIM at the sheet's constant-B fit (2% no-floor) → marginal std for [Lp,B,sigma]."""
    import _run_bform_campaign as camp
    import refactored_ucb_library as lib
    import contextlib, io
    lib.NF270_CF_RESIDUAL_SCALE_FRACTION = 0.02
    lib.NF270_CF_RESIDUAL_FLOOR_MM = None
    theta = json.loads((_study() / rid / "result_single.json").read_text())["parameters"]
    ds = camp.load(rid); mode = ds["mode"]
    with contextlib.redirect_stdout(io.StringIO()):
        doe = lib.calc_FIM(ds, mode, theta=dict(theta), B_form="single",
                           workflow_family="DATA3", nfe=140)
    FIM = np.asarray(doe["FIM"], float) if "FIM" in doe else None
    if FIM is None:
        return None
    try:
        cov = np.linalg.inv(FIM)
        std = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        return None
    return {"std": std[:3].tolist(), "sigma_fit": float(theta.get("sigma", 1.0))}


def main():
    study = _study()
    out = {}
    for g in GROUPS:
        sheets = [s for s in pc.POOL_GROUPS[g]
                  if (study / s / "result_single.json").exists()]
        recs = {}
        for s in sheets:
            r = _fim_std(s)
            if r:
                recs[s] = r
                print(f"[fim] {g}/{s.split('_')[0]}: std Lp={r['std'][0]:.3g} B={r['std'][1]:.3g} σ={r['std'][2]:.3g}  (σ_fit={r['sigma_fit']:.2f})", flush=True)
        if not recs:
            continue
        # pooled (shared Lp & σ): 1/sqrt(Σ 1/std^2); B kept per-sheet (concentration-dependent)
        pooled = {}
        for k, name in enumerate(PNAMES):
            inv = [1.0 / r["std"][k] ** 2 for r in recs.values() if r["std"][k] > 0]
            pooled[name] = float(1.0 / np.sqrt(sum(inv))) if inv else None
        out[g] = {"individual": {s.split('_')[0]: recs[s]["std"] for s in recs},
                  "pooled_shared": pooled, "n": len(recs)}
        print(f"[fim] {g} POOLED(shared Lp,σ): Lp={pooled['Lp']:.3g} σ={pooled['sigma']:.3g}  (N={len(recs)})", flush=True)

    (study / "pooling").mkdir(parents=True, exist_ok=True)
    (study / "pooling" / "fim_estimability.json").write_text(json.dumps(out, indent=2))

    # figure: per group, σ std error — individual bars + pooled bar
    fig, axes = plt.subplots(1, len(out), figsize=(5.2 * len(out), 4.4))
    if len(out) == 1:
        axes = [axes]
    for ax, (g, d) in zip(axes, out.items()):
        labels = list(d["individual"]) + ["POOLED"]
        vals = [d["individual"][s][2] for s in d["individual"]] + [d["pooled_shared"]["sigma"]]
        colors = ["0.6"] * (len(labels) - 1) + ["tab:red"]
        ax.bar(range(len(labels)), vals, color=colors)
        ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("σ std error  [-]"); ax.set_title(f"{g}  (N={d['n']})", fontsize=10)
        best_ind = min(d["individual"][s][2] for s in d["individual"])
        ax.text(0.5, 0.92, f"pooled {d['pooled_shared']['sigma']:.3g}  vs best individual {best_ind:.3g}  ({best_ind/d['pooled_shared']['sigma']:.1f}× tighter)",
                transform=ax.transAxes, ha="center", fontsize=7.5, color="tab:red")
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Stacking improves estimability — FIM σ std error shrinks ~1/√N with shared-σ pooling (membrane-salt property).  "
                 "σ at the σ=1 bound → local std is a lower bound; see shared-σ profile for the robust view.", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(study / "pooling" / "fim_estimability.png", dpi=150); plt.close(fig)
    print(f"[fim] -> {study/'pooling'/'fim_estimability.png'}")


if __name__ == "__main__":
    main()
