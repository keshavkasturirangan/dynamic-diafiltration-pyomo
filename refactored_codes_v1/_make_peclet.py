#!/usr/bin/env python3
"""DATA3 Peclet-number analysis — replicates DATA2_main.pdf Figures 7 & 8.

WHY: the DATA2 transport derivation shows the solute-permeability regime is set by the
Peclet number Pe = Jw*L/Dm (eq 25):
  * Pe -> 0  (diffusion-dominated): B = Dm*H/L = CONSTANT          (eq 28)
  * Pe ~ 1+  (convection + diffusion): B becomes CONCENTRATION-DEPENDENT,
             B = Jw * sum_i beta_i c^i  (the Taylor form, eqs 29-31)
So locating Pe is the mechanistic explanation for WHY B depends on concentration.

METHOD (DATA2 Fig 8): from each experiment's fitted model, forward-simulate the series
(Jw, c_in,f=cIn, c_in,p~=cH, Js).  Assume the partition varies linearly with interface conc,
H_f=h0+h1*cIn, H_p=h0+h1*cH (eq 41); substitute into the exact eq-25 solution to get eq 42,
  Js/Jw = h0*[ (cIn*e^Pe - cH)/(e^Pe-1) ] + h1*[ (cIn^2*e^Pe - cH^2)/(e^Pe-1) ],
which is LINEAR in (h0,h1) for a GIVEN Pe.  Sweep Pe in [1e-3, 1e2]; at each Pe linear-regress
(h0,h1) and record the MSE.  The Pe that minimizes MSE = the transport regime of the experiment.

OUTPUT (bform_study/peclet/):
  peclet_mse.png        MSE(Pe) per headline group (+ regressed h0,h1 vs Pe) -> DATA2 Fig 8
  peclet_jsfit.png      Js vs cIn: simulated vs the convection-diffusion fit -> DATA2 Fig 7C
  peclet.json          per group: Pe*, MSE(Pe*), h0,h1 at Pe*, regime label
Usage: python3 _make_peclet.py
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, contextlib, io
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

OUT = camp.OUT / "peclet"
OUT.mkdir(parents=True, exist_ok=True)
HEAD = {"NaCl diluting": "MC3.07.22.24_SNaCl", "NaCl concentrating": "MC2.05.07.24_NaCl",
        "CaCl2": "MC2.05.07.24_CaCl2", "LaCl3": "MC2.05.21.24_LaCl3"}
# B(c) form used to generate the (Jw,cIn,cH,Js) series — prefer the richest converged form
FORM_PREF = [(2, "poly2"), (1, "poly1"), ("single", "single")]
COLOR = {"NaCl diluting": "#1f77b4", "NaCl concentrating": "#4aa3df", "CaCl2": "#2ca02c", "LaCl3": "#d62728"}


def series(rid):
    """Forward-sim the best converged B(c) fit -> concatenated interior (cIn, cH, Jw, Js)."""
    for form, tag in FORM_PREF:
        rp = camp.OUT / rid / f"result_{tag}.json"
        if not rp.exists():
            continue
        r = json.loads(rp.read_text())
        if "parameters" not in r:
            continue
        ds = camp.load(rid)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                fit, sim, _ = lib.solve_model(ds, ds["mode"], theta=r["parameters"], sim_opt=True,
                                              B_form=form, workflow_family="DATA3", nfe=120,
                                              solver_max_cpu_time=120)
        except Exception:
            continue
        if not isinstance(sim, dict):
            continue
        nv0 = int(ds["data_config"].get("n_v0", 1))
        cIn, cH, Jw, Js = [], [], [], []
        for k, v in sim.items():
            if (k + 1) < nv0:
                continue
            cIn += list(np.asarray(v["cIn"], float)); cH += list(np.asarray(v["cH"], float))
            Jw += list(np.asarray(v["Jw"], float)); Js += list(np.asarray(v["Js"], float))
        return tag, np.array(cIn), np.array(cH), np.array(Jw), np.array(Js)
    return None


def pe_sweep(cIn, cH, Jw, Js, Pes):
    """At each Pe, OLS-regress y=Js/Jw on [x0,x1] (eq 42); return MSE(Pe), h0(Pe), h1(Pe)."""
    y = Js / np.where(np.abs(Jw) > 1e-12, Jw, np.nan)
    m = np.isfinite(y) & np.isfinite(cIn) & np.isfinite(cH)
    y, ci, ch = y[m], cIn[m], cH[m]
    mse, H0, H1 = [], [], []
    for Pe in Pes:
        e = np.exp(min(Pe, 700.0)); d = e - 1.0
        if abs(d) < 1e-300:
            d = 1e-300
        x0 = (ci * e - ch) / d
        x1 = (ci ** 2 * e - ch ** 2) / d
        X = np.column_stack([x0, x1])
        try:
            h, *_ = np.linalg.lstsq(X, y, rcond=None)
            res = y - X @ h
            mse.append(float(np.mean(res ** 2))); H0.append(float(h[0])); H1.append(float(h[1]))
        except Exception:
            mse.append(np.nan); H0.append(np.nan); H1.append(np.nan)
    return np.array(mse), np.array(H0), np.array(H1)


def plot_regimes(summary, out_png=None):
    """Single-axis regime figure: each salt's operating Pe* on a LOG Pe axis, with the
    diffusion region (Pe<1) and convection region (Pe>1) shaded and a Pe=1 divider.
    `summary` is the per-group dict built in main() (keys = HEAD labels, with 'Pe_star')."""
    if out_png is None:
        out_png = OUT / "peclet_regimes.png"
    # nicer salt labels for the plotted points
    SALT = {"NaCl diluting": "diluting NaCl", "NaCl concentrating": "conc NaCl",
            "CaCl2": "conc CaCl₂", "LaCl3": "conc LaCl₃"}
    fig, ax = plt.subplots(figsize=(10, 3.6))
    ax.set_xscale("log")
    xlo, xhi = 1e-4, 1e3
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(0, 1)
    # shaded regions: diffusion (Pe<1) and convection (Pe>1)
    ax.axvspan(xlo, 1.0, color="#cfe3f5", alpha=0.6, zorder=0)   # diffusion
    ax.axvspan(1.0, xhi, color="#fde0d0", alpha=0.6, zorder=0)   # convection
    ax.axvline(1.0, color="k", lw=1.6, zorder=1)
    ax.text(10 ** (0.5 * (np.log10(xlo) + 0)), 0.95, "diffusion region\n(Pe < 1)",
            ha="center", va="top", fontsize=10, color="#1f4e79")
    # anchor the convection header just right of Pe=1 (clear of the Pe≈30-80 points)
    ax.text(10 ** 0.4, 0.95, "convection region\n(Pe > 1)",
            ha="center", va="top", fontsize=10, color="#a5421a")
    # plot each salt's Pe* as a point; stagger label heights (and nudge x for
    # points that share a Pe) so labels never overlap.
    items = [(lbl, summary[lbl]) for lbl in HEAD if lbl in summary and np.isfinite(summary[lbl].get("Pe_star", np.nan))]
    # ladder of (y, va, dx_dec) label slots, cycled per point
    SLOTS = [(0.70, "bottom", 0.0), (0.30, "top", 0.0),
             (0.88, "bottom", 0.0), (0.12, "top", 0.0)]
    seen = {}                              # round(log10(pe),2) -> count, to spread coincident points
    for i, (lbl, d) in enumerate(items):
        pe = d["Pe_star"]
        col = COLOR[lbl]
        key = round(np.log10(pe), 2)
        k = seen.get(key, 0); seen[key] = k + 1
        # stagger the MARKERS (not just the labels) when two salts share a Pe* —
        # e.g. CaCl₂ and LaCl₃ both sit at Pe*≈36, so without this they print one
        # on top of the other and read as a single (duplicated) point.
        y_pt = 0.5 + 0.13 * k
        x_dot = pe * (10 ** (0.05 * (k if k % 2 == 0 else -k)))
        ax.plot(x_dot, y_pt, "o", ms=13, color=col, mec="k", mew=1.2, zorder=3)
        yt, va, _ = SLOTS[i % len(SLOTS)]
        x_lab = pe * (10 ** (0.20 * (k if k % 2 == 0 else -k)))   # fan out shared-x labels
        ax.annotate(f"{SALT.get(lbl, lbl)}\nPe*={pe:.2g}", xy=(x_dot, y_pt), xytext=(x_lab, yt),
                    ha="center", va=va, fontsize=10, fontweight="bold", color=col,
                    arrowprops=dict(arrowstyle="-", color=col, lw=1.0), zorder=3)
    ax.set_yticks([])
    ax.set_xlabel("Péclet number  Pe* = J_w·L/D_m  (log scale)", fontsize=11)
    ax.set_title("Transport regime by Péclet number", fontsize=13, fontweight="bold")
    ax.grid(axis="x", which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print("PECLET REGIMES ->", out_png, flush=True)
    return out_png


def main():
    Pes = np.logspace(-3, 2, 90)
    summary = {}
    fig1, axes1 = plt.subplots(1, 4, figsize=(17, 4.2))   # MSE(Pe) per group (DATA2 Fig 8A)
    fig2, axes2 = plt.subplots(1, 4, figsize=(17, 4.2))   # Js vs cIn fit (DATA2 Fig 7C)
    for axm, axj, (lbl, rid) in zip(axes1, axes2, HEAD.items()):
        s = series(rid)
        if s is None:
            axm.set_title(f"{lbl}\n(no fit)"); axj.set_title(f"{lbl}\n(no fit)"); continue
        tag, cIn, cH, Jw, Js = s
        mse, H0, H1 = pe_sweep(cIn, cH, Jw, Js, Pes)
        # Pe* = global min of MSE (prefer the convection-diffusion branch Pe>=0.1 if it is near-best)
        kmin = int(np.nanargmin(mse))
        Pe_star = Pes[kmin]
        regime = "convection+diffusion (B depends on c)" if Pe_star >= 0.3 else "diffusion-dominated (B ~ const)"
        summary[lbl] = {"rid": rid, "series_form": tag, "Pe_star": float(Pe_star),
                        "mse_star": float(mse[kmin]), "h0": float(H0[kmin]), "h1": float(H1[kmin]),
                        "regime": regime}
        # --- Fig 8A: MSE vs Pe ---
        axm.semilogx(Pes, mse, color=COLOR[lbl], lw=2)
        axm.axvline(Pe_star, color="k", ls="--", lw=1)
        axm.plot(Pe_star, mse[kmin], "r*", ms=14)
        axm.set_title(f"{lbl}\nPe*={Pe_star:.2g}  ({'conv-diff' if Pe_star>=0.3 else 'diffusion'})", fontsize=9)
        axm.set_xlabel("Peclet number  Pe = Jw·L/Dm"); axm.grid(alpha=0.3)
        if lbl == list(HEAD)[0]:
            axm.set_ylabel("regression MSE of (h₀,h₁)\n(eq 42)")
        # --- Fig 7C: Js vs cIn, simulated vs convection-diffusion fit at Pe* ---
        e = np.exp(min(Pe_star, 700.0)); d = e - 1.0
        x0 = (cIn * e - cH) / d; x1 = (cIn ** 2 * e - cH ** 2) / d
        Js_fit = (H0[kmin] * x0 + H1[kmin] * x1) * Jw
        o = np.argsort(cIn)
        axj.scatter(cIn, Js, s=6, color="0.5", label="simulated Jₛ", alpha=0.5)
        axj.plot(cIn[o], Js_fit[o], color=COLOR[lbl], lw=2, label="conv-diff fit (Pe*)")
        axj.set_title(f"{lbl}", fontsize=9); axj.set_xlabel("interfacial conc  c_in,f  [mM]")
        axj.legend(fontsize=7); axj.grid(alpha=0.3)
        if lbl == list(HEAD)[0]:
            axj.set_ylabel("solute flux  Jₛ")
        print(f"[peclet] {lbl:20s} ({rid}, {tag}): Pe*={Pe_star:.3g}  {regime}", flush=True)
    fig1.suptitle("DATA3 Peclet analysis (≙ DATA2 Fig 8A): regression MSE vs Pe — a minimum at Pe≳1 means "
                  "convection+diffusion, i.e. B is concentration-dependent (not the Pe→0 constant-B limit)", fontsize=10)
    fig1.tight_layout(rect=[0, 0, 1, 0.93]); fig1.savefig(OUT / "peclet_mse.png", dpi=150); plt.close(fig1)
    fig2.suptitle("DATA3 (≙ DATA2 Fig 7C): simulated Jₛ vs the exact convection-diffusion solution at Pe*", fontsize=10)
    fig2.tight_layout(rect=[0, 0, 1, 0.93]); fig2.savefig(OUT / "peclet_jsfit.png", dpi=150); plt.close(fig2)
    json.dump(summary, open(OUT / "peclet.json", "w"), indent=2)
    plot_regimes(summary)   # one-axis Pe* regime summary (peclet_regimes.png)
    print("PECLET DONE ->", OUT / "peclet_mse.png", flush=True)


if __name__ == "__main__":
    main()
