#!/usr/bin/env python3
"""
peclet_demo.py — a minimal, self-contained demo of how Pe* is found.

This is the TOY version of the `pe_sweep` step inside _make_peclet.py (lines 79-98).
It (1) builds a synthetic (cIn, cH, Jw, Js) dataset at a KNOWN true Peclet number,
(2) sweeps Pe and, at each Pe, least-squares-fits the partition coefficients (h0,h1),
(3) records the regression MSE for every Pe, and (4) reports
    Pe* = argmin_Pe MSE(Pe)
showing it recovers the true Pe. Run:  python3 peclet_demo.py
"""
import numpy as np

# ------------------------------------------------------------------ THE MODEL
# DATA2 convection-diffusion partition relation (eq 42). With a partition that is
# linear in interface concentration, H = h0 + h1*c, the exact solution of the
# convection-diffusion equation gives, for a GIVEN Peclet number Pe:
#
#     Js/Jw = h0 * x0(Pe) + h1 * x1(Pe)
#       x0(Pe) = (cIn   * e^Pe - cH  ) / (e^Pe - 1)
#       x1(Pe) = (cIn^2 * e^Pe - cH^2) / (e^Pe - 1)
#
# So for a fixed Pe the two basis columns x0, x1 are just numbers, and the model
# is LINEAR in the unknowns (h0, h1) -> ordinary least squares fits them in one shot.
def basis(Pe, cIn, cH):
    e = np.exp(np.minimum(Pe, 700.0))          # cap exponent to avoid overflow
    d = e - 1.0
    x0 = (cIn * e - cH) / d
    x1 = (cIn ** 2 * e - cH ** 2) / d
    return np.column_stack([x0, x1])           # design matrix, shape (N, 2)


# --------------------------------------------- 1) SYNTHETIC DATA AT A KNOWN Pe
rng = np.random.default_rng(0)
Pe_true, h0_true, h1_true = 0.6, 0.40, 0.010   # the ground truth we will try to recover
# NOTE: Pe~0.6 sits in the IDENTIFIABLE range. Above Pe~3 the basis x0,x1 -> cIn,cIn^2
# regardless of Pe, so the MSE valley goes flat — that is why a convection-dominated
# experiment only pins "Pe >> 1" loosely, never a sharp value (the real-data finding).
N = 60
cIn = np.linspace(2.0, 25.0, N)                # interface concentration [mM]
cH = np.full(N, 1.0)                           # permeate-side concentration [mM].
# IMPORTANT: cH must carry information that is NOT just proportional to cIn. If cH = a*cIn
# then Js/Jw collapses to a pure polynomial in cIn that ANY Pe can fit (the basis at low Pe
# is already {cIn, cIn^2}) — Pe becomes unidentifiable. Real permeate data breaks that.
Jw = np.full(N, 1.0)                           # water flux (constant here for clarity)
ratio_true = basis(Pe_true, cIn, cH) @ np.array([h0_true, h1_true])   # true Js/Jw
Js = ratio_true * (1.0 + 0.003 * rng.standard_normal(N)) * Jw         # +0.3% noise -> "measured" Js


# -------------------------------- 2-3) AT EACH Pe: FIT (h0,h1) BY OLS, GET MSE
def fit_mse(Pe):
    y = Js / Jw                                # regression target  Js/Jw
    X = basis(Pe, cIn, cH)                     # basis evaluated at this Pe
    h, *_ = np.linalg.lstsq(X, y, rcond=None)  # best (h0,h1) at this Pe (one OLS solve)
    resid = y - X @ h                          # residuals
    return float(np.mean(resid ** 2)), h       # MSE(Pe), fitted coeffs


Pes = np.logspace(-3, 2, 200)                  # the "dial" settings we try (1e-3 .. 1e2)
mse = np.array([fit_mse(Pe)[0] for Pe in Pes]) # the MSE-vs-Pe curve


# --------------------------------------------------------- 4) Pe* = argmin MSE
k = int(np.argmin(mse))                        # index of the smallest MSE
Pe_star = Pes[k]                               # <-- argmin: the Pe at the bottom
mse_star, (h0_hat, h1_hat) = fit_mse(Pe_star)

print(f"true Pe        = {Pe_true}")
print(f"Pe*  (argmin)  = {Pe_star:.4g}    [min of MSE over {len(Pes)} trial Pe]")
print(f"MSE(Pe*)       = {mse_star:.3e}")
print(f"fitted h0,h1   = {h0_hat:.4f}, {h1_hat:.4f}   (true {h0_true}, {h1_true})")

# optional: save the MSE-vs-Pe curve (this is the 'companion MSE figure')
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.figure(figsize=(6, 4))
    plt.semilogx(Pes, mse, lw=2)
    plt.axvline(Pe_star, ls="--", c="k")
    plt.plot(Pe_star, mse_star, "r*", ms=15, label=f"Pe* = {Pe_star:.2f}")
    plt.axvline(Pe_true, ls=":", c="0.5", label=f"true Pe = {Pe_true}")
    plt.xlabel("Peclet number Pe  (the dial we sweep)")
    plt.ylabel("regression MSE of (h0,h1)")
    plt.title("Pe* = argmin of MSE(Pe)")
    plt.legend(); plt.tight_layout()
    plt.savefig("peclet_demo.png", dpi=140)
    print("saved peclet_demo.png  (the MSE-vs-Pe curve, with Pe* marked)")
except Exception as e:
    print("(plot skipped:", e, ")")
