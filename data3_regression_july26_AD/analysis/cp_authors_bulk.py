"""
CP Phase 2 (docs/prompts/cp_phase2_regressions.md): apply concentration
polarization on top of the DATA3 authors' (non-CP) bulk interfacial
concentration, for the single-salt point-estimate/uncertainty regressions.

Mechanism (per dataset with an authors' processed sheet):
    c_int_CP = c_p + (c_int_bulk - c_p) * exp(Jw / k_salt)
(permeate assumed unpolarized, so c_p is unchanged -- same thin-film model
Phase 1 already validated in analysis/preprocessing/preprocess_nacl.py).

`k_salt` is imported from Phase 1's `mass_transfer_coefficient()` +
per-salt `Salt.solution_D` -- the single source of truth; NOT recomputed
here. This module only shifts the predictor (c_int); the response
(deltaC_m, built from c_p and Jw) is untouched by design.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "analysis" / "preprocessing"))
import preprocess_nacl as pp  # noqa: E402

SALT_NACL = pp.SALT_NACL
SALT_CACL2 = pp.SALT_CACL2
SALT_LACL3 = pp.SALT_LACL3

# k_salt [m/s], computed once from Phase 1's single source of truth.
K_SALT = {
    "NaCl": pp.mass_transfer_coefficient(D=SALT_NACL.solution_D),
    "CaCl2": pp.mass_transfer_coefficient(D=SALT_CACL2.solution_D),
    "LaCl3": pp.mass_transfer_coefficient(D=SALT_LACL3.solution_D),
}


def apply_cp(c_int_bulk, c_p, Jw, salt: str):
    """Return the CP-corrected interfacial concentration for `salt` in
    {"NaCl", "CaCl2", "LaCl3"}. c_p is returned unchanged by the thin-film
    model (permeate unpolarized) -- only c_int shifts."""
    k = K_SALT[salt]
    c_int_bulk = np.asarray(c_int_bulk, dtype=float)
    c_p = np.asarray(c_p, dtype=float)
    Jw = np.asarray(Jw, dtype=float)
    return c_p + (c_int_bulk - c_p) * np.exp(Jw / k)


if __name__ == "__main__":
    print("k_salt [m/s] (Phase 1 single source of truth):")
    for salt, k in K_SALT.items():
        print(f"  {salt:<6} k = {k:.6e}")
