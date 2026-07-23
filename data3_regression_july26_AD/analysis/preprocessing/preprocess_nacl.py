"""
Salt-general raw-to-processed preprocessing pipeline (NaCl 1:1, CaCl2 2:1,
LaCl3 3:1).

Converts a RAW diafiltration campaign sheet (time series of permeate mass,
retentate/permeate conductivity, plus a per-vial ICP-OES block) into the
PROCESSED per-vial analysis columns consumed by the regressions
(analysis/step2_nacl/, analysis/nacl_all_datasets/, analysis/cacl2/,
analysis/lacl3/) -- interfacial concentration (B), ionic strength (D),
rejection (F), permeate interfacial concentration (H), water flux (J), and
the solute permeability coefficient (K).

Module name kept as `preprocess_nacl` for backward compatibility with
existing importers (process_raw_runs.py, validate.py, make_figures.py, the
NaCl regression/appendix scripts) -- all of those call the pipeline with no
`salt=` argument and so keep getting the exact NaCl behavior they always
have (see the NaCl regression-test in the module's __main__ block and
docs/prompts/cp_preprocessing_extend.md, Phase 1). The `Salt` dataclass
below is what makes the pipeline salt-general; CaCl2/LaCl3 callers
(process_raw_multisalt.py) pass `salt=SALT_CACL2`/`SALT_LACL3` explicitly.

Recipe (see docs/reports/main.tex, Sec. "Data preprocessing", label
sec:preproc, for the same equations written up for the report):

  (a) Water flux from a SMOOTHED slope of cumulative permeate mass m(t):
        Jw = (1/(rho*A_m)) * dm/dt                                   [m/s]
  (b) Bulk concentration from raw conductivity via a linear calibration
      fit per-experiment from that run's own vial (conductivity, ICP) pairs:
        c_bulk = (kappa - intercept) / slope
  (c) Interfacial (retentate-side) concentration via thin-film CP:
        c_int = c_p + (c_bulk_ret - c_p) * exp(Jw / k)
      k = 0.23 * v0^0.57 * D^0.67 / (nu^0.24 * b^0.43), v0 = omega*b/2,
      omega = 2*pi*rpm/60. D is the SOLUTION-phase ambipolar (Nernst-Hartley)
      diffusivity of the salt (per-salt, see Salt.solution_D below), NOT the
      membrane-phase D the regressions use for B. Permeate is assumed
      unpolarized: c_p = bulk permeate concentration directly (no boundary
      layer on that side). If apply_cp=False, bulk == interfacial for both
      sides (matches the "(2)"-sheet / non-CP convention used by Step 2/
      Step 3 and the original CaCl2/LaCl3 analyses).
  (d) Derived: I = 0.5*sum(z_i^2 c_i) -- salt-general via `Salt.i_multiplier`
      (I = 1*c_int for NaCl, 3*c_int for CaCl2, 6*c_int for LaCl3),
      R = 1 - c_p/c_int, B_coeff = Jw*c_p/(c_int-c_p) * 1e6  [um/s].
      Preprocessing emits SALT concentrations (c_int_mM, c_p_mM); the co-ion
      (Cl-) stoichiometry is a regression-model concern (Phase 2), not
      applied here.

WORKING CHOICES / FLAGS (see analysis/preprocessing/README.md for the full
discussion and validation results):
  1. A_m and rho are not documented upstream; A_m is *inferred* by matching
     the reference Jw column (see README). rho = 1000 g/L (water) assumed.
     Salt-independent (same stirred cell for all three salts).
  2. Jw smoothing method/window is undocumented; we use a centered
     rolling-window linear-regression slope of mass vs. time, with the
     window tuned against the reference J column (see validate.py).
  3. Conductivity calibration is fit per-experiment from vial (conductivity,
     ICP) pairs where available; ICP mg/L is converted to mM using the
     ELEMENTAL cation molar mass (Na+ 22.99, Ca2+ 40.08, La3+ 138.91 g/mol),
     since ICP-OES measures the elemental metal, not the molecular salt --
     see README for the validation that supports this choice (NaCl) and its
     direct extension (moles metal = moles salt) to CaCl2/LaCl3.
  4. CP vs non-CP is a caller-supplied flag (apply_cp); Step 2/Step 3, the
     original CaCl2/LaCl3 analyses, and analysis/nacl_all_datasets/ used the
     NON-CP convention (bulk == interfacial) -- CP is now computed in this
     preprocessing layer for all three salts (docs/prompts/
     cp_preprocessing_extend.md Phase 1); re-pointing the regressions at the
     CP-corrected CSVs is Phase 2, deliberately not done here.
  5. H/permeate-side smoothing (incl. the ~0.04 mM floor seen for the first
     ~14-25 rows of the reference sheets, most likely sensor dead-volume
     before permeate first reaches the conductivity cell) is approximated
     with a short floor + rolling-median smooth; not reproduced exactly.

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/preprocessing/preprocess_nacl.py   # smoke-test / demo (NaCl)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "data"

# ---------------------------------------------------------------------------
# Physical constants (working choices -- see module docstring / README)
# ---------------------------------------------------------------------------
RHO_G_PER_L = 1000.0                     # water density [g/L]
A_M_CM2 = 3.95                           # stirred-cell exposed membrane area [cm^2] (inferred)
A_M_M2 = A_M_CM2 * 1e-4                  # [m^2]

NU_M2_S = 1.003e-6                       # kinematic viscosity [m^2/s]
B_CELL_M = 0.0254                        # stirred-cell bore diameter [m]
RPM = 350.0

# ---------------------------------------------------------------------------
# Per-salt constants (Salt-general pipeline; docs/prompts/
# cp_preprocessing_extend.md Phase 1)
# ---------------------------------------------------------------------------
# Solution-phase (NOT membrane-phase) ion diffusivities, m^2/s -- these feed
# the CP boundary-layer mass-transfer coefficient k, distinct from the
# membrane-phase D's the regressions use for B (which are these values / 1000
# -- analysis/step2_nacl/regress_nacl.py's D_NA_M/D_CL_M, analysis/cacl2/
# regress_cacl2.py's D_CA_M/D_CL_M, analysis/lacl3/regress_lacl3.py's
# D_LA_M/D_CL_M -- one source of truth, see cp_preprocessing_extend.md).
D_NA_M2_S = 1.33e-9   # Na+ solution diffusivity
D_CA_M2_S = 0.79e-9   # Ca2+ solution diffusivity
D_LA_M2_S = 0.626e-9  # La3+ solution diffusivity
D_CL_M2_S = 2.03e-9   # Cl- solution diffusivity (shared anion, all three salts)

MOLAR_MASS_NA = 22.99     # g/mol, elemental Na+ (ICP-OES measures elemental metal)
MOLAR_MASS_CA = 40.08     # g/mol, elemental Ca2+
MOLAR_MASS_LA = 138.91    # g/mol, elemental La3+
MOLAR_MASS_NACL = 58.44   # g/mol, molecular NaCl (NOT used; kept for reference)


@dataclass(frozen=True)
class Salt:
    """Per-salt constants for a fully-dissociated binary electrolyte
    A_{z_anion} B_{z_cation} (cation A^{z_cation+}, anion B^{z_anion-},
    minimal-integer charge-balanced stoichiometry): NaCl (z_cation=1,
    z_anion=1), CaCl2 (2,1), LaCl3 (3,1), Na2SO4 (1,2 -- the FIRST salt in
    this project where the co-ion, SO4^2-, is the divalent species rather
    than the cation; z_anion=1 remains the default for the Cl- salts).

    General relations (charge balance fixes `cation_stoich` = z_anion,
    `anion_stoich` = z_cation cations/anions per formula unit):
        I = 0.5 * z_cation * z_anion * (z_cation + z_anion) * c_salt
        D_s = (z_cation+z_anion) * D_cation * D_anion
              / (z_cation*D_cation + z_anion*D_anion)
    which reduce to the familiar closed forms at z_anion=1:
        NaCl (1,1): I = c,  D_s = 2 D_Na D_Cl / (D_Na + D_Cl)
        CaCl2 (2,1): I = 3c, D_s = 3 D_Ca D_Cl / (2 D_Ca + D_Cl)
        LaCl3 (3,1): I = 6c, D_s = 4 D_La D_Cl / (3 D_La + D_Cl)
    and give, at z_cation=1, z_anion=2 (Na2SO4):
        I = 3c,  D_s = 3 D_Na D_SO4 / (D_Na + 2 D_SO4)
    `icp_molar_mass` is always the ELEMENTAL CATION's (ICP-OES measures the
    metal); `cation_stoich` (= z_anion) divides the ICP-derived cation
    concentration down to salt concentration (1 for the Cl- salts, since
    each carries exactly 1 cation per formula unit; 2 for Na2SO4, which
    carries 2 Na+ per formula unit).
    """
    name: str
    z_cation: int
    icp_molar_mass: float          # g/mol, elemental cation
    d_cation_m2_s: float           # solution-phase cation diffusivity [m^2/s]
    z_anion: int = 1
    d_anion_m2_s: float = D_CL_M2_S

    @property
    def cation_stoich(self) -> int:
        """Cations per formula unit (= z_anion, by charge balance)."""
        return self.z_anion

    @property
    def i_multiplier(self) -> float:
        za, zc = self.z_anion, self.z_cation
        return 0.5 * zc * za * (zc + za)

    @property
    def solution_D(self) -> float:
        za, zc = self.z_anion, self.z_cation
        return (zc + za) * self.d_cation_m2_s * self.d_anion_m2_s / (
            zc * self.d_cation_m2_s + za * self.d_anion_m2_s
        )


SALT_NACL = Salt(name="NaCl", z_cation=1, icp_molar_mass=MOLAR_MASS_NA, d_cation_m2_s=D_NA_M2_S)
SALT_CACL2 = Salt(name="CaCl2", z_cation=2, icp_molar_mass=MOLAR_MASS_CA, d_cation_m2_s=D_CA_M2_S)
SALT_LACL3 = Salt(name="LaCl3", z_cation=3, icp_molar_mass=MOLAR_MASS_LA, d_cation_m2_s=D_LA_M2_S)
# Na2SO4: SO4^2- is the co-ion (divalent), Na+ the counter-ion -- the first
# salt in this project where the divalent species is the ANION, not the
# cation (see docs/prompts/copolymer_nacl_na2so4.md). Solution diffusivities
# D_Na=1.33e-9, D_SO4=1.065e-9 m^2/s per that prompt (D_Na matches D_NA_M2_S
# above; D_SO4 is new).
D_SO4_M2_S = 1.065e-9
SALT_NA2SO4 = Salt(name="Na2SO4", z_cation=1, icp_molar_mass=MOLAR_MASS_NA,
                    d_cation_m2_s=D_NA_M2_S, z_anion=2, d_anion_m2_s=D_SO4_M2_S)

# Ambipolar (Nernst-Hartley) NaCl solution diffusivity -- kept as a
# module-level constant (== SALT_NACL.solution_D, verified equal in the
# NaCl regression test) since it is `mass_transfer_coefficient`'s default and
# predates the Salt generalization.
D_NACL_M2_S = 2 * D_NA_M2_S * D_CL_M2_S / (D_NA_M2_S + D_CL_M2_S)

ICP_MOLAR_MASS = MOLAR_MASS_NA  # working choice -- see README validation

DEFAULT_JW_WINDOW = 15           # rows, centered rolling window (odd)
DEFAULT_HFLOOR_ROWS = 20         # rows treated as permeate-side sensor dead time
DEFAULT_HFLOOR_MM = 0.04         # mM floor applied during dead time
DEFAULT_H_SMOOTH_WINDOW = 7      # rolling-median window for permeate-side smoothing


def mass_transfer_coefficient(D=D_NACL_M2_S, rpm=RPM, b=B_CELL_M, nu=NU_M2_S):
    """Stirred-cell mass-transfer coefficient k [m/s] (Leveque-type
    correlation used throughout this project, docs/reports/main.tex sec:cp).
    Pass `D=salt.solution_D` for CaCl2/LaCl3 (NaCl's ambipolar D is the
    default, for backward compatibility)."""
    omega = 2 * math.pi * rpm / 60.0
    v0 = omega * b / 2.0
    k = 0.23 * v0 ** 0.57 * D ** 0.67 / (nu ** 0.24 * b ** 0.43)
    return k


# ---------------------------------------------------------------------------
# Raw sheet loading
# ---------------------------------------------------------------------------
@dataclass
class RawExperiment:
    label: str
    ts: pd.DataFrame           # time, mass, pressure, ret_temp, ret_cond, perm_temp, perm_cond, vial_swap
    vials: pd.DataFrame        # label, cond, icp_mgL
    n_declared: int
    notes: str


def load_raw_sheet(path, sheet_name) -> RawExperiment:
    """Load one raw diafiltration sheet: the time-series (row 3 header, data
    from row 4) and the per-vial ICP/conductivity block (cols O:T, rows
    4-17).

    Uses `Worksheet.iter_rows()` (sequential access) rather than per-cell
    `ws.cell(row=r, column=c)` calls: openpyxl's read-only mode is
    optimized for sequential iteration, and per-cell random access can be
    catastrophically slow on some workbooks (observed: >6 minutes on a
    ~2000-row copolymer sheet that `iter_rows` reads in under a second --
    see docs/prompts/copolymer_nacl_na2so4.md). Output is unchanged (same
    column mapping, same break condition on a blank column-A cell) --
    verified byte-identical against the pre-optimization version on the
    existing NF270 raw sheets used throughout this project."""
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet_name]

    n_declared = ws.cell(row=1, column=2).value
    notes = ws.cell(row=1, column=5).value or ""

    rows = []
    for row in ws.iter_rows(min_row=4, min_col=1, max_col=8, values_only=True):
        if row[0] is None:
            break
        rows.append(list(row))
    ts = pd.DataFrame(
        rows,
        columns=["time_s", "mass_g", "pressure_psi", "ret_temp", "ret_cond",
                 "perm_temp", "perm_cond", "vial_swap"],
    ).astype(
        {"time_s": float, "mass_g": float, "pressure_psi": float,
         "ret_temp": float, "ret_cond": float, "perm_temp": float,
         "perm_cond": float, "vial_swap": "Int64"}
    )

    vial_rows = []
    for row in ws.iter_rows(min_row=4, max_row=17, min_col=15, max_col=20, values_only=True):
        label = row[0]        # col O
        if label is None:
            continue
        cond = row[1]          # col P
        sample_vol = row[2]    # col Q, mL
        acid_vol = row[3]      # col R, mL
        icp_mgL = row[5]       # col T (DILUTED reading)
        vial_rows.append([label, cond, sample_vol, acid_vol, icp_mgL])
    vials = pd.DataFrame(
        vial_rows, columns=["label", "cond", "sample_vol_mL", "acid_vol_mL", "icp_mgL"]
    )

    wb.close()
    return RawExperiment(label=sheet_name, ts=ts, vials=vials,
                          n_declared=int(n_declared), notes=str(notes))


# ---------------------------------------------------------------------------
# (a) Water flux
# ---------------------------------------------------------------------------
def compute_jw(time_s, mass_g, window=DEFAULT_JW_WINDOW, rho=RHO_G_PER_L, A_m=A_M_M2):
    """Smoothed water flux Jw [m/s] from a centered rolling-window linear
    regression slope of mass(t) (dm/dt), converted via Jw = dm/dt/(rho*A_m).

    rho is given in g/L; converted to g/m^3 (1 m^3 = 1000 L) so that
    Jw = (dm/dt [g/s]) / (rho [g/m^3] * A_m [m^2]) comes out in [m/s] ==
    [m^3/m^2/s].
    """
    time_s = np.asarray(time_s, dtype=float)
    mass_g = np.asarray(mass_g, dtype=float)
    n = len(time_s)
    half = window // 2
    dmdt = np.full(n, np.nan)
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        if hi - lo < 2:
            continue
        t_win = time_s[lo:hi]
        m_win = mass_g[lo:hi]
        slope, _ = np.polyfit(t_win, m_win, 1)
        dmdt[i] = slope
    rho_g_per_m3 = rho * 1000.0  # g/L -> g/m^3
    Jw = dmdt / (rho_g_per_m3 * A_m)
    return Jw


# ---------------------------------------------------------------------------
# (b) Conductivity calibration
# ---------------------------------------------------------------------------
@dataclass
class Calibration:
    slope: float
    intercept: float
    n_points: int
    r2: float
    source: str  # "fit" or "given"


def fit_calibration(vials: pd.DataFrame, icp_molar_mass=ICP_MOLAR_MASS,
                     min_points=3, cation_stoich=1) -> Calibration | None:
    """Fit a linear calibration cond = intercept + slope*c_mM from an
    experiment's own vial (conductivity, ICP) pairs. Returns None if fewer
    than `min_points` vials have both conductivity and ICP populated (e.g.
    the 05.27.26/05.07.24 raw sheets, which record no per-vial conductivity
    at all -- see README).

    IMPORTANT: "ICP Salt 1 (mg/L)" is the concentration of the DILUTED ICP
    sample, not the original vial -- it must be scaled back up by the
    dilution factor (sample_vol + acid_vol)/sample_vol before converting to
    mM. Skipping this (as an earlier version of this pipeline did) silently
    understates concentration by ~40-200x and produces a wildly wrong,
    physically implausible calibration slope (~12000 instead of ~O(100));
    see README for the full account of this bug and how it was caught
    (comparing against the Feed/Diafiltrate concentrations stated in each
    run's notes, e.g. "1 mM NaCl Feed, 150 mM NaCl Diafiltrate").

    `cation_stoich` (= Salt.cation_stoich, default 1) divides the
    ICP-derived CATION concentration down to SALT concentration -- 1 for
    every Cl- salt (NaCl/CaCl2/LaCl3, exactly 1 cation per formula unit) but
    2 for Na2SO4 (2 Na+ per formula unit; ICP-OES measures elemental Na,
    confirmed against the "1 mM Na2SO4 feed" notes in the copolymer raw
    workbooks -- see docs/prompts/copolymer_nacl_na2so4.md)."""
    df = vials.dropna(subset=["cond", "icp_mgL", "sample_vol_mL", "acid_vol_mL"])
    if len(df) < min_points:
        return None
    dilution = (df["sample_vol_mL"] + df["acid_vol_mL"]) / df["sample_vol_mL"]
    c_mM = (df["icp_mgL"].to_numpy() * dilution.to_numpy()) / icp_molar_mass / cation_stoich
    cond = df["cond"].to_numpy()
    slope, intercept = np.polyfit(c_mM, cond, 1)
    pred = intercept + slope * c_mM
    ss_res = np.sum((cond - pred) ** 2)
    ss_tot = np.sum((cond - cond.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return Calibration(slope=slope, intercept=intercept, n_points=len(df),
                        r2=r2, source="fit")


def apply_calibration(cond, calib: Calibration):
    return (np.asarray(cond, dtype=float) - calib.intercept) / calib.slope


# ---------------------------------------------------------------------------
# (c) Concentration polarization
# ---------------------------------------------------------------------------
def compute_interfacial(c_bulk_ret, c_bulk_perm, Jw, apply_cp=True, k=None):
    """Return (c_int, c_p). If apply_cp, c_int is CP-corrected from the
    retentate bulk concentration via the thin-film model; c_p is taken as
    the (unpolarized) permeate bulk concentration in both cases. If not
    apply_cp, bulk == interfacial on both sides (matches the "(2)" sheet)."""
    c_bulk_ret = np.asarray(c_bulk_ret, dtype=float)
    c_bulk_perm = np.asarray(c_bulk_perm, dtype=float)
    if not apply_cp:
        return c_bulk_ret.copy(), c_bulk_perm.copy()
    if k is None:
        k = mass_transfer_coefficient()
    Jw = np.asarray(Jw, dtype=float)
    c_p = c_bulk_perm.copy()
    c_int = c_p + (c_bulk_ret - c_p) * np.exp(Jw / k)
    return c_int, c_p


def smooth_permeate_side(c_p, n_floor_rows=DEFAULT_HFLOOR_ROWS,
                          floor_value=DEFAULT_HFLOOR_MM,
                          window=DEFAULT_H_SMOOTH_WINDOW):
    """Approximate the permeate-side (H/O) smoothing: floor the first
    `n_floor_rows` (sensor dead-volume before permeate reaches the
    conductivity cell) at `floor_value`, then apply a rolling median. This
    is a heuristic approximation, not an exact reproduction -- see README
    flag 5."""
    c_p = pd.Series(np.asarray(c_p, dtype=float))
    smoothed = c_p.rolling(window, center=True, min_periods=1).median()
    smoothed.iloc[:n_floor_rows] = floor_value
    return smoothed.to_numpy()


# ---------------------------------------------------------------------------
# (d) Derived quantities
# ---------------------------------------------------------------------------
def compute_derived(c_int, c_p, Jw, i_multiplier=1.0):
    """i_multiplier = Salt.i_multiplier (1.0 for NaCl, reproducing I=c_int
    exactly -- 1.0*x is bit-identical to x -- 3.0 for CaCl2, 6.0 for LaCl3)."""
    c_int = np.asarray(c_int, dtype=float)
    c_p = np.asarray(c_p, dtype=float)
    Jw = np.asarray(Jw, dtype=float)
    I = i_multiplier * c_int
    R = 1.0 - c_p / c_int
    B_coeff = Jw * c_p / (c_int - c_p) * 1e6  # um/s
    return I, R, B_coeff


# ---------------------------------------------------------------------------
# Full pipeline
# ---------------------------------------------------------------------------
@dataclass
class ProcessResult:
    df: pd.DataFrame
    calibration: Calibration | None
    flags: list = field(default_factory=list)


def process_experiment(raw: RawExperiment, apply_cp=True, jw_window=DEFAULT_JW_WINDOW,
                        icp_molar_mass=None, calibration: Calibration | None = None,
                        min_calib_points=3, salt: Salt = SALT_NACL) -> ProcessResult:
    """Run the full raw->processed pipeline for one experiment.

    `salt` (default SALT_NACL, for backward compatibility with all existing
    callers) supplies the ICP molar mass (unless overridden via
    `icp_molar_mass`), the valence-aware ionic-strength multiplier, and the
    CP mass-transfer coefficient's solution diffusivity.

    If `calibration` is not supplied, it is fit from the sheet's own vial
    data; if that fit fails (fewer than `min_calib_points` usable vials),
    the caller-visible `flags` list records it and the returned df has NaN
    concentration columns.
    """
    flags = []
    ts = raw.ts
    if icp_molar_mass is None:
        icp_molar_mass = salt.icp_molar_mass

    if calibration is None:
        calibration = fit_calibration(raw.vials, icp_molar_mass=icp_molar_mass,
                                       min_points=min_calib_points,
                                       cation_stoich=salt.cation_stoich)
        if calibration is None:
            flags.append(
                f"{raw.label}: fewer than {min_calib_points} vials with both "
                "conductivity and ICP populated -- cannot fit a calibration; "
                "concentration columns are NaN (provisional/skip)."
            )

    Jw = compute_jw(ts["time_s"], ts["mass_g"], window=jw_window)

    if calibration is not None:
        c_bulk_ret = apply_calibration(ts["ret_cond"], calibration)
        c_bulk_perm = apply_calibration(ts["perm_cond"], calibration)
        c_int, c_p_raw = compute_interfacial(c_bulk_ret, c_bulk_perm, Jw, apply_cp=apply_cp)
        c_p = smooth_permeate_side(c_p_raw)
        # c_int should use the smoothed c_p too, for consistency with the CP
        # formula's use of c_p (matches how H feeds back into K in the
        # reference sheet's live formula).
        if apply_cp:
            k = mass_transfer_coefficient(D=salt.solution_D)
            c_int = c_p + (c_bulk_ret - c_p) * np.exp(np.asarray(Jw) / k)
        I, R, B_coeff = compute_derived(c_int, c_p, Jw, i_multiplier=salt.i_multiplier)
    else:
        c_int = c_p = I = R = B_coeff = np.full(len(ts), np.nan)

    df = pd.DataFrame({
        "c_int_mM": c_int,          # col B
        "ionic_strength_mM": I,     # col D
        "rejection": R,             # col F
        "c_p_mM": c_p,              # col H
        "Jw_m3_m2_s": Jw,           # col J
        "B_um_s": B_coeff,          # col K
        "time_s": ts["time_s"].to_numpy(),   # col M
        "ret_cond_uS_cm": ts["ret_cond"].to_numpy(),  # col N
        "perm_cond_uS_cm": ts["perm_cond"].to_numpy(),  # col O (raw, unsmoothed)
    })
    return ProcessResult(df=df, calibration=calibration, flags=flags)


if __name__ == "__main__":
    # Smoke test: process the (2)-sheet's source raw experiment (non-CP,
    # embedded calibration) and print a short summary.
    raw = load_raw_sheet(DATA_DIR / "NF270_MC5.xlsx", "05.27.26_NaCl")
    given_calib = Calibration(slope=76.685, intercept=-63.706, n_points=0,
                               r2=float("nan"), source="given (embedded formula)")
    result = process_experiment(raw, apply_cp=False, calibration=given_calib)
    print(f"Processed {raw.label}: n={len(result.df)} rows")
    print(result.df.head())
    if result.flags:
        print("FLAGS:", result.flags)
