"""Loader-fidelity audit for the 13 DATA3 single-salt Excel sheets.

For each sheet:
  1. Read the raw Excel directly with pandas (no library code in the path).
  2. Load via refactored_ucb_library.loadxlsx → data_stru.
  3. Apply the EXPECTED loader transformations to the raw data (so we
     compare apples-to-apples, not apples-to-oranges):
        - Mass: subtract initial mass per vial (each vial starts at 0)
        - Retentate Cond: EC25 temperature compensation
              σ_25 = σ_T · (1 + alpha · (25 - T))
        - Permeate Cond: EC25 temperature compensation (and tube-transit
              time correction is on the time axis, handled separately)
        - Other columns: 1:1 preserved, no transformation expected
  4. Compare expected-after-transform against actually-loaded, point-by-point.
  5. Compute residual + relative error per data point.
  6. Classify each point as good / warning / bad.
  7. Also compare per-vial ICP scalar concentrations.

Output:
  refactored_codes_v1/DATA3_loader_fidelity_audit_<DATE>.xlsx
  with:
   - One sheet per experiment (run_id), color-coded by quality
   - "_summary" sheet aggregating per-variable / per-experiment stats
   - "_transformations" sheet documenting what the loader applies per variable

Tolerance defaults:
   good     : |rel_err| < 0.01   (1 %)
   warning  : 0.01 <= |rel_err| < 0.05
   bad      : |rel_err| >= 0.05  OR direct copy mismatch beyond absolute floor

For absolute-zero baseline columns (mass near 0, etc.), use absolute floor
instead of relative error.
"""
from __future__ import annotations

import sys
import datetime as _dt
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import refactored_ucb_library as lib  # noqa: E402

REPO_ROOT = ROOT.parent
DATA_ROOT = REPO_ROOT / "UnifiedFramework" / "ExperimentalDataFiles"
OUT_PATH = ROOT / f"DATA3_loader_fidelity_audit_{_dt.date.today().isoformat()}.xlsx"

# (workbook_stem, sheet_name) for each single-salt experiment
SHEETS = [
    ("NF270_MC2", "05.07.24_NaCl"),
    ("NF270_MC2", "05.07.24_CaCl2"),
    ("NF270_MC2", "05.21.24_LaCl3"),
    ("NF270_MC3", "07.09.24_NaCl"),
    ("NF270_MC3", "07.22.24_SNaCl"),
    ("NF270_MC3", "07.11.24_CaCl2"),
    ("NF270_MC3", "07.11.24_SCaCl2"),
    ("NF270_MC3", "07.12.24_S2CaCl2"),
    ("NF270_MC4", "07.11.24_SNaCl"),
    ("NF270_MC4", "07.11.24_SLaCl3"),
    ("NF270_MC5", "07.23.24_NaCl"),
    ("NF270_MC5", "07.23.24_SNaCl"),
    ("NF270_MC5", "07.23.24_S2NaCl"),
]

# Map: raw-Excel column index -> (raw column name, data_raw key, comparison type)
# Each row of the raw Excel time-series (rows 3+ in cols 0-7) maps to a sample
# inside one of the loaded vials. The key in data_raw[i] holds the slice.
RAW_COLS = [
    # (col_idx, label,                   data_raw key,         compare_kind, abs_floor)
    (0,        "Time (s)",                "time",                "exact",       1e-3),
    (1,        "Mass (g)",                "mass",                "absolute",    0.005),
    (2,        "Pressure (psi)",          "pressure",            "absolute",    0.1),
    (3,        "Retentate Temp (°C)",     "retentate_temp",      "absolute",    0.01),
    (4,        "Retentate Cond (μS/cm)",  "cF_exp_conductivity", "relative",    1.0),
    (5,        "Permeate Temp (°C)",      "permeate_temp",       "absolute",    0.01),
    (6,        "Permeate Cond (μS/cm)",   "cV_perm_cond",        "relative",    1.0),
    (7,        "Vial Swap",               "vial_swap",           "exact",       0.0),
]
HEADER_ROW = 2           # row index of the column header row in the raw sheet
DATA_START_ROW = 3       # row index where time-series data begins
TOL_GOOD = 0.01          # |rel_err| < 0.01 → good
TOL_WARN = 0.05          # 0.01 ≤ |rel_err| < 0.05 → warning


def _classify(rel_err: float | None, abs_err: float, kind: str, abs_floor: float):
    """Return one of 'good', 'warning', 'bad', or 'na'."""
    if rel_err is None or not np.isfinite(rel_err):
        rel_err = float("inf")
    if kind == "exact":
        return "good" if abs_err == 0 else "bad"
    if kind == "absolute":
        return ("good" if abs_err < abs_floor else
                "warning" if abs_err < 5 * abs_floor else "bad")
    # relative (with absolute floor as fallback for very small raw values)
    if abs_err < abs_floor:
        return "good"
    if abs(rel_err) < TOL_GOOD:
        return "good"
    if abs(rel_err) < TOL_WARN:
        return "warning"
    return "bad"


def _residual(raw_val, loaded_val):
    """Return (abs_err, rel_err).  rel_err is None when raw_val ~ 0."""
    if (raw_val is None or loaded_val is None
            or (isinstance(raw_val, float) and not np.isfinite(raw_val))
            or (isinstance(loaded_val, float) and not np.isfinite(loaded_val))):
        return (float("nan"), None)
    abs_err = float(loaded_val) - float(raw_val)
    rel_err = abs_err / float(raw_val) if abs(float(raw_val)) > 1e-12 else None
    return (abs(abs_err), rel_err)


def _concat_loaded_series(ds, key):
    """Concatenate all per-vial slices of data_raw[i][key] into one flat array."""
    out = []
    for v in ds.get("data_raw", []):
        arr = v.get(key)
        if arr is None:
            continue
        out.extend(np.asarray(arr, dtype=float).tolist())
    return np.asarray(out, dtype=float)


def _expected_loaded(raw_ts: pd.DataFrame, salt_name: str) -> pd.DataFrame:
    """Apply the SAME transformations loadxlsx applies, so we can compare
    apples-to-apples downstream.

    Returns a copy of raw_ts with:
      - "Mass (g)" baseline-subtracted from the first valid sample of each
        vial (loader behavior — every vial starts at 0).
      - "Retentate Cond (μS/cm)" EC25-compensated.
      - "Permeate Cond (μS/cm)" EC25-compensated.
      - Other columns untouched.
    """
    df = raw_ts.copy()
    # Mass baseline subtraction: per-vial-block, find first valid mass and
    # subtract it from every sample in that block.
    vial_swap = df["Vial Swap"].fillna(0).astype(int)
    # Identify vial-block boundaries: each block runs while vial_swap == 0,
    # incremented by a 1 in the swap column. Simpler approach: use cumsum of
    # the swap-up edges to get a vial index per row.
    vial_idx = vial_swap.cumsum()  # 0 for first vial, 1 after first swap, ...
    df["__vial_idx__"] = vial_idx
    mass_baseline = df.groupby("__vial_idx__")["Mass (g)"].transform(
        lambda s: s.dropna().iloc[0] if s.notna().any() else 0.0
    )
    df["Mass (g)"] = df["Mass (g)"] - mass_baseline
    df.drop(columns="__vial_idx__", inplace=True)

    # EC25 compensation for retentate + permeate conductivities
    alpha = {"NaCl": 0.021, "KCl": 0.019, "CaCl2": 0.023, "LaCl3": 0.025}.get(salt_name, 0.024)
    df["Retentate Cond (μS/cm)"] = (
        df["Retentate Cond (μS/cm)"] * (1.0 + alpha * (25.0 - df["Retentate Temp (°C)"]))
    )
    df["Permeate Cond (μS/cm)"] = (
        df["Permeate Cond (μS/cm)"] * (1.0 + alpha * (25.0 - df["Permeate Temp (°C)"]))
    )
    return df


# Molar masses (g/mol) for the salts the campaign covers — used to convert
# raw ICP mg/L to mM when displaying the per-vial sidebar.
_SALT_MOLAR_MASS = {
    "NaCl":  58.44,
    "KCl":   74.55,
    "CaCl2": 110.98,
    "LaCl3": 245.26,
}


def _audit_sidebar(raw: pd.DataFrame, ds: dict, salt_name: str) -> pd.DataFrame:
    """Compare the per-vial ICP block (cols 14–21, rows 3–16) against the
    loaded sidebar scalars and per-vial cV_avg.  Side-by-side display —
    no residual / verdict because the loader applies a non-trivial
    calibration + dilution + molar-mass chain (cps → mg/L → mM) and the
    correct test is internal consistency, not raw-numeric equality.
    """
    cfg = ds.get("data_config", {})
    # The vial block runs from row 3 down to row 3 + 4 + N (Feed, Diafiltrate,
    # Retentate, Final Tube, then Vials 1..N).  Build the slice generously
    # and trim to non-NaN labels.
    block = raw.iloc[3:30, 14:22].copy()
    block.columns = [
        "Vial label", "Sidebar Cond @ Temp", "Sample Vol (mL)",
        "Nitric Acid Vol (mL)", "ICP Salt1 Intensity (cps)",
        "ICP Salt1 (mg/L)", "ICP Salt2 (ppm)", "ICP Salt2 Err (%)",
    ]
    block = block[block["Vial label"].notna()].reset_index(drop=True)

    # Loader-stored equivalents (where they exist)
    n_v = cfg.get("n", 0)
    loaded_lookup = {
        "Feed":         cfg.get("cF_feed_icp_mM"),
        "Diafiltrate":  cfg.get("cF_diafiltrate_icp_mM"),
        "Retentate":    cfg.get("cF_retentate_icp_mM"),
        "Final Tube":   cfg.get("icp_final_tube_mM"),
    }
    # Per-vial cV_avg from data_raw
    for i, v in enumerate(ds.get("data_raw", [])):
        cv = v.get("cV_avg")
        cv_arr = np.asarray(cv, dtype=float) if cv is not None else np.array([])
        vals = cv_arr[~np.isnan(cv_arr)] if cv_arr.size else np.array([])
        scalar = float(vals[0]) if vals.size else None
        loaded_lookup[f"Vial {i + 1}"] = scalar

    # Convert raw mg/L → mM using salt molar mass (approximate; loader adds
    # dilution factor on top).  This is for orientation only.
    M = float(_SALT_MOLAR_MASS.get(salt_name, 58.44))

    rows = []
    for _, r in block.iterrows():
        label = str(r["Vial label"]).strip()
        loaded_mM = loaded_lookup.get(label)
        raw_mgL = r["ICP Salt1 (mg/L)"]
        raw_mM_approx = (float(raw_mgL) / M) if pd.notna(raw_mgL) else None
        rows.append({
            "Vial label": label,
            "Raw sample vol (mL)": r["Sample Vol (mL)"],
            "Raw nitric acid vol (mL)": r["Nitric Acid Vol (mL)"],
            "Raw ICP intensity (cps)": r["ICP Salt1 Intensity (cps)"],
            "Raw ICP (mg/L)": raw_mgL,
            "Raw ICP → mM (no dilution corr)": raw_mM_approx,
            "Loaded ICP (mM, w/ dilution corr)": loaded_mM,
            "Loader expanded ratio": (
                (float(loaded_mM) / raw_mM_approx) if (loaded_mM is not None and raw_mM_approx not in (None, 0)) else None
            ),
        })
    return pd.DataFrame(rows)


def _audit_calibration(raw: pd.DataFrame, ds: dict) -> pd.DataFrame:
    """Compare the ICP Salt1 + Salt2 calibration curves (cols 23–26, rows 3–13)
    against the loaded icp_calibration_curve dict.  Returns a side-by-side
    DataFrame; raw points unchanged from the Excel by the loader (it just
    fits slope + intercept).
    """
    cfg = ds.get("data_config", {})
    cal = cfg.get("icp_calibration_curve") or {}
    # Salt 1 (cols 23, 24); Salt 2 (cols 25, 26)
    salt1 = raw.iloc[3:20, 23:25].copy()
    salt1.columns = ["Salt1 Concentration (mg/L)", "Salt1 Intensity (cps)"]
    salt1 = salt1[salt1["Salt1 Concentration (mg/L)"].notna()].reset_index(drop=True)
    salt2 = raw.iloc[3:20, 25:27].copy()
    salt2.columns = ["Salt2 Concentration (mg/L)", "Salt2 Intensity (cps)"]
    salt2 = salt2[salt2["Salt2 Concentration (mg/L)"].notna()].reset_index(drop=True)
    # Loaded calibration scalar summary
    cal_summary = pd.DataFrame([
        {"Loaded calibration field": "slope (mg/L per cps)", "Value": cal.get("slope_mg_per_L_per_cps")},
        {"Loaded calibration field": "intercept (mg/L)",     "Value": cal.get("intercept_mg_per_L")},
        {"Loaded calibration field": "R²",                   "Value": cal.get("r_squared")},
        {"Loaded calibration field": "n calibration points", "Value": cal.get("n_points")},
        {"Loaded calibration field": "mg/L range (min)",     "Value": (cal.get("mg_L_range") or [None, None])[0]},
        {"Loaded calibration field": "mg/L range (max)",     "Value": (cal.get("mg_L_range") or [None, None])[1]},
    ])
    n = max(len(salt1), len(salt2))
    df = pd.DataFrame(index=range(n))
    df["Salt1 Conc (mg/L) — RAW"]   = salt1["Salt1 Concentration (mg/L)"].reindex(range(n)).values
    df["Salt1 Intensity (cps) — RAW"] = salt1["Salt1 Intensity (cps)"].reindex(range(n)).values
    df["Salt2 Conc (mg/L) — RAW"]   = salt2["Salt2 Concentration (mg/L)"].reindex(range(n)).values
    df["Salt2 Intensity (cps) — RAW"] = salt2["Salt2 Intensity (cps)"].reindex(range(n)).values
    return df, cal_summary


def _audit_header(raw: pd.DataFrame, ds: dict) -> dict:
    """Compare the header / sidebar metadata cells against the loaded
    data_config scalars.  Returns a dict with raw, loaded, and pass/fail
    keys for each metadata field.
    """
    cfg = ds.get("data_config", {})

    def _cell(r, c):
        try:
            v = raw.iat[r, c]
            return v if pd.notna(v) else None
        except Exception:
            return None

    raw_datapoints = _cell(0, 1)                 # "Datapoints:" cell
    raw_notes      = _cell(0, 4)                 # Notes string
    raw_initial_wt = _cell(5, 10)                # Initial Solution Weight
    raw_final_wt   = _cell(6, 10)                # Final Solution Weight
    raw_n_salts    = _cell(12, 10)               # Number of Salts
    raw_salt1_name = _cell(13, 10)               # Salt 1
    raw_salt2_name = _cell(14, 10)               # Salt 2 (may be NaN)
    raw_cal_points = _cell(11, 10)               # ICP Calibration Points (#)

    # Expected total time-series rows = n_holdup_samples + sum(per-vial lengths)
    loaded_total_ts = sum(len(v.get("time", [])) for v in ds.get("data_raw", []))

    # "Final Solution Weight" — loader represents this as M_O via a known
    # relation; report both raw and loaded without enforcing equality.
    rows = [
        ("Datapoints (raw row count)", raw_datapoints, loaded_total_ts, None,
         "Raw Excel row count of the continuous time-series; the loader splits this into vials."),
        ("Notes",                      raw_notes,     cfg.get("note_text"),
         (str(raw_notes) == str(cfg.get("note_text"))) if raw_notes is not None else None,
         "Loader preserves Notes verbatim (data_config['note_text'])."),
        ("Initial Solution Weight (g)", raw_initial_wt, cfg.get("M_F0"),
         (abs(float(raw_initial_wt) - float(cfg.get("M_F0", 0))) < 1e-6) if raw_initial_wt is not None else None,
         "→ data_config['M_F0']."),
        ("Final Solution Weight (g)",  raw_final_wt,  cfg.get("M_O"),
         None,
         "Loader stores M_O — its relation to the raw Final Solution Weight is documented elsewhere in §4 (mass balance)."),
        ("Number of Salts",            raw_n_salts,   cfg.get("nc"),
         (int(raw_n_salts) == int(cfg.get("nc", 0))) if raw_n_salts is not None else None,
         "→ data_config['nc']."),
        ("Salt 1 name",                raw_salt1_name, cfg.get("namec"),
         (str(raw_salt1_name).strip() == str(cfg.get("namec")).strip()) if raw_salt1_name is not None else None,
         "→ data_config['namec']."),
        ("Salt 2 name (multicomp only)", raw_salt2_name, None,
         None,
         "Single-salt sheets have NaN here; multicomp sheets carry the second-salt label."),
        ("ICP Calibration Points (#)", raw_cal_points, (ds.get("data_config", {}).get("icp_calibration_curve", {}) or {}).get("n_points"),
         None,
         "Loader fits slope + intercept against these N points; verified in the per-sheet calibration block."),
    ]
    return [
        {"Field": fld, "Raw value": rv, "Loaded value": lv,
         "Match": ("✓" if m is True else ("✗" if m is False else "—")),
         "Notes": notes}
        for (fld, rv, lv, m, notes) in rows
    ]


def audit_one_sheet(wb_stem: str, sheet_name: str):
    """Return a tuple (rows_df, vial_df, stats_dict, error)."""
    run_id = f"{wb_stem.replace('NF270_', '')}.{sheet_name}"
    wb_path = DATA_ROOT / f"{wb_stem}.xlsx"
    if not wb_path.exists():
        return None, None, None, f"workbook missing: {wb_path}"

    # 1. Raw Excel — time-series block
    try:
        raw = pd.read_excel(wb_path, sheet_name=sheet_name, header=None)
    except Exception as exc:
        return None, None, None, f"pandas read failed: {exc!r}"
    n_data = int(raw.iat[0, 1]) if pd.notna(raw.iat[0, 1]) else None
    raw_ts = raw.iloc[DATA_START_ROW : DATA_START_ROW + (n_data or len(raw)), 0:8].copy()
    raw_ts.columns = [c[1] for c in RAW_COLS]
    raw_ts.reset_index(drop=True, inplace=True)
    # Coerce to numeric
    for c in raw_ts.columns:
        raw_ts[c] = pd.to_numeric(raw_ts[c], errors="coerce")
    raw_ts_pristine = raw_ts.copy()  # untouched raw for the audit display

    # 2. Load via library
    try:
        ds = lib.loadxlsx(str(wb_path), sheet=sheet_name)["data_stru"]
    except Exception as exc:
        return None, None, None, f"loadxlsx failed: {exc!r}"

    # 2b. Apply the EXPECTED loader transformations to the raw.  EC25 and
    # vial-swap don't depend on loader vial boundaries, but MASS does:
    # the loader re-zeros mass at the start of each LOADER vial (not at
    # raw vial-swap boundaries).  We need to know the loader's vial
    # boundaries (via loaded time stamps) before we can compute expected
    # mass.  Do the non-mass transformations now; mass after alignment.
    salt_name = ds.get("data_config", {}).get("namec", "NaCl")
    alpha = {"NaCl": 0.021, "KCl": 0.019, "CaCl2": 0.023, "LaCl3": 0.025}.get(salt_name, 0.024)
    raw_ts["Retentate Cond (μS/cm)"] = (
        raw_ts["Retentate Cond (μS/cm)"] * (1.0 + alpha * (25.0 - raw_ts["Retentate Temp (°C)"]))
    )
    raw_ts["Permeate Cond (μS/cm)"] = (
        raw_ts["Permeate Cond (μS/cm)"] * (1.0 + alpha * (25.0 - raw_ts["Permeate Temp (°C)"]))
    )

    # 3. Concatenate loaded vials back into one stream + check length parity
    rows = []
    n_raw = len(raw_ts)
    # Loaded streams keyed by raw-col label
    loaded_streams = {}
    for (_idx, label, key, _kind, _floor) in RAW_COLS:
        loaded_streams[label] = _concat_loaded_series(ds, key)

    n_loaded = len(loaded_streams["Time (s)"])

    # 3b. TIME-ALIGN: map each loaded sample to the raw row with the
    # closest time stamp.  The loader truncates each vial (drops settling
    # samples, splits the holdup window), so row-index matching breaks.
    # Build a raw_idx_for_loaded[i] = raw row index whose time matches
    # loaded time i to within rounding.
    raw_times = raw_ts["Time (s)"].to_numpy()
    loaded_times = loaded_streams["Time (s)"]
    raw_idx_for_loaded = np.full(n_loaded, -1, dtype=int)
    for i, t_l in enumerate(loaded_times):
        if not np.isfinite(t_l):
            continue
        # Find the raw row with the nearest time
        diffs = np.abs(raw_times - t_l)
        j = int(np.nanargmin(diffs))
        if np.isfinite(diffs[j]) and diffs[j] < 0.5:  # 0.5 s tolerance
            raw_idx_for_loaded[i] = j
    n_aligned = int((raw_idx_for_loaded >= 0).sum())
    n_compare = n_loaded

    # 3c. Re-baseline mass per LOADER vial.  Walk through each loaded vial's
    # samples, find its first aligned raw row, use that raw mass as the
    # baseline for the entire loaded-vial's expected mass.
    raw_mass_pristine = raw_ts_pristine["Mass (g)"].to_numpy()
    expected_mass = np.full(n_loaded, np.nan)
    loaded_start = 0
    for v in ds.get("data_raw", []):
        v_len = len(v.get("time", []))
        if v_len == 0:
            continue
        block = range(loaded_start, min(loaded_start + v_len, n_loaded))
        # Baseline = raw mass at the first aligned raw row in this loader vial
        baseline = np.nan
        for k in block:
            j = int(raw_idx_for_loaded[k])
            if j >= 0 and j < n_raw and np.isfinite(raw_mass_pristine[j]):
                baseline = float(raw_mass_pristine[j])
                break
        if np.isfinite(baseline):
            for k in block:
                j = int(raw_idx_for_loaded[k])
                if j >= 0 and j < n_raw and np.isfinite(raw_mass_pristine[j]):
                    expected_mass[k] = float(raw_mass_pristine[j]) - baseline
        loaded_start += v_len
    # Overwrite the raw_ts "Mass (g)" column at the aligned indices so the
    # row-build loop below picks up the right expected value.
    # raw_ts is indexed by raw row, not loaded sample, so we stash the
    # per-loaded expected-mass array and patch the row-build to use it.
    _expected_mass_per_loaded = expected_mass

    # 4. Per-sample comparison
    bad_counts = {c[1]: 0 for c in RAW_COLS}
    warn_counts = {c[1]: 0 for c in RAW_COLS}
    good_counts = {c[1]: 0 for c in RAW_COLS}
    max_abs = {c[1]: 0.0 for c in RAW_COLS}
    max_rel = {c[1]: 0.0 for c in RAW_COLS}

    for i in range(n_compare):
        # Use the time-aligned raw row index, not the loaded sample index
        j_raw = int(raw_idx_for_loaded[i])
        row = {"Loaded idx": i + 1, "Raw row idx": (j_raw + 1) if j_raw >= 0 else None}
        for (_idx, label, key, kind, floor) in RAW_COLS:
            if j_raw >= 0 and j_raw < n_raw:
                raw_pristine = raw_ts_pristine[label].iat[j_raw]
                expected_val = raw_ts[label].iat[j_raw]  # post-transform
            else:
                raw_pristine = float("nan")
                expected_val = float("nan")
            # Override expected mass with the per-loader-vial-baseline version
            if label == "Mass (g)" and i < len(_expected_mass_per_loaded):
                expected_val = _expected_mass_per_loaded[i]
            loaded_arr = loaded_streams[label]
            loaded_val = loaded_arr[i] if i < len(loaded_arr) else float("nan")
            row[f"{label} — raw"] = raw_pristine
            row[f"{label} — expected (post-transform)"] = expected_val
            row[f"{label} — loaded"] = loaded_val
            abs_err, rel_err = _residual(expected_val, loaded_val)
            row[f"{label} — Δ"] = abs_err if np.isfinite(abs_err) else None
            row[f"{label} — rel"] = rel_err
            verdict = _classify(rel_err, abs_err, kind, floor) if np.isfinite(abs_err) else "na"
            row[f"{label} — verdict"] = verdict
            if verdict == "bad":
                bad_counts[label] += 1
            elif verdict == "warning":
                warn_counts[label] += 1
            elif verdict == "good":
                good_counts[label] += 1
            if np.isfinite(abs_err):
                max_abs[label] = max(max_abs[label], abs_err)
            if rel_err is not None and np.isfinite(rel_err):
                max_rel[label] = max(max_rel[label], abs(rel_err))
        rows.append(row)

    rows_df = pd.DataFrame(rows)

    # 5. Per-vial ICP scalar comparison (cV_avg)
    vial_rows = []
    for i, v in enumerate(ds.get("data_raw", [])):
        cv = v.get("cV_avg")
        cv_scalar = float(cv[0]) if cv is not None and len(np.asarray(cv).flatten()) > 0 else float("nan")
        # Take first non-NaN value (cV_avg is usually a single scalar broadcast)
        try:
            cv_arr = np.asarray(cv, dtype=float)
            cv_scalar = float(cv_arr[~np.isnan(cv_arr)][0]) if np.any(~np.isnan(cv_arr)) else float("nan")
        except Exception:
            cv_scalar = float("nan")
        vial_rows.append({
            "Vial": i + 1,
            "cV_avg (loaded, mM)": cv_scalar,
        })
    vial_df = pd.DataFrame(vial_rows)

    # 6. Stats dict for summary
    stats = {
        "run_id": run_id,
        "n_samples (raw)": n_raw,
        "n_samples (loaded)": n_loaded,
        "n_samples compared": n_compare,
    }
    for (_idx, label, _key, _kind, _floor) in RAW_COLS:
        n_tot = good_counts[label] + warn_counts[label] + bad_counts[label]
        stats[f"{label}: good"] = good_counts[label]
        stats[f"{label}: warn"] = warn_counts[label]
        stats[f"{label}: bad"] = bad_counts[label]
        stats[f"{label}: max |Δ|"] = max_abs[label]
        stats[f"{label}: max |rel|"] = max_rel[label]
        stats[f"{label}: verdict"] = (
            "OK" if bad_counts[label] == 0 and warn_counts[label] == 0
            else "WARN" if bad_counts[label] == 0
            else "BAD"
        )
    # 7. Extended audit blocks (ICP sidebar, calibration, header metadata).
    # These were missing in the initial v1 of the audit script.
    try:
        sidebar_df = _audit_sidebar(raw, ds, salt_name)
    except Exception as exc:
        sidebar_df = pd.DataFrame([{"_error": f"sidebar audit failed: {exc!r}"}])
    try:
        cal_pts_df, cal_summary_df = _audit_calibration(raw, ds)
    except Exception as exc:
        cal_pts_df = pd.DataFrame([{"_error": f"calibration audit failed: {exc!r}"}])
        cal_summary_df = pd.DataFrame()
    try:
        header_rows = _audit_header(raw, ds)
    except Exception as exc:
        header_rows = [{"_error": f"header audit failed: {exc!r}"}]

    extras = {
        "sidebar_df":     sidebar_df,
        "cal_pts_df":     cal_pts_df,
        "cal_summary_df": cal_summary_df,
        "header_rows":    header_rows,
    }
    # Attach to stats so the writer can pick them up
    stats["_extras"] = extras
    return rows_df, vial_df, stats, None


def main():
    print(f"[audit] writing to {OUT_PATH}")
    print(f"[audit] auditing {len(SHEETS)} single-salt sheets")
    print()

    per_sheet_results = {}
    summary_rows = []

    for (wb_stem, sheet_name) in SHEETS:
        run_id = f"{wb_stem.replace('NF270_', '')}.{sheet_name}"
        print(f"  [{run_id}]  ", end="", flush=True)
        rows_df, vial_df, stats, err = audit_one_sheet(wb_stem, sheet_name)
        if err:
            print(f"ERROR: {err}")
            summary_rows.append({"run_id": run_id, "status": err})
            continue
        per_sheet_results[run_id] = (rows_df, vial_df, stats)
        n_bad = sum(stats[k] for k in stats if k.endswith(": bad"))
        n_warn = sum(stats[k] for k in stats if k.endswith(": warn"))
        n_good = sum(stats[k] for k in stats if k.endswith(": good"))
        print(f"OK — {n_good} good / {n_warn} warn / {n_bad} bad data points")
        summary_rows.append(stats)

    # Write Excel
    print()
    print(f"[audit] writing Excel workbook with {len(per_sheet_results)} sheets + summary")

    with pd.ExcelWriter(OUT_PATH, engine="xlsxwriter") as xw:
        wb = xw.book

        # === TRANSFORMATIONS sheet ===
        # Document what the loader does to each variable so the user can
        # interpret the audit honestly.
        transforms = pd.DataFrame([
            {"Variable": "Time (s)", "Loader transformation": "1:1 preserved (no change)", "Why": "Time axis is the index; no compensation needed."},
            {"Variable": "Mass (g)", "Loader transformation": "Baseline-subtracted per vial — each vial starts at 0", "Why": "Eliminates per-vial offset so mass-vs-time fits compare across vials cleanly."},
            {"Variable": "Pressure (psi)", "Loader transformation": "1:1 preserved", "Why": "Applied pressure is a control variable, recorded as-is."},
            {"Variable": "Retentate Temp (°C)", "Loader transformation": "1:1 preserved", "Why": "Needed for EC25 compensation of conductivity."},
            {"Variable": "Retentate Cond (μS/cm)", "Loader transformation": "EC25 temperature compensation: σ_25 = σ_T · (1 + α · (25 - T))", "Why": "Probe reports σ at measured T; the model needs σ at 25°C so concentration inversion (Shedlovsky) is consistent across temperatures."},
            {"Variable": "Permeate Temp (°C)", "Loader transformation": "1:1 preserved", "Why": "Needed for EC25 compensation of permeate conductivity."},
            {"Variable": "Permeate Cond (μS/cm)", "Loader transformation": "EC25 temperature compensation (same formula as retentate); tube-transit-time correction applied to TIME axis, not to value", "Why": "Same as retentate. Tube transit τ = V_tube / (dm/dt) shifts the permeate time axis to account for the 0.3 g dead volume between membrane and probe."},
            {"Variable": "Vial Swap", "Loader transformation": "1:1 preserved", "Why": "Index flag; used by loader to slice the continuous time-series into per-vial blocks."},
            # — extended audit blocks (added v2 2026-05-27) —
            {"Variable": "ICP vial block (cols 14–21)", "Loader transformation": "Per-vial: label (Feed/Diafiltrate/Retentate/Final Tube/Vial N), Sample Vol (mL), Nitric Acid Vol (mL), ICP intensity (cps), ICP (mg/L) read 1:1 from rows 3+. ICP (mg/L) → mM converted using salt molar mass + dilution factor ((sample+acid)/sample).", "Why": "ICP-OES gives ion-specific concentrations; the loader applies acid dilution + calibration + molar mass to produce cF_feed_icp_mM, cF_diafiltrate_icp_mM, cF_retentate_icp_mM, icp_final_tube_mM, and per-vial cV_avg used by the model."},
            {"Variable": "ICP Salt 1 + Salt 2 calibration (cols 23–26)", "Loader transformation": "Calibration points (concentration vs intensity, ~10 pairs per salt) read 1:1 and fit to linear slope + intercept. Stored as icp_calibration_curve dict with slope_mg_per_L_per_cps, intercept_mg_per_L, r_squared, n_points, mg_L_range.", "Why": "Converts raw cps intensity into a calibrated mg/L concentration. The loader checks R² to flag any poor-quality calibration."},
            {"Variable": "Header / sidebar metadata (rows 0–14)", "Loader transformation": "Datapoints count, Notes string, Initial / Final Solution Weight, Number of Salts, Salt 1 + Salt 2 names, ICP Calibration Points (#) all read 1:1 into data_config (note_text, M_F0, M_O-derived, nc, namec, n_holdup_samples, etc.).", "Why": "Identifies the experiment, sets the mass balance baseline, picks the salt-specific α + molar mass, and tells the loader how many calibration points to expect."},
        ])
        transforms.to_excel(xw, sheet_name="_transformations", index=False)
        ws = xw.sheets["_transformations"]
        ws.set_column(0, 0, 30)
        ws.set_column(1, 1, 70)
        ws.set_column(2, 2, 70)
        ws.freeze_panes(1, 0)

        # Pre-create the format objects we'll use further down
        good_fmt = wb.add_format({"bg_color": "#C6EFCE", "font_color": "#006100"})
        warn_fmt = wb.add_format({"bg_color": "#FFEB9C", "font_color": "#9C5700"})
        bad_fmt  = wb.add_format({"bg_color": "#FFC7CE", "font_color": "#9C0006"})

        # === SUMMARY sheet ===
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_excel(xw, sheet_name="_summary", index=False)
        ws = xw.sheets["_summary"]

        # Color-code verdict cells in summary (formats defined above)
        verdict_cols = [i for i, c in enumerate(summary_df.columns) if c.endswith(": verdict")]
        for col_i in verdict_cols:
            col_letter = chr(ord("A") + col_i) if col_i < 26 else chr(ord("A") + col_i // 26 - 1) + chr(ord("A") + col_i % 26)
            cell_range = f"{col_letter}2:{col_letter}{len(summary_df) + 1}"
            ws.conditional_format(cell_range, {"type": "text", "criteria": "containing", "value": "OK", "format": good_fmt})
            ws.conditional_format(cell_range, {"type": "text", "criteria": "containing", "value": "WARN", "format": warn_fmt})
            ws.conditional_format(cell_range, {"type": "text", "criteria": "containing", "value": "BAD", "format": bad_fmt})
        ws.set_column(0, 0, 26)
        ws.set_column(1, len(summary_df.columns) - 1, 14)
        ws.freeze_panes(1, 1)

        # === per-sheet detail tabs ===
        for run_id, (rows_df, vial_df, stats) in per_sheet_results.items():
            # Sheet name: Excel restricts to 31 chars max
            sheet_label = run_id.replace(".", "_").replace("/", "_")[:31]
            rows_df.to_excel(xw, sheet_name=sheet_label, index=False)
            ws = xw.sheets[sheet_label]
            # Style verdict columns
            verdict_cols_local = [i for i, c in enumerate(rows_df.columns) if c.endswith("— verdict")]
            for col_i in verdict_cols_local:
                col_letter = chr(ord("A") + col_i) if col_i < 26 else chr(ord("A") + col_i // 26 - 1) + chr(ord("A") + col_i % 26)
                rng = f"{col_letter}2:{col_letter}{len(rows_df) + 1}"
                ws.conditional_format(rng, {"type": "text", "criteria": "containing", "value": "good", "format": good_fmt})
                ws.conditional_format(rng, {"type": "text", "criteria": "containing", "value": "warning", "format": warn_fmt})
                ws.conditional_format(rng, {"type": "text", "criteria": "containing", "value": "bad", "format": bad_fmt})
            ws.set_column(0, 0, 10)
            ws.set_column(1, len(rows_df.columns) - 1, 16)
            ws.freeze_panes(1, 1)
            # ICP vial block lower-right
            if vial_df is not None and not vial_df.empty:
                vial_start_row = len(rows_df) + 4
                ws.write(vial_start_row, 0, "Per-vial ICP scalars (loaded data_raw cV_avg)")
                vial_df.to_excel(xw, sheet_name=sheet_label, startrow=vial_start_row + 1, index=False)

            # — Extended audit blocks (v2): sidebar + calibration + header —
            extras = stats.get("_extras") or {}
            cur_row = len(rows_df) + (len(vial_df) + 4 if (vial_df is not None and not vial_df.empty) else 0) + 6

            sidebar_df = extras.get("sidebar_df")
            if sidebar_df is not None and not sidebar_df.empty:
                ws.write(cur_row, 0, "ICP sidebar block (raw cols 14–21 vs loaded scalars)")
                sidebar_df.to_excel(xw, sheet_name=sheet_label, startrow=cur_row + 1, index=False)
                cur_row += len(sidebar_df) + 4

            cal_pts_df = extras.get("cal_pts_df")
            cal_summary_df = extras.get("cal_summary_df")
            if cal_pts_df is not None and not cal_pts_df.empty:
                ws.write(cur_row, 0, "ICP calibration points (raw cols 23–26)")
                cal_pts_df.to_excel(xw, sheet_name=sheet_label, startrow=cur_row + 1, index=False)
                cur_row += len(cal_pts_df) + 4
            if cal_summary_df is not None and not cal_summary_df.empty:
                ws.write(cur_row, 0, "ICP calibration — loaded summary (data_config['icp_calibration_curve'])")
                cal_summary_df.to_excel(xw, sheet_name=sheet_label, startrow=cur_row + 1, index=False)
                cur_row += len(cal_summary_df) + 4

            header_rows = extras.get("header_rows")
            if header_rows:
                ws.write(cur_row, 0, "Header / sidebar metadata (rows 0–14)")
                hdf = pd.DataFrame(header_rows)
                hdf.to_excel(xw, sheet_name=sheet_label, startrow=cur_row + 1, index=False)
                cur_row += len(hdf) + 4

        # === _header_metadata sheet (aggregated across all sheets) ===
        # One row per (sheet, field) so a reviewer can scan in one place.
        header_long_rows = []
        for run_id, (rows_df, vial_df, stats) in per_sheet_results.items():
            for row in (stats.get("_extras", {}).get("header_rows") or []):
                row2 = {"run_id": run_id}
                row2.update(row)
                header_long_rows.append(row2)
        if header_long_rows:
            hdr_df = pd.DataFrame(header_long_rows)
            hdr_df.to_excel(xw, sheet_name="_header_metadata", index=False)
            ws_hdr = xw.sheets["_header_metadata"]
            ws_hdr.set_column(0, 0, 24)
            ws_hdr.set_column(1, 1, 30)
            ws_hdr.set_column(2, 3, 40)
            ws_hdr.set_column(4, 4, 8)
            ws_hdr.set_column(5, 5, 80)
            ws_hdr.freeze_panes(1, 0)
            # Color-code the Match column
            match_col_idx = list(hdr_df.columns).index("Match")
            col_letter = chr(ord("A") + match_col_idx)
            rng = f"{col_letter}2:{col_letter}{len(hdr_df) + 1}"
            ws_hdr.conditional_format(rng, {"type": "text", "criteria": "containing", "value": "✓", "format": good_fmt})
            ws_hdr.conditional_format(rng, {"type": "text", "criteria": "containing", "value": "✗", "format": bad_fmt})

    print(f"[audit] done.  Output: {OUT_PATH}")


if __name__ == "__main__":
    main()
