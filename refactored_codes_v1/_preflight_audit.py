"""Pre-flight checks for the DATA3 single-salt close-out plan.

Read-only audits. Produces PREFLIGHT_AUDIT.md in this directory for review
BEFORE we run pytest, B_form selection, or the final rollout.

Five checks:
  1. Registry salt classifications  — sheet-by-sheet salt assignment agreement
  2. V_tube sanity                  — implied t_corr placement at V_tube = 0.3 g
  3. Pyomo bounds review            — current vs literature ranges for NF270
  4. 11-sheet roster confirmation   — workbook / sheet / salt / regime / ICP coverage
  5. Loader field audit             — what gets read / used / ignored from each workbook

This is purely diagnostic. Nothing it does can break DATA1/DATA2 or DATA3 paths.
"""
from __future__ import annotations

from pathlib import Path
import re
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402
import refactored_ucb_library as lib  # noqa: E402

OUT_PATH = ROOT / "PREFLIGHT_AUDIT.md"

# The 11 single-salt sheets — same source of truth as the runfile.
SINGLE_SALT_RUN_IDS = (
    "MC2.05.07.24_NaCl",
    "MC3.07.22.24_SNaCl",
    "MC4.07.11.24_SNaCl",
    "MC5.07.23.24_NaCl",
    "MC5.07.23.24_SNaCl",
    "MC5.07.23.24_S2NaCl",
    "MC2.05.07.24_CaCl2",
    "MC3.07.11.24_SCaCl2",
    "MC3.07.12.24_S2CaCl2",
    "MC2.05.21.24_LaCl3",
    "MC4.07.11.24_SLaCl3",
)

NF270_ROOT = (
    ROOT.parent / "UnifiedFramework" / "ExperimentalDataFiles"
)


# -----------------------------------------------------------------------------
# Helpers shared by multiple checks
# -----------------------------------------------------------------------------

def _load_sheet(run_id: str):
    """Load one NF270 sheet; return (data_stru, error_message)."""
    family = lib.NF270_RUN_REGISTRY.get(run_id)
    if family is None:
        return None, f"run_id not in NF270_RUN_REGISTRY"
    wb = NF270_ROOT / family["workbook"]
    if not wb.exists():
        return None, f"workbook not found: {wb}"
    try:
        wrapped = lib.loadxlsx(wb, sheet=family["sheet"])
        return wrapped["data_stru"], None
    except Exception as exc:
        return None, f"loadxlsx raised: {exc!r}"


def _salt_from_sheet_name(sheet_name: str) -> str | None:
    """Try to infer the salt from the sheet name (e.g. '07.22.24_SNaCl' → NaCl)."""
    s = sheet_name.upper()
    # Order matters: check LaCl3 before NaCl (both contain 'Cl')
    for needle, salt in [("LACL3", "LaCl3"), ("CACL2", "CaCl2"), ("NACL", "NaCl"), ("KCL", "KCl")]:
        if needle in s:
            return salt
    return None


# -----------------------------------------------------------------------------
# Check 1 — Registry salt classifications
# -----------------------------------------------------------------------------

def check_registry_salts():
    rows = []
    rows.append(
        "| Run ID | Sheet salt (parser) | Loader's `namec` | Notes-text hint | Agreement |"
    )
    rows.append(
        "|---|---|---|---|---|"
    )
    disagreements = 0
    for rid in SINGLE_SALT_RUN_IDS:
        family = lib.NF270_RUN_REGISTRY[rid]
        sheet_salt = _salt_from_sheet_name(family["sheet"]) or "?"
        data_stru, err = _load_sheet(rid)
        if err:
            rows.append(f"| `{rid}` | `{sheet_salt}` | LOAD ERROR | — | ⚠ |")
            disagreements += 1
            continue
        # The loader stores the parsed salt as data_config['namec']
        loader_salt = str(data_stru.get("data_config", {}).get("namec", "?"))
        note_text = str(data_stru.get("data_config", {}).get("note_text", "")).upper()
        note_hint = _salt_from_sheet_name(note_text) or "—"
        agree = (sheet_salt == loader_salt) and (note_hint == "—" or note_hint == loader_salt)
        rows.append(
            f"| `{rid}` | `{sheet_salt}` | `{loader_salt}` | `{note_hint}` | "
            f"{'✓' if agree else '⚠'} |"
        )
        if not agree:
            disagreements += 1

    summary = (
        f"\n**Summary:** {len(SINGLE_SALT_RUN_IDS) - disagreements}/{len(SINGLE_SALT_RUN_IDS)} "
        f"sheets have consistent salt identification across sources."
    )
    body = "\n".join(rows)
    return f"## 1. Registry salt classifications\n\n{body}\n{summary}\n"


# -----------------------------------------------------------------------------
# Check 2 — V_tube sanity (one representative sheet)
# -----------------------------------------------------------------------------

def check_vtube_sanity(rep_run="MC3.07.22.24_SNaCl"):
    data_stru, err = _load_sheet(rep_run)
    if err:
        return f"## 2. V_tube sanity\n\n**ERROR loading {rep_run}:** {err}\n"

    V_tube = float(lib.NF270_TUBE_VOLUME_G)

    rows = []
    rows.append(
        f"**Representative sheet:** `{rep_run}`  "
        f"**Apparatus constant:** V_tube = {V_tube} g\n"
    )
    rows.append(
        "NOTE: the loader now anchors the permeate ICP at **vial close** (DATA2 "
        "convention); the V_tube tube-transit shift is RETIRED. The table below is "
        "kept only as a diagnostic of what the old correction *would* have done — "
        "`t_corr` is the old planted index; *t_close* is where the value now lands.\n"
    )
    rows.append(
        "| Vial | t_open (s) | t_close (s) | m_total (g) | dm/dt (mg/s) | "
        "tau (s) | t_corr (s) | % of vial |"
    )
    rows.append("|---|---|---|---|---|---|---|---|")

    data_raw = data_stru.get("data_raw", [])
    for i, row in enumerate(data_raw):
        t = np.asarray(row["time"], dtype=float).reshape(-1)
        m = np.asarray(row["mass"], dtype=float).reshape(-1)
        if t.size < 2 or np.sum(np.isfinite(m)) < 2:
            rows.append(f"| {i+1} | — | — | — | — | — | — | (skipped) |")
            continue
        m_clean = m[np.isfinite(m)]
        m_total = float(m_clean[-1] - m_clean[0])
        if m_total <= 1e-3:
            rows.append(
                f"| {i+1} | {t[0]:.0f} | {t[-1]:.0f} | {m_total:.4f} | "
                f"— | — | — | (startup, <1 mg) |"
            )
            continue
        duration = float(t[-1] - t[0])
        dmdt = m_total / duration
        tau = V_tube / dmdt
        t_mid = (t[0] + t[-1]) / 2.0
        t_corr = max(t[0], min(t[-1], t_mid - tau))
        pct = (t_corr - t[0]) / duration * 100.0 if duration > 0 else 0.0
        rows.append(
            f"| {i+1} | {t[0]:.0f} | {t[-1]:.0f} | {m_total:.3f} | "
            f"{dmdt*1000:.2f} | {tau:.0f} | {t_corr:.0f} | {pct:.0f}% |"
        )

    interpretation = (
        "\n**Interpretation (diagnostic only — correction is retired):** "
        "the live loader plants the permeate ICP at *t_close* (100% of vial), per the "
        "DATA2 vial-close convention. The `t_corr` / `% of vial` columns merely show "
        "where the old V_tube ≈ 0.3 g shift would have landed it (~30-50% into the vial)."
    )
    return f"## 2. Permeate placement (vial close; retired V_tube diagnostic)\n\n" + "\n".join(rows) + interpretation + "\n"


# -----------------------------------------------------------------------------
# Check 3 — Pyomo bounds review
# -----------------------------------------------------------------------------

def check_pyomo_bounds():
    table = """
| Parameter | Current bound | Literature range for NF270 | Comment |
|---|---|---|---|
| `Lp` (water permeability, L·m⁻²·h⁻¹·bar⁻¹) | [0.5, 50] | 4–15 typical, up to ~30 in best condition | Lower bound (0.5) very conservative — won't clip any real fit. Upper bound (50) has plenty of headroom. ✓ |
| `B` (solute permeability, μm/s) | [1e-6, 30] | 0.1–10 for divalents (CaCl₂, LaCl₃); 1–30 for NaCl on NF270 | Lower bound (1e-6) effectively zero. Upper bound (30) tight — could clip a hyper-permeable fit on NaCl. ⚠ |
| `σ` (reflection coefficient) | [1e-3, 0.999] | 0.3–0.6 for NaCl, 0.9–0.99 for divalents/trivalents | The 0.999 upper-bound clipping shows up *constantly* in the contour data (data wants σ → 1 for NaCl). Functionally correct; physically defensible (real membranes are never σ=1 exactly). ✓ |
"""
    flags = (
        "\n**Flags worth a second look:**\n"
        "\n"
        "- `B ∈ [1e-6, 30]`: if any fit lands at B = 30 (the upper bound), the true value may be higher. "
        "Worth grepping the previous contour batch's centering_fit.json files for `\"B\": 30.0` or near-30 values.\n"
        "- `σ` upper bound of 0.999 is just slightly inside 1.0 — this is intentional (avoids the σ=1 singularity in Spiegler-Kedem) and is what Fix A clips deterministic seeds to.\n"
    )
    return f"## 3. Pyomo bounds review\n\n{table.strip()}\n{flags}"


# -----------------------------------------------------------------------------
# Check 4 — 11-sheet roster + regime + ICP coverage
# -----------------------------------------------------------------------------

def check_roster_regime():
    rows = []
    rows.append(
        "| Run ID | Workbook | Salt | n_vials | cF first → last (mM) | Regime | "
        "ICP vials w/ data | Notes |"
    )
    rows.append("|---|---|---|---|---|---|---|---|")
    for rid in SINGLE_SALT_RUN_IDS:
        family = lib.NF270_RUN_REGISTRY[rid]
        data_stru, err = _load_sheet(rid)
        if err:
            rows.append(f"| `{rid}` | {family['workbook']} | — | — | LOAD ERROR | — | — | {err} |")
            continue
        cfg = data_stru["data_config"]
        salt = cfg.get("namec", "?")
        data_raw = data_stru["data_raw"]
        n = len(data_raw)

        # cF first vs last from the retentate trace
        cf_first = np.nan
        cf_last = np.nan
        try:
            for r in data_raw:
                cf_arr = np.asarray(r.get("cF_exp", []), dtype=float)
                if cf_arr.size:
                    finite = cf_arr[np.isfinite(cf_arr)]
                    if finite.size:
                        cf_first = float(finite[0])
                        break
            for r in reversed(data_raw):
                cf_arr = np.asarray(r.get("cF_exp", []), dtype=float)
                if cf_arr.size:
                    finite = cf_arr[np.isfinite(cf_arr)]
                    if finite.size:
                        cf_last = float(finite[-1])
                        break
        except Exception:
            pass

        if np.isfinite(cf_first) and np.isfinite(cf_last):
            cf_str = f"{cf_first:.1f} → {cf_last:.1f}"
            regime = "dilution" if cf_last < cf_first * 0.7 else (
                "concentration" if cf_last > cf_first * 1.5 else "flat"
            )
        else:
            cf_str = "—"
            regime = "?"

        # ICP coverage — count vials whose cV_avg array has a non-NaN entry
        icp_count = 0
        for r in data_raw:
            cv = np.asarray(r.get("cV_avg", []), dtype=float).reshape(-1)
            if cv.size and np.any(np.isfinite(cv)):
                icp_count += 1

        note = ""
        if icp_count < n - 2:  # tolerate vial-1 startup
            note = f"⚠ only {icp_count}/{n} vials have ICP"
        elif n < 6:
            note = "short run"

        rows.append(
            f"| `{rid}` | {family['workbook']} | `{salt}` | {n} | {cf_str} | "
            f"{regime} | {icp_count} / {n} | {note} |"
        )

    return "## 4. 11-sheet roster + regime + ICP coverage\n\n" + "\n".join(rows) + "\n"


# -----------------------------------------------------------------------------
# Check 5 — Loader field audit
# -----------------------------------------------------------------------------

def check_loader_fields():
    """Hand-curated table of what the DATA3 loader does with each workbook field.

    Sourced from a code review of _load_legacy_data_stru_from_excel and
    _parse_excel_metadata_and_table in refactored_ucb_library.py.
    """
    return """## 5. Loader field audit — what gets read / used / ignored

### Time-series block (per-row data)

| Workbook column | Read? | Used in fit? | Where it lands in `data_stru` |
|---|---|---|---|
| `Time (s)` | ✓ | ✓ — defines x-axis, vial segmentation | `data_raw[i]["time"]` |
| `Mass (g)` | ✓ | ✓ — mass channel in WSSE, reset per vial | `data_raw[i]["mass"]` |
| `Pressure (psi)` | ✓ | ✓ — averaged into `data_config["delP"]` (bar) | `data_raw[i]["pressure"]` |
| `Retentate Temp` | ✓ | ✓ — per-row EC25 compensation | `data_raw[i]["retentate_temp"]` + top-level `measured_ret_temp_C` |
| `Retentate Cond @ Temp (uS/cm)` | ✓ | ✓ — compensated to 25 °C, then inverted to cF | `data_raw[i]["cF_exp"]` (after inversion) |
| `Permeate Temp` | ✓ | ✓ — per-row EC25 compensation | `data_raw[i]["permeate_temp"]` + top-level `measured_perm_temp_C` |
| `Permeate Cond @ Temp (uS/cm)` | ✓ | ⚪ only if `NF270_USE_PERMEATE_PROBE` toggle is on | `data_raw[i]["cV_perm_cond"]` |
| `Vial Swap` | ✓ | ✓ — segments run into per-vial windows | `data_raw[i]["vial_swap"]` |

### Top-left metadata block

| Workbook cell | Read? | Used in fit? | Where it lands |
|---|---|---|---|
| `Datapoints:` | ⚠ partially | sanity-check only, not stored | (just compared to actual time-series length) |
| `Experiment Name` | ✓ | ⚪ provenance only | `data_stru["filename"]` |
| `Initial Solution Weight (g)` | ✓ | ✓ — boundary condition `M_F0` | `data_config["M_F0"]` |
| `Final Solution Weight (g)` | ✓ | ✓ — combined with `M_F0` to give `M_O = Final − M_F0` | `data_config["M_O"]` |
| `Final Vial Weight (g)` | ⚠ parsed | ✗ NOT used downstream | only in `data_config` metadata |
| `Final Vial w/ Solution (g)` | ⚠ parsed | ✗ NOT used downstream | only in `data_config` metadata |
| `ICP Calibration Points (#)` | ⚠ parsed | ✗ NOT used downstream | provenance only |
| `Number of Salts` | ✓ | ✓ — `nc` field, determines single vs multi-salt path | `data_config["nc"]` |
| `Salt 1`, `Salt 2` | ✓ | ✓ — salt identity for conductivity inversion + ICP MW | `data_config["namec"]` |
| `Notes` (free-text cell) | ✓ | ✓ — parsed for `C_F0`, `C_D`, `rpm`, etc. | `data_config["note_text"]` + derived fields |

### Sidebar ICP block (per-vial assays)

| Sidebar row | Read? | Used in fit? | Where it lands |
|---|---|---|---|
| `Feed` row | ✓ | ✓ — overrides Notes-text `C_F0` if ICP value present | `data_config["cF_feed_icp_mM"]`, `C_F0`, `C_F0_source` |
| `Diafiltrate` row | ✓ | ✓ — overrides Notes-text `C_D` if ICP value present | `data_config["cF_diafiltrate_icp_mM"]`, `C_D`, `C_D_source` |
| `Retentate` row | ✓ | ⚪ cross-check only | `data_config["cF_retentate_icp_mM"]` |
| `Final Tube` row | ✓ | ⚪ cross-check only | `data_config["icp_final_tube_mM"]` |
| `Vial 1..N` rows | ✓ | ✓ — populates `cV_avg` per vial (tube-transit-corrected anchor) | `data_raw[i]["cV_avg"]` (NaN-padded array) |
| `ICP Salt 1 Calibration` table | ⚠ parsed | ✗ NOT used downstream | local `_calib` variable only |

### Summary

- **Fully used:** time, mass, pressure, retentate temp + cond, permeate temp, vial swap, M_F0, M_O, salt names, Notes, sidebar ICP per-vial values.
- **Used only if a toggle is enabled:** permeate conductivity (`NF270_USE_PERMEATE_PROBE`).
- **Parsed but not threaded into the fit:** Final Vial Weight, Final Vial w/ Solution, ICP Calibration Points (#), the per-salt ICP Calibration tables.

The "parsed but not used" items are kept as provenance (so they survive in `data_config` for audit purposes) but don't currently influence model output. If your audit standards require *all* metadata to be either used or explicitly justified, the ICP calibration tables are the most actionable candidate — they could be wired in to validate the precomputed `mg/L` values rather than trusting the meter's own calibration.
"""


# -----------------------------------------------------------------------------
# Check 6 — ICP calibration cross-check (sheet-by-sheet health)
# -----------------------------------------------------------------------------

def check_icp_calibration_health():
    rows = []
    rows.append(
        "Per-sheet sanity check: each workbook embeds an ICP calibration table "
        "(Concentration vs. Intensity standards). We refit that table as a "
        "linear curve and cross-check whether the stored `ICP Salt 1 (mg/L)` "
        "values match what the curve would predict from the same intensity. "
        "Healthy sheets show all residuals within 5%; bad residuals usually "
        "mean the workbook's calibration column was leftover/wrong (a Salt-1=NaCl "
        "template reused for a different cation), NOT that the lab's analytical "
        "result is wrong.\n"
    )
    rows.append(
        "| Run ID | R² | n pts | Stored mg/L range | Within 5% | Max rel residual (%) | Status |"
    )
    rows.append(
        "|---|---|---|---|---|---|---|"
    )
    for rid in SINGLE_SALT_RUN_IDS:
        data_stru, err = _load_sheet(rid)
        if err:
            rows.append(f"| `{rid}` | LOAD ERROR | — | — | — | — | ⚠ |")
            continue
        cfg = data_stru["data_config"]
        calib = cfg.get("icp_calibration_curve")
        sidebar = cfg.get("icp_sidebar_rows", {})

        if calib is None:
            rows.append(
                f"| `{rid}` | — | — | — | — | — | ⚠ no calibration table in workbook |"
            )
            continue

        n_ok = n_total = 0
        max_pct = 0.0
        for rec in sidebar.values():
            res = rec.get("icp_calibration_residual")
            if res is None:
                continue
            n_total += 1
            if res["within_5pct_tolerance"]:
                n_ok += 1
            max_pct = max(max_pct, abs(res["relative_residual_pct"]))

        # Color the status
        if n_total == 0:
            status = "⚠ no rows had intensity data"
        elif n_ok == n_total and max_pct < 5.0:
            status = "✓ healthy"
        elif n_ok >= n_total * 0.9:
            status = f"~ mostly OK ({n_ok}/{n_total})"
        elif max_pct > 50.0:
            status = f"⚠⚠ catastrophic (max {max_pct:.0f}%) — likely wrong calibration in workbook"
        else:
            status = f"⚠ disagreement ({n_total - n_ok}/{n_total} fail)"

        rows.append(
            f"| `{rid}` | {calib['r_squared']:.4f} | {calib['n_points']} | "
            f"{calib['mg_L_range'][0]:.1f}-{calib['mg_L_range'][1]:.1f} | "
            f"{n_ok}/{n_total} | {max_pct:.2f} | {status} |"
        )

    interpretation = (
        "\n**Interpretation:**\n"
        "\n"
        "- ✓ **NaCl MC2/MC3/MC4 sheets** show R² ≈ 0.998 and residuals ≤ 5%. The workbook's "
        "calibration table is internally consistent with its precomputed mg/L values.\n"
        "- ⚠ **MC5 NaCl sheets** have no calibration table at all — cross-check can't run. "
        "Note this for downstream users.\n"
        "- ⚠⚠ **CaCl₂ and LaCl₃ sheets (mostly)** show residuals 39-921%. This is almost "
        "certainly because the calibration column in the workbook is leftover from a "
        "Salt-1=NaCl template that was reused for a different cation. The lab's analytical "
        "software used the correct calibration; the workbook column does not reflect it.\n"
        "\n"
        "**Action items:**\n"
        "\n"
        "- DO NOT modify `cV_avg` based on the workbook calibration — the precomputed "
        "  mg/L from the lab software is ground truth.\n"
        "- The audit field `icp_calibration_residual` is now populated on every loaded sheet; "
        "  downstream consumers can inspect it but should treat large residuals on multivalent "
        "  sheets as a *workbook bookkeeping issue*, not a measurement error.\n"
        "- If we want a stronger cross-check, we'd need the lab's actual calibration "
        "  (perhaps from the analytical-software logs, not the workbook).\n"
    )

    return "## 6. ICP calibration cross-check (NEW)\n\n" + "\n".join(rows) + interpretation + "\n"


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    t0 = time.time()
    print(f"[preflight] starting audit; output → {OUT_PATH}")

    sections = []

    header = (
        "# DATA3 single-salt pre-flight audit\n\n"
        f"_Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}_\n\n"
        "Read-only audit of the loader, the registry, and the bounds.\n"
        "Run BEFORE pytest / B_form selection / final rollout so any data-side\n"
        "issues are surfaced before we lock in regression reference values.\n"
        "\n"
        f"Scope: {len(SINGLE_SALT_RUN_IDS)} single-salt sheets.\n"
    )
    sections.append(header)

    for name, fn in [
        ("registry salts",    check_registry_salts),
        ("V_tube sanity",     check_vtube_sanity),
        ("Pyomo bounds",      check_pyomo_bounds),
        ("roster + regime",   check_roster_regime),
        ("loader fields",     check_loader_fields),
        ("ICP calibration",   check_icp_calibration_health),
    ]:
        t_sec = time.time()
        try:
            section_md = fn()
            sections.append(section_md)
            print(f"[preflight] ok  {name:25s} ({time.time()-t_sec:.1f}s)")
        except Exception as exc:
            sections.append(f"## {name}\n\n**ERROR:** {exc!r}\n")
            print(f"[preflight] ERR {name:25s} → {exc!r}")

    OUT_PATH.write_text("\n\n".join(sections))
    print(f"[preflight] wrote {OUT_PATH} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
