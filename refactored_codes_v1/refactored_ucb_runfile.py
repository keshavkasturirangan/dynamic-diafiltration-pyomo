#!/usr/bin/env python3
"""
Standalone CLI for the UCB diafiltration refactor.

Tree-shaped workflow
====================
                                           __ Mass plots
                                          /_ Concentration plots
   ROOTS                TRUNK            /__ Parameter tables       BRANCHES
   ─────                ─────────        ____ FIM heatmap           ────────
   DATA1 (.mat) ──┐                     /__ DoE recommendations
   DATA2 (.mat) ──┤    model_construct  ──> render_*
   DATA3 (.xlsx)──┤    +  ParmEst       ──> plot_*
   Custom file ──┘    +  multistart    ──> tables, contours, ...
                       +  FIM
                       +  Pyomo.DoE

   ROOTS  are loaders.        loadmat / loadxlsx  (library)
   TRUNK  is the science.     solve_model, solve_mass_balance_only,
                               estimate_parameters, calc_FIM,
                               compute_doe_metrics, design_next_experiment.
   BRANCHES are outputs.      plot_sim_comparison, render_*_parameters_table,
                               plot_contour, render_concentration_*, etc.

All processing and plotting lives in refactored_ucb_library.py. This script
is just the orchestrator — it walks down the tree, asks four questions
(ROOT, SUBSET, TRUNK, BRANCHES), then dispatches into the library to
produce the chosen artifacts.

The library is self-contained — no patches to apply. The only sibling
dependency is conductivity_paper.py for conductivity-to-concentration
conversion (do NOT modify it).

Environment overrides
---------------------
  DIAFILTRATION_DATA1_ROOT  — DATA1 .mat directory
  DIAFILTRATION_DATA2_ROOT  — DATA2 .mat directory
  DIAFILTRATION_NF270_ROOT  — NF270 XLSX directory (NF270_MC2..MC5.xlsx)
"""
from __future__ import annotations

import os
import sys
import zipfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import refactored_ucb_library as ucb


# =============================================================================
# Constants — data roots and well-known subsets
# =============================================================================

REPO_ROOT = Path(__file__).resolve().parents[1]

DATA1_ROOT = Path(
    os.environ.get("DIAFILTRATION_DATA1_ROOT",
                    REPO_ROOT / "legacy" / "data1_matlab" / "data")
).expanduser().resolve()

DATA2_ROOT = Path(
    os.environ.get("DIAFILTRATION_DATA2_ROOT",
                    REPO_ROOT / "legacy" / "data1_matlab" / "data_library")
).expanduser().resolve()

_NF270_DEFAULT = Path(
    "/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo"
    "/UnifiedFramework/ExperimentalDataFiles"
)
NF270_ROOT = Path(
    os.environ.get("DIAFILTRATION_NF270_ROOT", _NF270_DEFAULT)
).expanduser().resolve()

DATA3_CUSTOM_FIGURES = (
    REPO_ROOT / "UnifiedFramework" / "DATA3" / "figures" / "data3_option3"
)

DATA1_FAST_SUBSET = ("figure_2", "figure_3", "data_analysis")

# DATA2_LEGACY_WHITELIST mirrors the old runfile's DATA2_PIPELINE_ONLY_FIGURES
# tuple exactly. Picking "fast subset" for DATA2 reproduces the old behavior
# byte-for-byte (same materialize_all call, same FigureSpec subset).
# The new default (subset=None) runs the FULL DATA2 manifest, which is a
# strict superset — every figure that used to come out still comes out with
# the same numbers; additional figures from the rest of the manifest are
# produced too. The numeric results for the 4 legacy figures are unchanged.
DATA2_FAST_SUBSET = (
    "calibration_plots",
    "pressure_changes",
    "startup_barplot",
    "model_error_visualization",
)
NF270_FAST_SUBSET = ("MC2.05.08.24_E1", "MC4.07.12.24_E13",
                     "tables.parameters_table")

# Single-salt subset: 11 sheets identified as green-flagged AND single-salt
# (NaCl-only, CaCl2-only, or LaCl3-only). Source of truth:
# nf270_single_salt_filter.json.
NF270_SINGLE_SALT_RUNS = (
    # NaCl (6 sheets)
    "MC2.05.07.24_NaCl",
    "MC3.07.22.24_SNaCl",
    "MC4.07.11.24_SNaCl",
    "MC5.07.23.24_NaCl",
    "MC5.07.23.24_SNaCl",
    "MC5.07.23.24_S2NaCl",
    # CaCl2 (3 sheets)
    "MC2.05.07.24_CaCl2",
    "MC3.07.11.24_SCaCl2",
    "MC3.07.12.24_S2CaCl2",
    # LaCl3 (2 sheets)
    "MC2.05.21.24_LaCl3",
    "MC4.07.11.24_SLaCl3",
)


# =============================================================================
# Small helpers
# =============================================================================

def _prompt(message: str, default: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    value = input(f"{message}{suffix}: ").strip()
    return value or default


def _prompt_yes_no(message: str, default: str = "n") -> bool:
    return _prompt(message + " (y/n)", default).lower().startswith("y")


def _print_outputs(label: str, results: list[dict]) -> None:
    ok = [r for r in results if r.get("status") == "ok"]
    err = [r for r in results if r.get("status") == "error"]
    print(f"\n[{label}] manifest entries attempted: {len(results)} "
          f"(ok={len(ok)}, error={len(err)})")
    if err:
        print("\nFailures:")
        for entry in err:
            print(f"  ! {entry['name']}: {entry['error']}")
    n_paths = sum(len(r.get("paths", [])) for r in ok)
    print(f"\nFiles written: {n_paths}")
    if n_paths and n_paths <= 40:
        for entry in ok:
            for path in entry.get("paths", []):
                print(f"  + {path}")
    elif n_paths:
        print(f"  (too many to list inline; see save_dir)")


def _print_coverage(campaign: str, save_dir: Path) -> None:
    fn = getattr(ucb, "report_paper_coverage", None)
    fmt = getattr(ucb, "format_paper_coverage", None)
    if fn is None or fmt is None:
        return
    try:
        report = fn(campaign, save_dir)
        print(fmt(report))
    except Exception:
        pass


def _ensure_data2_library() -> None:
    sentinel = DATA2_ROOT / "data_stru-dataset270511.123.mat"
    if sentinel.exists():
        return
    archive_path = REPO_ROOT / "legacy" / "data1_matlab" / "inputdata_mat_files.zip"
    if not archive_path.exists():
        raise FileNotFoundError(
            f"DATA2 source folder is missing at {DATA2_ROOT} and archive was "
            f"not found at {archive_path}.")
    DATA2_ROOT.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path, "r") as zip_ref:
        zip_ref.extractall(DATA2_ROOT.parent)


# =============================================================================
# ROOTS — pick a data source
# =============================================================================

def pick_root() -> str:
    print("\n" + "=" * 60)
    print("ROOTS  —  pick a data source")
    print("=" * 60)
    print("  1. DATA1  (.mat — published paper figures)")
    print("  2. DATA2  (.mat — published paper figures)")
    print("  3. DATA3  (.xlsx — NF270 experimental campaign)")
    print("  4. Custom (.mat or .xlsx — single-file run)")
    choice = _prompt("Root", "3").strip().lower()
    if choice in {"1", "data1"}:  return "DATA1"
    if choice in {"2", "data2"}:  return "DATA2"
    if choice in {"3", "data3", "nf270"}:  return "DATA3"
    return "CUSTOM"


def pick_subset(root: str):
    """Within a chosen root, pick which sheets / files to run."""
    if root == "DATA1":
        if _prompt_yes_no("Use fast subset?", "n"):
            return DATA1_FAST_SUBSET
        return None
    if root == "DATA2":
        if _prompt_yes_no("Use fast subset?", "n"):
            return DATA2_FAST_SUBSET
        return None
    if root == "DATA3":
        print("\nSubset within NF270 campaign:")
        print("  1. All 26 green-flagged sheets")
        print("  2. Single-salt only (11 sheets — 6 NaCl, 3 CaCl2, 2 LaCl3)")
        print("  3. One specific sheet (you'll be prompted)")
        print("  4. Fast smoke-test (3 entries)")
        choice = _prompt("Subset", "2").strip()
        if choice == "1":
            return None
        if choice == "3":
            return _prompt_one_nf270_run()
        if choice == "4":
            return NF270_FAST_SUBSET
        return NF270_SINGLE_SALT_RUNS
    return None


def _prompt_one_nf270_run() -> tuple[str, ...] | None:
    if not hasattr(ucb, "NF270_RUN_REGISTRY"):
        return None
    runs = list(ucb.NF270_RUN_REGISTRY.keys())
    print(f"\n  available {len(runs)} runs (first 10 shown):")
    for run_id in runs[:10]:
        print(f"    {run_id}")
    if len(runs) > 10:
        print(f"    ...and {len(runs)-10} more")
    pick = _prompt("Enter exact run_id", runs[0]).strip()
    if pick in runs:
        return (pick,)
    print(f"  [warn] {pick!r} not found; falling back to first single-salt sheet.")
    return (NF270_SINGLE_SALT_RUNS[0],)


# =============================================================================
# TRUNK — pick the workflow
# =============================================================================

# Trunk recipes map a single name to (multistart, fim_uncertainty,
# doe_after_fit, mass_litmus_test).
TRUNK_RECIPES = {
    "simulate":        dict(multistart=False, fim=False, doe=False, litmus=False, sim_only=True),
    "fit":             dict(multistart=False, fim=False, doe=False, litmus=False, sim_only=False),
    "fit_multistart":  dict(multistart=True,  fim=False, doe=False, litmus=False, sim_only=False),
    "fit_FIM":         dict(multistart=True,  fim=True,  doe=False, litmus=False, sim_only=False),
    "fit_FIM_DoE":     dict(multistart=True,  fim=True,  doe=True,  litmus=False, sim_only=False),
    "mass_litmus":     dict(multistart=False, fim=False, doe=False, litmus=True,  sim_only=False),
}


def pick_trunk(root: str = "DATA3") -> dict:
    """Pick the trunk workflow.

    Default differs by root:
      DATA1 / DATA2 → 'manifest_default' (reproduce published-paper figures
                       exactly with no override)
      DATA3 / NF270 / Custom → 'fit_multistart' (the recommended path)
    """
    print("\n" + "=" * 60)
    print("TRUNK  —  pick the workflow")
    print("=" * 60)
    if root in ("DATA1", "DATA2"):
        print("  0. Manifest defaults  (reproduce paper figures exactly — RECOMMENDED)")
    print("  1. Simulate only       (forward-only, no fit)")
    print("  2. Fit                 (single-shot ParmEst)")
    print("  3. Fit + multistart    (recommended for DATA3)")
    print("  4. Fit + multistart + FIM uncertainty")
    print("  5. Fit + multistart + FIM + Pyomo.DoE / model discrimination")
    print("  6. Mass-litmus test    (sigma=0, B dropped, closed-form np regression — DATA3 only)")
    default = "0" if root in ("DATA1", "DATA2") else "3"
    choice = _prompt("Trunk", default).strip()
    name = {"0": "manifest_default", "1": "simulate", "2": "fit",
            "3": "fit_multistart", "4": "fit_FIM", "5": "fit_FIM_DoE",
            "6": "mass_litmus"}.get(choice, "manifest_default" if root in ("DATA1", "DATA2") else "fit_multistart")
    if name == "manifest_default":
        return {"name": "manifest_default",
                "multistart": None, "fim": None, "doe": False,
                "litmus": False, "sim_only": False}
    recipe = dict(TRUNK_RECIPES[name])
    recipe["name"] = name
    return recipe


# =============================================================================
# BRANCHES — pick which outputs to generate
# =============================================================================

def pick_branches() -> set[str]:
    print("\n" + "=" * 60)
    print("BRANCHES  —  pick outputs (comma-separated, or 'all')")
    print("=" * 60)
    print("  m   Mass-vs-time plots (one per vial / file)")
    print("  c   Concentration-vs-time plots")
    print("  p   Parameter table")
    print("  f   FIM heatmap / sigma-sensitivity contour")
    print("  d   DoE next-experiment recommendations")
    print("  all all of the above")
    raw = _prompt("Branches", "all").strip().lower()
    if raw in ("all", "*", ""):
        return {"m", "c", "p", "f", "d"}
    parts = [p.strip() for p in raw.replace(",", " ").split() if p.strip()]
    valid = {"m", "c", "p", "f", "d"}
    chosen = {p for p in parts if p in valid}
    return chosen or {"m", "c"}


# =============================================================================
# Dispatch — the runner
# =============================================================================

def _dispatch_paper_campaign(root: str, subset, trunk: dict, branches: set[str]) -> None:
    """DATA1 / DATA2 — uses materialize_all dispatcher (manifest-driven).

    IMPORTANT — reproducibility guarantee:
    -------------------------------------
    DATA1 and DATA2 produce published-paper figures whose numerical
    results must match prior runs exactly. The manifest already encodes
    the right multistart / FIM / B-form / mode settings per FigureSpec.
    This dispatch deliberately does NOT pass any extra_opts overrides
    when the user picks the default trunk recipe — materialize_all is
    called with exactly the same arguments the older runfile used.

    The mass-litmus toggle is also explicitly NOT activated for DATA1 /
    DATA2 even if the user picks trunk=mass_litmus, because the toggle
    is NF270-specific (NF270_MASS_LITMUS_TEST_ACTIVE) and would only
    affect a sheet whose loader returns NF270-style data_stru. We still
    print a warning in that case.
    """
    if root == "DATA1":
        save_dir = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                    / "paper_artifacts" / "data1" / "notebook_figures")
        data_root = DATA1_ROOT
    else:
        save_dir = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                    / "paper_artifacts" / "data2" / "notebook_figures")
        data_root = DATA2_ROOT
        _ensure_data2_library()
    save_dir.mkdir(parents=True, exist_ok=True)

    if trunk.get("litmus"):
        print(f"\n[note] mass-litmus is NF270-specific. Ignoring for {root}; "
              f"running the manifest default workflow instead.")

    # Reproducibility guarantee: pass nothing through extra_opts unless the
    # user explicitly deviated from the manifest defaults. The manifest's
    # baked-in RunRequest settings are what produced the published figures.
    materialize_kwargs = {
        "campaign":  root,
        "save_dir":  save_dir,
        "data_root": data_root,
        "only":      subset,
    }
    extra_opts = {}
    if trunk["name"] == "simulate":
        extra_opts["sim_only_override"] = True
    if trunk["name"] in ("fit",):
        # Force off multistart for the user-requested single-shot fit.
        extra_opts["multistart_override"] = False
    if trunk["name"] in ("fit_FIM_DoE",):
        # Append DoE on top of the manifest defaults.
        extra_opts["compute_doe"] = True
    if extra_opts:
        materialize_kwargs["extra_opts"] = extra_opts

    results = ucb.materialize_all(**materialize_kwargs)
    _print_outputs(root, results)
    _print_coverage(root, save_dir)


def _dispatch_nf270(subset, trunk: dict, branches: set[str]) -> None:
    """DATA3 / NF270 — also through materialize_all."""
    if not hasattr(ucb, "NF270_RUN_REGISTRY"):
        print("\n[error] NF270_RUN_REGISTRY is not defined in the library.")
        return
    if not NF270_ROOT.exists():
        print(f"\n[error] NF270 data root does not exist: {NF270_ROOT}")
        print("        Set DIAFILTRATION_NF270_ROOT to the directory containing")
        print("        NF270_MC2.xlsx, NF270_MC3.xlsx, NF270_MC4.xlsx, NF270_MC5.xlsx.")
        return

    save_dir = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                / "paper_artifacts" / "nf270" / "campaign_figures")
    save_dir.mkdir(parents=True, exist_ok=True)
    print(f"  data root:  {NF270_ROOT}")
    print(f"  save dir:   {save_dir}")
    print(f"  registered: {len(ucb.NF270_RUN_REGISTRY)} green-flagged sheets")

    extra_opts = {
        "multistart_override": trunk.get("multistart", False),
        "fim_override":        trunk.get("fim", False),
        "branches":            sorted(branches),
    }
    if trunk.get("doe"):
        extra_opts["compute_doe"] = True
    if trunk.get("litmus"):
        extra_opts["mass_litmus_test"] = True

    # Toggle the module-level mass-litmus flag for the duration of this run.
    prior_litmus = getattr(ucb, "NF270_MASS_LITMUS_TEST_ACTIVE", False)
    if trunk.get("litmus"):
        ucb.NF270_MASS_LITMUS_TEST_ACTIVE = True
        print(f"  trunk:      mass-litmus test (sigma=0, closed-form Lp regression)")
    try:
        results = ucb.materialize_all(
            campaign="NF270", save_dir=save_dir, data_root=NF270_ROOT,
            only=subset, extra_opts=extra_opts,
        )
    finally:
        if hasattr(ucb, "NF270_MASS_LITMUS_TEST_ACTIVE"):
            ucb.NF270_MASS_LITMUS_TEST_ACTIVE = prior_litmus

    _print_outputs("NF270", results)
    _print_coverage("NF270", save_dir)


def _dispatch_custom(trunk: dict, branches: set[str]) -> None:
    """Single-file run — accepts .mat or .xlsx."""
    raw = input("\nPaste one experiment file path (.mat or .xlsx): ").strip()
    if not raw:
        print("No file selected.")
        return
    file_path = Path(raw).expanduser().resolve()
    if not file_path.exists():
        print(f"[error] not found: {file_path}")
        return

    is_xlsx = file_path.suffix.lower() in {".xlsx", ".xls"}
    selector = None
    if is_xlsx:
        sel = _prompt("Excel sheet selector (sheet name or index)", "0")
        selector = int(sel) if sel.isdigit() else sel

    mode = _prompt("Model mode (DATA / Lag / Overflow)",
                    "Lag" if is_xlsx else "DATA")
    workflow_family = _prompt("Workflow family (DATA1/DATA2/DATA3)",
                               "DATA3" if is_xlsx else "DATA2").strip().upper()
    B_form = _prompt("B form (single / pervial / convection / 0..3)", "single")

    use_parmest = (trunk.get("name") != "simulate")
    uq_method = "fim" if trunk.get("fim") else None
    if trunk.get("litmus"):
        # Closed-form path doesn't go through Pyomo.
        data_stru = ucb.loadxlsx(file_path, sheet=selector)["data_stru"] if is_xlsx \
                    else ucb.loadmat(file_path)["data_stru"]
        fit_stru, sim_stru = ucb.solve_mass_balance_only(data_stru, LOUD=True)
        print(f"\nFit:  Lp = {fit_stru['parameters']['Lp']:.4f} L/m²/hr/bar")
        print(f"      RMSE = {fit_stru['rmse_g']:.4f} g over {fit_stru['n_points']} points")
        if "m" in branches or "c" in branches:
            ucb.plot_sim_comparison(data_stru, sim_stru,
                                     plot_pred=True, lg=True, LOUD=True)
        return

    results = ucb.run_pipeline(
        file_path,
        mode=mode, workflow_family=workflow_family, B_form=B_form,
        selector=selector, use_parmest=use_parmest,
        uncertainty_method=uq_method,
    )
    if is_xlsx and workflow_family == "DATA3":
        figs = ucb.run_data3_time_series_plots(
            results.to_dict(), save_dir=DATA3_CUSTOM_FIGURES, show=False)
        if figs:
            print("\nDATA3 plots:")
            for f in figs:
                print(f"  + {f}")

    if "p" in branches:
        json_dir = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "custom_runs"
        json_dir.mkdir(parents=True, exist_ok=True)
        ucb.report_parameters(results,
                                save_path=json_dir / f"{file_path.stem}_parameters.json")
        if results.uncertainty:
            ucb.report_uncertainty(results,
                                     save_path=json_dir / f"{file_path.stem}_uncertainty.json")

    print("\n--- StageResults summary ---")
    for field_name in ("data_file", "parameters", "is_batch"):
        v = getattr(results, field_name, None)
        if v not in (None, {}, []):
            print(f"  {field_name}: {v!r}")


# =============================================================================
# Main — walk the tree
# =============================================================================

def main() -> None:
    print(__doc__.split("Environment")[0])

    root = pick_root()
    subset = pick_subset(root)
    trunk = pick_trunk(root)
    branches = pick_branches()

    print("\n" + "=" * 60)
    print(f"DISPATCHING")
    print("=" * 60)
    print(f"  root      = {root}")
    print(f"  subset    = {subset if subset else 'ALL'}")
    print(f"  trunk     = {trunk['name']}  "
          f"(multistart={trunk['multistart']}, FIM={trunk['fim']}, "
          f"DoE={trunk['doe']}, litmus={trunk['litmus']})")
    print(f"  branches  = {sorted(branches)}")

    if root in ("DATA1", "DATA2"):
        _dispatch_paper_campaign(root, subset, trunk, branches)
    elif root == "DATA3":
        _dispatch_nf270(subset, trunk, branches)
    else:
        _dispatch_custom(trunk, branches)


if __name__ == "__main__":
    main()
