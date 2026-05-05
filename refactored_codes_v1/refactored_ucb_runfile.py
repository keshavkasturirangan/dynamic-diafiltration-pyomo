#!/usr/bin/env python3
"""
Standalone CLI for the UCB diafiltration refactor (full-paper edition).

Replaces the previous runfile. Differences:
- Option 1 (DATA1): runs the full DATA1 manifest (no `only=` filter) so
  every figure listed in the DATA1 main + SI papers is produced.
- Option 2 (DATA2): runs the full DATA2 manifest (drops the previous
  4-entry whitelist). Includes cross-verification per-case fits if the
  cross_verification_patch_block has been applied.
- Both options call `report_paper_coverage` from the
  paper_coverage_patch_block at the end and print a checklist showing
  every paper figure + table the campaign is supposed to produce, with
  found/missing status.
- A "fast subset" prompt is offered for users who only want the
  cheap-to-render figures (the previous default).

The library is imported as `refactored_ucb_library` - which can be
either the standalone v2 file (renamed) or the original library with
the cross_verification + paper_coverage patches applied.
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


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA1_ROOT = Path(
    os.environ.get(
        "DIAFILTRATION_DATA1_ROOT",
        REPO_ROOT / "legacy" / "data1_matlab" / "data",
    )
).expanduser().resolve()
DATA2_ROOT = Path(
    os.environ.get(
        "DIAFILTRATION_DATA2_ROOT",
        REPO_ROOT / "legacy" / "data1_matlab" / "data_library",
    )
).expanduser().resolve()
DATA3_CUSTOM_FIGURES = REPO_ROOT / "UnifiedFramework" / "DATA3" / "figures" / "data3_option3"

# The previous runfile's whitelist - kept around as the "fast subset"
# for users who don't want the full paper reproduction every time.
DATA2_FAST_SUBSET = (
    "calibration_plots",
    "pressure_changes",
    "model_error_visualization",
    "publication_figures",
)
DATA1_FAST_SUBSET = (
    "figure_2",
    "figure_3",
    "data_analysis",
)


# ---- helpers ----------------------------------------------------------------

def _ensure_data2_library() -> None:
    sentinel = DATA2_ROOT / "data_stru-dataset270511.123.mat"
    if sentinel.exists():
        return
    archive_path = REPO_ROOT / "legacy" / "data1_matlab" / "inputdata_mat_files.zip"
    if not archive_path.exists():
        raise FileNotFoundError(
            f"DATA2 source folder is missing at {DATA2_ROOT} and archive was not "
            f"found at {archive_path}."
        )
    DATA2_ROOT.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path, "r") as zip_ref:
        zip_ref.extractall(DATA2_ROOT.parent)


def _prompt(message: str, default: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    value = input(f"{message}{suffix}: ").strip()
    return value or default


def _prompt_yes_no(message: str, default: str = "n") -> bool:
    return _prompt(message + " (y/n)", default).lower().startswith("y")


def _print_outputs(label: str, results: list[dict]) -> None:
    ok = [r for r in results if r.get("status") == "ok"]
    err = [r for r in results if r.get("status") == "error"]
    print(f"\n[{label}] manifest entries attempted: {len(results)}  "
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
    """Print the paper-figure coverage checklist if the patch is applied."""
    fn = getattr(ucb, "report_paper_coverage", None)
    fmt = getattr(ucb, "format_paper_coverage", None)
    if fn is None or fmt is None:
        print(f"\n[note] paper_coverage_patch_block is not loaded into the library, "
              f"so no per-figure checklist is available. Apply the patch with "
              f"`python apply_paper_coverage_patch.py` to enable this output.")
        return
    report = fn(campaign, save_dir)
    print(fmt(report))


# ---- option 1: DATA1 --------------------------------------------------------

def _run_data1() -> None:
    print("\nRunning DATA1 paper reproduction (full main + SI)...")
    save_dir = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                / "paper_artifacts" / "data1" / "notebook_figures")
    save_dir.mkdir(parents=True, exist_ok=True)

    fast = _prompt_yes_no("Use fast subset (skip slow contour panels)?", "n")
    only = DATA1_FAST_SUBSET if fast else None
    results = ucb.materialize_all(
        campaign="DATA1",
        save_dir=save_dir,
        data_root=DATA1_ROOT,
        only=only,
    )
    _print_outputs("DATA1", results)
    _print_coverage("DATA1", save_dir)


# ---- option 2: DATA2 --------------------------------------------------------

def _run_data2() -> None:
    print("\nRunning DATA2 paper reproduction (full main + SI)...")
    _ensure_data2_library()
    save_dir = (REPO_ROOT / "UnifiedFramework" / "DATA3" / "results"
                / "paper_artifacts" / "data2" / "notebook_figures")
    save_dir.mkdir(parents=True, exist_ok=True)

    fast = _prompt_yes_no(
        "Use fast subset (skip cross-verification + model_variations)?", "n"
    )
    only = DATA2_FAST_SUBSET if fast else None
    results = ucb.materialize_all(
        campaign="DATA2",
        save_dir=save_dir,
        data_root=DATA2_ROOT,
        only=only,
    )
    _print_outputs("DATA2", results)
    _print_coverage("DATA2", save_dir)


# ---- option 3: custom one-file run ------------------------------------------

def _choose_custom_file() -> Path | None:
    raw = input("\nPaste one experiment file path (.mat or .xlsx): ").strip()
    if not raw:
        return None
    return Path(raw).expanduser().resolve()


def _run_custom() -> None:
    file_path = _choose_custom_file()
    if file_path is None:
        print("No file was selected.")
        return

    selector: object = None
    if file_path.suffix.lower() in {".xlsx", ".xls"}:
        print("Detected Excel input; loading via the Excel data path.")
        selector_raw = _prompt("Excel sheet selector (sheet name or index)", "0")
        if selector_raw.strip():
            selector = int(selector_raw) if selector_raw.strip().isdigit() else selector_raw.strip()
        else:
            selector = 0

    print("\nModel choices:")
    mode = _prompt("Model mode", "DATA")
    default_family = "DATA3" if file_path.suffix.lower() in {".xlsx", ".xls"} else "DATA2"
    workflow_family = _prompt("Model family (DATA1/DATA2/DATA3)", default_family).strip().upper()
    B_form = _prompt("B form (single/exp/linear)", "single").strip().lower()

    print("\nFitting choices:")
    use_parmest = _prompt_yes_no("Use Pyomo ParmEst?", "n")
    cached_str = _prompt("Path to cached fit_stru.mat (blank=fit fresh)", "")
    cached_fit_path = Path(cached_str).expanduser().resolve() if cached_str.strip() else None

    print("\nUncertainty choices:")
    uq_method = _prompt("Uncertainty method (fim/cov_est/none)", "fim").strip().lower()
    if uq_method == "none":
        uq_method = None

    results = ucb.run_pipeline(
        file_path,
        mode=mode,
        workflow_family=workflow_family,
        B_form=B_form,
        selector=selector,
        use_parmest=use_parmest,
        uncertainty_method=uq_method,
        cached_fit_path=cached_fit_path,
    )

    if file_path.suffix.lower() in {".xlsx", ".xls"} and workflow_family == "DATA3":
        # DATA3 time-series via the legacy renderer.
        figure_outputs = ucb.run_data3_time_series_plots(
            results.to_dict(), save_dir=DATA3_CUSTOM_FIGURES, show=False
        )
        if figure_outputs:
            print("\nCreated DATA3 plots:")
            for item in figure_outputs:
                print(f"  + {item}")

    save_results = _prompt_yes_no(
        "\nSave numeric results (parameters + uncertainty) to JSON?", "y"
    )
    if save_results:
        json_dir = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "custom_runs"
        json_dir.mkdir(parents=True, exist_ok=True)
        params_json = json_dir / f"{file_path.stem}_parameters.json"
        uq_json = json_dir / f"{file_path.stem}_uncertainty.json"
        ucb.report_parameters(results, save_path=params_json)
        if results.uncertainty:
            ucb.report_uncertainty(results, save_path=uq_json)
        print(f"  + parameters: {params_json}")
        if results.uncertainty:
            print(f"  + uncertainty: {uq_json}")

    print("\n--- StageResults summary ---")
    for field_name in ("data_file", "parameters", "is_batch"):
        value = getattr(results, field_name, None)
        if value not in (None, {}, []):
            print(f"  {field_name}: {value!r}")
    if results.uncertainty:
        scalar_keys = {k: results.uncertainty[k]
                        for k in ("trace", "det", "method")
                        if k in results.uncertainty}
        print(f"  uncertainty (scalars): {scalar_keys!r}")


# ---- mode picker ------------------------------------------------------------

def _choose_mode() -> str:
    print("\nChoose a workflow:")
    print("  1. DATA1 paper reproduction (main + SI figures + Table 1)")
    print("  2. DATA2 paper reproduction (main + SI figures + tables)")
    print("  3. Custom one-file run")
    choice = _prompt("Workflow", "1").strip().lower()
    if choice in {"1", "data1"}:
        return "DATA1"
    if choice in {"2", "data2"}:
        return "DATA2"
    return "CUSTOM"


def main() -> None:
    workflow = _choose_mode()
    if workflow == "DATA1":
        _run_data1()
    elif workflow == "DATA2":
        _run_data2()
    else:
        _run_custom()


if __name__ == "__main__":
    main()
