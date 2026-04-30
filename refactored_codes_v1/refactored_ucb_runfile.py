#!/usr/bin/env python3
"""Simple runner for the refactored diafiltration workflow.

The runner stays thin:
- choose DATA1 or DATA2 to recreate the paper-style plots from the notebooks
- choose custom to run one experimental file through the stage-based workflow

The refactored folder is self-contained:
- `refactored_ucb_library.py` holds the workflow functions
- `conductivity_paper.py` sits beside it and handles conductivity-to-
  concentration conversion
"""

from __future__ import annotations

import os
from pathlib import Path
import sys
import zipfile
import matplotlib.pyplot as plt

# Add the script folder to Python's import path so Spyder can run it directly.
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from refactored_ucb_library import (
    run_data2_notebook_workflow,
    run_data1_notebook_workflow,
    run_workflow,
)


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

DATA1_MAIN_FIGURES = [
    "figure_2.png",
    "figure_4.png",
    "figure_5.png",
    "mass-dat501.1.png",
    "concentration-dat501.1.png",
    "mass-dat511.12.png",
    "concentration-dat511.12.png",
    "concentration_range.png",
    "sigma_sensitivity-mass.png",
    "sigma_sensitivity-reten_conc.png",
    "sigma_sensitivity-perme_conc.png",
    "contour_fixsig-mass.png",
    "contour_fixsig-retentate_conc.png",
    "contour_fixsig-permeate_conc.png",
    "contour_fixB-mass.png",
    "contour_fixB-retentate_conc.png",
    "contour_fixB-permeate_conc.png",
]

DATA1_SI_FIGURES = [
    "figure_s2.png",
    "data1_si_reduced_diafiltration.png",
    "data1_si_diafiltration_B.png",
    "data1_si_reduced_filtration.png",
    "data1_si_filtration_sigma.png",
    "data1_si_filtration_B.png",
    "figure_s3.png",
]

DATA2_SI_FIGURES = [
    "calib_curve.png",
    "mass-dat270611.121.png",
    "concentration-dat270611.121.png",
    "mass-dat270711.121.png",
    "concentration-dat270711.121.png",
    "mass-dat270511.221.png",
    "concentration-dat270511.221.png",
    "mass-dat270511.321.png",
    "concentration-dat270511.321.png",
    "mass-dat270511.421.png",
    "concentration-dat270511.421.png",
    "mass-dat270511.921.png",
    "concentration-dat270511.921.png",
    "mass-dat270511.521.png",
    "concentration-dat270511.521.png",
    "mass-dat270511.621.png",
    "concentration-dat270511.621.png",
    "mass-dat270511.721.png",
    "concentration-dat270511.721.png",
    "mass-dat270511.821.png",
    "concentration-dat270511.821.png",
    "Bpervial.png",
    "Js_Jw_cin.png",
    "Js_predict.png",
    "Js_predict0.png",
    "Js_predict1.png",
    "Jw_predict.png",
]


def _ensure_data2_library() -> None:
    """Extract the bundled DATA2 source folder when it is missing on a machine."""
    sentinel = DATA2_ROOT / "data_stru-dataset270511.123.mat"
    if sentinel.exists():
        return

    archive_path = REPO_ROOT / "legacy" / "data1_matlab" / "inputdata_mat_files.zip"
    if not archive_path.exists():
        raise FileNotFoundError(
            f"DATA2 source folder is missing at {DATA2_ROOT} and archive was not found at {archive_path}."
        )

    DATA2_ROOT.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path, "r") as zip_ref:
        zip_ref.extractall(DATA2_ROOT.parent)


def _find_generated_figure(filename: str, outputs: list[str], save_dir: Path) -> Path | None:
    """Resolve one generated figure by basename across the known output folders."""
    output_paths = [Path(item) for item in outputs]
    for path in output_paths:
        if path.name == filename and path.exists():
            return path

    candidates = [
        save_dir / filename,
        REPO_ROOT / "UnifiedFramework" / "DATA3" / "figures" / filename,
        REPO_ROOT / "figures" / filename,
        REPO_ROOT / filename,
        Path.cwd() / filename,
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def _find_generated_figures(pattern: str, outputs: list[str], save_dir: Path) -> list[Path]:
    """Resolve one generated figure or wildcard pattern across the known output folders."""
    if "*" not in pattern:
        fig_path = _find_generated_figure(pattern, outputs, save_dir)
        return [fig_path] if fig_path is not None else []

    output_paths = [Path(item) for item in outputs]
    matches = [path for path in output_paths if path.match(pattern) and path.exists()]
    if matches:
        return sorted(dict.fromkeys(matches))

    search_roots = [
        save_dir,
        REPO_ROOT / "UnifiedFramework" / "DATA3" / "figures",
        REPO_ROOT / "figures",
        REPO_ROOT,
        Path.cwd(),
    ]
    found: list[Path] = []
    for root in search_roots:
        found.extend(sorted(root.glob(pattern)))
    return sorted(dict.fromkeys(path for path in found if path.exists()))


def _show_paper_figures(outputs: list[str], save_dir: Path, figure_names: list[str], title: str) -> None:
    """Show one published-paper figure set in matplotlib windows."""
    shown = 0
    for filename in figure_names:
        for fig_path in _find_generated_figures(filename, outputs, save_dir):
            try:
                image = plt.imread(fig_path)
            except Exception:
                continue
            shown += 1
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.imshow(image)
            ax.axis("off")
            ax.set_title(f"{title} {shown}: {fig_path.name}", fontsize=10)
    if shown:
        plt.show()


def _prompt(prompt: str, default: str) -> str:
    value = input(f"{prompt} [{default}]: ").strip()
    return value or default


def _choose_mode() -> str:
    # Show the user the three simple ways to run the code.
    print("\nChoose a simple workflow:")
    print("  1. DATA1 paper plots")
    print("  2. DATA2 paper plots")
    print("  3. Custom one-file run")
    # Ask for a choice and default to DATA1 if they just press Enter.
    choice = _prompt("Workflow", "1").strip().lower()
    if choice in {"1", "data1"}:
        return "DATA1"
    if choice in {"2", "data2"}:
        return "DATA2"
    return "CUSTOM"


def _run_data1() -> None:
    """Recreate the paper-style plots for DATA1."""
    print("\nRunning DATA1 paper reproduction...")
    outputs = run_data1_notebook_workflow(
        data_root=DATA1_ROOT,
        save_dir=REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "data1" / "notebook_figures",
        show=True,
    )
    print("\nCreated outputs:")
    for item in outputs:
        print(f"  - {item}")


def _run_data2() -> None:
    """Recreate the paper-style plots for DATA2."""
    print("\nRunning DATA2 paper reproduction...")
    _ensure_data2_library()
    save_dir = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "data2" / "notebook_figures"
    outputs = run_data2_notebook_workflow(
        data_root=DATA2_ROOT,
        save_dir=save_dir,
        show=False,
        fast_mode=True,
    )
    _show_paper_figures(outputs, save_dir, DATA2_SI_FIGURES, "DATA2 SI")
    print("\nCreated outputs:")
    for item in outputs:
        print(f"  - {item}")


def _choose_custom_file() -> Path | None:
    # Ask for one MATLAB file path so the flow stays close to the staged example.
    raw = input("\nPaste one .mat file path: ").strip()
    if not raw:
        return None
    return Path(raw).expanduser().resolve()


def _run_custom() -> None:
    """Run one experimental file through the staged workflow."""
    file_path = _choose_custom_file()
    if file_path is None:
        print("No file was selected.")
        return

    # Ask for the key choices, keeping the prompt list short and explicit.
    mode = _prompt("Model mode", "DATA")
    workflow_family = _prompt("Model family", "DATA1 or DATA2").strip().upper()
    use_parmest = _prompt("Use ParmEst? (y/n)", "n").lower().startswith("y")
    uncertainty_method = _prompt("Uncertainty method (fim/cov_est)", "fim").strip().lower()

    results = run_workflow(
        file_path,
        mode=mode,
        workflow_family=workflow_family,
        use_parmest=use_parmest,
        uncertainty_method=uncertainty_method,
    )
    print("\nCreated results:")
    for key in results:
        if not key.startswith("_"):
            print(f"  - {key}")


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
