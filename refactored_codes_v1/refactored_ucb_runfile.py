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

from pathlib import Path
import sys

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
DATA1_ROOT = REPO_ROOT / "legacy" / "data1_matlab" / "data"


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
    outputs = run_data2_notebook_workflow(
        data_root=REPO_ROOT / "legacy" / "data1_matlab" / "data_library",
        save_dir=REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "data2" / "notebook_figures",
        show=True,
        fast_mode=True,
    )
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
