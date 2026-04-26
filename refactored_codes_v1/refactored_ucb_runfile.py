#!/usr/bin/env python3
"""Simple runner for the refactored diafiltration workflow.

The runner is intentionally plain:
- choose DATA1 or DATA2 to recreate the paper-style plots
- choose custom to run one experimental file through the model
"""

from __future__ import annotations

from pathlib import Path
import sys

# Add the script folder to Python's import path so Spyder can run it directly.
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from refactored_ucb_library import (
    calc_FIM,
    loadmat,
    _normalize_conductivity_measurements,
    run_campaign,
    solve_model,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


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
    # Let the library handle the full DATA1 paper workflow.
    outputs = run_campaign("DATA1", save_dir=REPO_ROOT / "UnifiedFramework" / "DATA3" / "figures" / "data1")
    print("\nCreated outputs:")
    for item in outputs:
        print(f"  - {item}")


def _run_data2() -> None:
    """Recreate the paper-style plots for DATA2."""
    print("\nRunning DATA2 paper reproduction...")
    # Let the library handle the full DATA2 paper workflow.
    outputs = run_campaign("DATA2", save_dir=REPO_ROOT / "UnifiedFramework" / "DATA3" / "figures" / "data2")
    print("\nCreated outputs:")
    for item in outputs:
        print(f"  - {item}")


def _choose_custom_files() -> list[Path]:
    # Let the user paste one or more file paths, separated by commas.
    raw = input("\nPaste one or more .mat file paths separated by commas: ").strip()
    if not raw:
        return []
    # Clean up the user input and turn each item into a real Path object.
    return [Path(item.strip()).expanduser().resolve() for item in raw.split(",") if item.strip()]


def _run_custom() -> None:
    """Run one or more experimental files through the existing model code."""
    files = _choose_custom_files()
    if not files:
        print("No files were selected.")
        return

    # Ask the user for the basic modeling choices.
    mode = _prompt("Model mode", "DATA")
    workflow_family = _prompt("Model family", "DATA1 or DATA2").strip().upper()
    run_fim = _prompt("Also compute FIM? (y/n)", "n").lower().startswith("y")

    for file_path in files:
        # Load the MATLAB file into the nested Python dictionary used by the model.
        data_stru = loadmat(str(file_path)).get("data_stru")
        if data_stru is None:
            print(f"Skipping {file_path}: no data_stru entry was found.")
            continue

        # Convert conductivity to concentration when the file asks for it.
        _normalize_conductivity_measurements(data_stru)

        # Fit the first-principles model to this data file.
        print(f"\nFitting {file_path.name}...")
        fit_stru, sim_stru, sim_inter = solve_model(data_stru, mode, workflow_family=workflow_family)
        print(f"Finished fit for dataset {data_stru.get('dataset')}")
        if fit_stru and "parameters" in fit_stru:
            print(f"Parameters: {fit_stru['parameters']}")

        if run_fim:
            # Compute the Fisher Information Matrix if the user asked for it.
            print("Computing FIM...")
            fim = calc_FIM(data_stru, mode, theta=fit_stru["parameters"] if fit_stru else None, workflow_family=workflow_family)
            print(f"FIM trace: {fim.get('trace')}")
            print(f"FIM det: {fim.get('det')}")


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
