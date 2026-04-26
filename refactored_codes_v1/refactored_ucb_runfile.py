#!/usr/bin/env python3
"""Easy-to-follow runner for the refactored diafiltration workflow.

The idea is simple:
1. Recreate DATA1 or DATA2 paper figures.
2. Or run one custom file.

The science stays in `refactored_ucb_library.py` and the unified DATA3
reproduction script.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Sequence

import pandas as pd

try:
    from .refactored_ucb_library import DataLoader, ModelOptions, UQEngine, WorkflowFamily, run_paper_reproduction
except ImportError:  # pragma: no cover - allows running as a plain script
    from refactored_ucb_library import DataLoader, ModelOptions, UQEngine, WorkflowFamily, run_paper_reproduction


REPO_ROOT = Path(__file__).resolve().parents[1]


def _prompt_text(prompt: str, default: str) -> str:
    value = input(f"{prompt} [{default}]: ").strip()
    return value or default


def _prompt_yes_no(prompt: str, default: bool = False) -> bool:
    hint = "Y/n" if default else "y/N"
    value = input(f"{prompt} ({hint}): ").strip().lower()
    if not value:
        return default
    return value in {"y", "yes", "1", "true"}


def _prompt_choice(prompt: str, choices: Sequence[str], default: str) -> str:
    choice_text = ", ".join(f"{i + 1}={choice}" for i, choice in enumerate(choices))
    value = input(f"{prompt} ({choice_text}) [{default}]: ").strip()
    if not value:
        return default
    if value.isdigit():
        idx = int(value) - 1
        if 0 <= idx < len(choices):
            return choices[idx]
    if value in choices:
        return value
    print(f"Unrecognized choice {value!r}; using {default!r}.")
    return default


def _select_files_from_list(paths: Sequence[Path]) -> List[Path]:
    print("\nAvailable files:")
    for i, path in enumerate(paths, start=1):
        print(f"  {i}. {path}")
    raw = input(
        "\nType numbers like 1,3 or paste file paths separated by commas.\n"
        "Press Enter to use all listed files: "
    ).strip()
    if not raw:
        return list(paths)
    if all(token.strip().isdigit() for token in raw.split(",")):
        chosen: List[Path] = []
        for token in raw.split(","):
            idx = int(token.strip()) - 1
            if 0 <= idx < len(paths):
                chosen.append(paths[idx])
        return chosen
    chosen_paths: List[Path] = []
    for token in raw.split(","):
        token = token.strip()
        if token:
            chosen_paths.append(Path(token).expanduser().resolve())
    return chosen_paths


def _parse_file_args(raw_files: Sequence[str]) -> List[Path]:
    return [Path(item).expanduser().resolve() for item in raw_files]


def _load_tabular_file(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    return pd.read_csv(path)


def _default_data1_files(loader: DataLoader) -> List[Path]:
    registry = loader.data1_registry()
    return [loader.resolve(path) for path in registry["core_mat"]]


def _make_options(family: WorkflowFamily, mode: str, output_dir: Path) -> ModelOptions:
    return ModelOptions(
        model_family=family,
        mode=mode,
        output_dir=output_dir,
        plot_pred=True,
        lg=True,
        LOUD=False,
        sigma_fixed=(family == WorkflowFamily.DATA1),
    )


def _paper_workflow_mode() -> str:
    print("\nQuick start:")
    print("  1. Recreate DATA1 paper plots")
    print("  2. Recreate DATA2 paper plots")
    print("  3. Custom one-file run")
    choice = _prompt_text("Quick start", "1").strip().lower()
    if choice in {"1", "data1"}:
        return "DATA1"
    if choice in {"2", "data2"}:
        return "DATA2"
    return "CUSTOM"


def _choose_workflow_profile() -> WorkflowFamily:
    print("\nChoose the workflow family:")
    print("  1. DATA1")
    print("  2. DATA2")
    choice = _prompt_text("Workflow family", "1").strip()
    return WorkflowFamily.DATA1 if choice in {"1", "DATA1", "data1"} else WorkflowFamily.DATA2


def _run_paper_shortcut(paper_family: str, output_dir: Optional[Path] = None) -> None:
    """Hand off to the existing unified reproduction script."""
    print("\nRunning the paper reproduction workflow:")
    print(f"  paper: {paper_family}")
    result = run_paper_reproduction(paper_family, repo_root=REPO_ROOT, output_root=output_dir)
    print(f"  script: {result['script']}")
    print(f"  output root: {result['output_root']}")


def _choose_input_files(loader: DataLoader, family: WorkflowFamily) -> List[Path]:
    if family == WorkflowFamily.DATA1:
        print("\nDATA1 preset files are available, but you can still type your own paths.")
        return _select_files_from_list(_default_data1_files(loader))
    raw = input("\nPaste one or more DATA2 file paths separated by commas: ").strip()
    if not raw:
        raise ValueError("No input files were provided.")
    return _parse_file_args([token.strip() for token in raw.split(",") if token.strip()])


def _choose_run_depth() -> str:
    print("\nWhat do you want to do?")
    print("  1. Estimate parameters only")
    print("  2. Estimate parameters + uncertainty quantification")
    print("  3. Estimate parameters + suggest new experiments")
    print("  4. Full workflow (fit + UQ + suggestions + plots)")
    choice = _prompt_text("Run type", "1").strip()
    if choice in {"1", "fit", "fit_only"}:
        return "fit_only"
    if choice in {"2", "uq", "uncertainty", "uncertainty_quantification"}:
        return "uq"
    if choice in {"3", "suggest", "suggestions", "design"}:
        return "suggest"
    return "full"


def _maybe_convert_conductivity(loader: DataLoader, input_path: Path, output_dir: Path) -> Optional[Path]:
    is_table = input_path.suffix.lower() in {".csv", ".xlsx", ".xls"}
    if not is_table:
        return None

    if not _prompt_yes_no("Does this file contain conductivity that should be converted to concentration?", default=False):
        return None

    column = _prompt_text("Conductivity column name", "Conductivity")
    salt_name = _prompt_text("Salt name for paper model autofill", "NaCl")
    model = _prompt_choice("Conductivity model", ("variant_shedlovsky", "msa"), "variant_shedlovsky")
    output_units = _prompt_choice("Output units", ("mM", "M"), "mM")
    return _convert_conductivity_file(
        loader,
        input_path,
        output_dir=output_dir,
        column=column,
        model=model,
        salt_name=salt_name,
        output_units=output_units,
    )


def _convert_conductivity_file(
    loader: DataLoader,
    file_path: Path,
    *,
    output_dir: Path,
    column: str,
    model: str,
    salt_name: str,
    output_units: str,
) -> Path:
    # Load the spreadsheet or CSV as a plain table.
    frame = _load_tabular_file(file_path)
    # Ask the shared library to convert conductivity into concentration.
    converted = loader.convert_conductivity_frame(
        frame,
        conductivity_column=column,
        model=model,
        salt_name=salt_name,
        output_units=output_units,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{file_path.stem}_converted.csv"
    converted.to_csv(out_path, index=False)
    print("\nConductivity conversion complete:")
    print(f"  input: {file_path}")
    print(f"  output: {out_path}")
    return out_path


def _run_one(
    engine: UQEngine,
    file_path: Path,
    options: ModelOptions,
    *,
    run_fim: bool = False,
    run_suggest: bool = False,
    run_plot: bool = True,
) -> None:
    experiment = engine.build_experiment(file_path, options=options, label=file_path.stem)
    experiment.load()

    # Fit the first-principles Pyomo model to the chosen file.
    experiment.fit()

    summary = experiment.summarize()

    print("\nFinished:")
    print(f"  file: {file_path}")
    print(f"  dataset: {summary['dataset']}")
    if experiment.fit_stru and 'parameters' in experiment.fit_stru:
        print(f"  parameters: {experiment.fit_stru['parameters']}")
    if run_fim:
        # Compute sensitivity / Fisher information after the fit is available.
        fim_result = experiment.calc_fim()
        print(f"  FIM trace: {fim_result.get('trace')}")
        print(f"  FIM determinant: {fim_result.get('det')}")
        if fim_result.get("std") is not None:
            print(f"  parameter std: {fim_result.get('std')}")
        if run_suggest:
            suggestion_result = _suggest_new_experiments(experiment, fim_result, options.output_dir)
            print(f"  suggestion file: {suggestion_result.get('path')}")
    if run_plot:
        # Plot the matched model vs. data.
        plot_result = experiment.plot()
        print(f"  output: {plot_result}")


def _suggest_new_experiments(
    experiment,
    fim_result: dict,
    output_dir: Path,
) -> dict:
    """Write a short, human-readable next-step summary from the UQ results."""
    lines = [
        "Suggested next experiments",
        "--------------------------",
        f"Dataset: {experiment.dataset}",
        f"Label: {experiment.label}",
        f"Model family: {experiment.options.model_family.value}",
        "",
    ]

    std = fim_result.get("std") or []
    eig_dir = fim_result.get("eig_dir") or []
    if std:
        ranked = sorted(enumerate(std), key=lambda item: item[1], reverse=True)
        lines.append("Parameters with the largest estimated uncertainty:")
        for idx, value in ranked[:3]:
            name = eig_dir[idx] if idx < len(eig_dir) else f"parameter_{idx + 1}"
            lines.append(f"  - {name}: std = {value}")
        lines.append("")
    else:
        lines.append("No parameter standard deviations were returned by the FIM step.")
        lines.append("")

    lines.extend(
        [
            "Simple next-step ideas:",
            "  - Add data in the region that makes the uncertain parameters move the most.",
            "  - Repeat the fit after adding new conductivity, concentration, or flow measurements.",
            "  - Use the largest-uncertainty directions to pick the next test condition.",
        ]
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{experiment.label}-suggestions.txt"
    out_path.write_text("\n".join(lines) + "\n")
    return {"path": str(out_path), "lines": lines}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the simple diafiltration workflow.")
    parser.add_argument(
        "--paper",
        choices=("DATA1", "DATA2"),
        default=None,
        help="Shortcut to recreate the DATA1 or DATA2 paper figures.",
    )
    parser.add_argument(
        "--file",
        dest="files",
        action="append",
        default=[],
        help="Input file path. Use more than once for multiple files.",
    )
    parser.add_argument(
        "--mode",
        default=None,
        help="Model mode to use. DATA1 usually uses DATA; DATA2 usually uses Lag.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Folder where plots and outputs will be saved.",
    )
    parser.add_argument(
        "--run-fim",
        action="store_true",
        help="Also compute the FIM after fitting.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    loader = DataLoader(repo_root=REPO_ROOT, data_root=REPO_ROOT)
    paper_family = args.paper or _paper_workflow_mode()
    if paper_family in {"DATA1", "DATA2"}:
        output_dir_text = args.output_dir or str(REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "reproduction")
        _run_paper_shortcut(paper_family, Path(output_dir_text).expanduser().resolve())
        return

    engine = UQEngine(loader)
    family = _choose_workflow_profile()

    if args.files:
        file_paths = _parse_file_args(args.files)
    else:
        file_paths = _choose_input_files(loader, family)

    if not file_paths:
        raise ValueError("No input files were selected.")

    default_mode = "DATA" if family == WorkflowFamily.DATA1 else "Lag"
    mode = args.mode or _prompt_text("Model mode", default_mode)
    output_dir_text = args.output_dir or _prompt_text("Output folder", str(REPO_ROOT / "figures"))
    output_dir = Path(output_dir_text).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_depth = _choose_run_depth()
    do_fim = args.run_fim or run_depth in {"uq", "suggest", "full"}
    do_suggest = run_depth in {"suggest", "full"}
    do_plot = run_depth == "full"

    print("\nSummary of your choices:")
    print(f"  family: {family.value}")
    print(f"  mode: {mode}")
    print(f"  output folder: {output_dir}")
    print(f"  run depth: {run_depth}")
    print(f"  compute FIM: {do_fim}")
    print(f"  suggest new experiments: {do_suggest}")
    print(f"  plot results: {do_plot}")
    print("  files:")
    for file_path in file_paths:
        print(f"    - {file_path}")

    # Build one shared options object, then run each chosen file the same way.
    options = _make_options(family, mode, output_dir)
    for file_path in file_paths:
        converted = _maybe_convert_conductivity(loader, file_path, output_dir)
        if converted is not None:
            continue
        _run_one(
            engine,
            file_path,
            options,
            run_fim=do_fim,
            run_suggest=do_suggest,
            run_plot=do_plot,
        )


if __name__ == "__main__":
    main()
