#!/usr/bin/env python3
"""Simple user-facing runner for the refactored diafiltration workflow.

The goal of this file is to give non-coders a small set of easy choices:

1. Choose DATA1 or DATA2
2. Choose one or more input files
3. Choose an output folder
4. Run the fit and plotting workflow

All scientific logic stays in `refactored_ucb_library.py`.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Sequence

import pandas as pd

try:
    from .refactored_ucb_library import DataLoader, ModelOptions, UQEngine, WorkflowFamily
except ImportError:  # pragma: no cover - allows running as a plain script
    from refactored_ucb_library import DataLoader, ModelOptions, UQEngine, WorkflowFamily


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


def _run_one(engine: UQEngine, file_path: Path, options: ModelOptions) -> None:
    experiment = engine.build_experiment(file_path, options=options, label=file_path.stem)
    experiment.load()
    experiment.fit()
    plot_result = experiment.plot()
    summary = experiment.summarize()
    print("\nFinished:")
    print(f"  file: {file_path}")
    print(f"  dataset: {summary['dataset']}")
    if experiment.fit_stru and 'parameters' in experiment.fit_stru:
        print(f"  parameters: {experiment.fit_stru['parameters']}")
    print(f"  output: {plot_result}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the simple diafiltration workflow.")
    parser.add_argument(
        "--family",
        choices=("DATA1", "DATA2"),
        default=None,
        help="Workflow family to run. If omitted, you will be prompted.",
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
        "--fim",
        action="store_true",
        help="Also compute the FIM after fitting.",
    )
    parser.add_argument(
        "--convert-conductivity",
        action="store_true",
        help="Convert a conductivity CSV/XLSX to concentration and save a new CSV.",
    )
    parser.add_argument(
        "--conductivity-file",
        default=None,
        help="CSV/XLSX file to convert when --convert-conductivity is used.",
    )
    parser.add_argument(
        "--conductivity-column",
        default="Conductivity",
        help="Column name that stores conductivity values in the tabular file.",
    )
    parser.add_argument(
        "--conductivity-model",
        choices=("variant_shedlovsky", "msa"),
        default="variant_shedlovsky",
        help="Paper model used for conductivity conversion.",
    )
    parser.add_argument(
        "--salt-name",
        default="NaCl",
        help="Salt name used for paper-model parameter autofill.",
    )
    parser.add_argument(
        "--output-units",
        choices=("mM", "M"),
        default="mM",
        help="Units for the converted concentration column.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    loader = DataLoader(repo_root=REPO_ROOT, data_root=REPO_ROOT)
    engine = UQEngine(loader)

    if args.convert_conductivity:
        conduct_file = args.conductivity_file
        if not conduct_file:
            # If the user does not type a path, offer a sensible default.
            conduct_file = _prompt_text("Conductivity file to convert", str(REPO_ROOT / "UnifiedFramework" / "ExperimentalDataFiles" / "NF270_MC2.xlsx"))
        conduct_path = Path(conduct_file).expanduser().resolve()
        output_dir_text = args.output_dir or _prompt_text("Output folder", str(REPO_ROOT / "figures"))
        output_dir = Path(output_dir_text).expanduser().resolve()
        # In conversion mode we stop after writing the new concentration CSV.
        _convert_conductivity_file(
            loader,
            conduct_path,
            output_dir=output_dir,
            column=args.conductivity_column,
            model=args.conductivity_model,
            salt_name=args.salt_name,
            output_units=args.output_units,
        )
        return

    family_name = args.family or _prompt_choice("Choose workflow family", ("DATA1", "DATA2"), "DATA1")
    family = WorkflowFamily(family_name)

    if args.files:
        file_paths = _parse_file_args(args.files)
    elif family == WorkflowFamily.DATA1:
        # DATA1 users can pick from a short list of known MAT files.
        print("\nDATA1 preset files are available, but you can still type your own paths.")
        file_paths = _select_files_from_list(_default_data1_files(loader))
    else:
        # DATA2 users usually paste the exact file path they want to analyze.
        raw = input("\nPaste one or more DATA2 file paths separated by commas: ").strip()
        if not raw:
            raise ValueError("No input files were provided.")
        file_paths = _parse_file_args([token.strip() for token in raw.split(",") if token.strip()])

    if not file_paths:
        raise ValueError("No input files were selected.")

    default_mode = "DATA" if family == WorkflowFamily.DATA1 else "Lag"
    mode = args.mode or _prompt_text("Model mode", default_mode)
    output_dir_text = args.output_dir or _prompt_text("Output folder", str(REPO_ROOT / "figures"))
    output_dir = Path(output_dir_text).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_fim = args.fim or _prompt_yes_no("Compute FIM too?", default=False)

    print("\nSummary of your choices:")
    print(f"  family: {family.value}")
    print(f"  mode: {mode}")
    print(f"  output folder: {output_dir}")
    print(f"  compute FIM: {run_fim}")
    print("  files:")
    for file_path in file_paths:
        print(f"    - {file_path}")

    options = _make_options(family, mode, output_dir)
    for file_path in file_paths:
        _run_one(engine, file_path, options)

        if run_fim:
            experiment = engine.build_experiment(file_path, options=options, label=file_path.stem)
            experiment.load()
            experiment.fit()
            fim_result = experiment.calc_fim()
            print(f"  FIM trace: {fim_result.get('trace')}")
            print(f"  FIM determinant: {fim_result.get('det')}")


if __name__ == "__main__":
    main()
