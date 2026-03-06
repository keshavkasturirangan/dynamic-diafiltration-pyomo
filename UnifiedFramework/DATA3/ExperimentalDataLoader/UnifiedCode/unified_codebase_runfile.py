#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""User-facing runner with explicit DATA1, DATA2, and DATA3 branches."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional

try:
    from .unified_codebase_library import (
        ExperimentMode,
        ModelOptions,
        RunMode,
        UnifiedPipelineConfigV24,
        run_unified_pipeline_v24,
    )
except ImportError:
    from unified_codebase_library import (
        ExperimentMode,
        ModelOptions,
        RunMode,
        UnifiedPipelineConfigV24,
        run_unified_pipeline_v24,
    )

REPO_ROOT = Path(__file__).resolve().parents[4]

PRESET_DATA1_PATHS = [
    REPO_ROOT / "DATA1_matlab" / "data_library" / "data_stru-dataset501.1.mat",
    REPO_ROOT / "DATA1_matlab" / "data_library" / "data_stru-dataset501.11.mat",
    REPO_ROOT / "DATA1_matlab" / "data_library" / "data_stru-dataset511.11.mat",
    REPO_ROOT / "DATA1_matlab" / "data_library" / "data_stru-dataset511.12.mat",
]
PRESET_DATA2_PATHS = [
    REPO_ROOT / "DATA1_matlab" / "data_library" / "data_stru-dataset270511.123.mat",
    REPO_ROOT / "DATA1_matlab" / "data_library" / "data_stru-dataset270611.123.mat",
]
PRESET_DATA3_PATH = REPO_ROOT / "UnifiedFramework" / "ExperimentalDataFiles" / "NF270_MC2.xlsx"
PRESET_DATA3_SELECTOR = "05.07.24_NaCl"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run unified pipeline with a preset profile.\n"
            "Profiles: DATA1, DATA2, DATA3"
        )
    )
    parser.add_argument(
        "--profile",
        choices=("DATA1", "DATA2", "DATA3"),
        default="DATA1",
        help="Preset run profile.",
    )
    parser.add_argument(
        "--file-path",
        default=None,
        help="Optional override input file path. Recommended for DATA3 custom experiments.",
    )
    parser.add_argument(
        "--selector",
        default=None,
        help="Optional XLSX sheet selector override.",
    )
    parser.add_argument(
        "--run-doe",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Enable/disable DoE/FIM stage. Defaults: DATA1/DATA2=False, DATA3=True.",
    )
    parser.add_argument(
        "--calc-cov",
        action="store_true",
        help="Enable covariance computation in ParmEst.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Enable quick plotting during load stage.",
    )
    parser.add_argument(
        "--convert-to-concentration",
        action="store_true",
        help="Convert conductivity to concentration during load stage.",
    )
    parser.add_argument(
        "--nfe",
        type=int,
        default=30,
        help="Finite elements for discretization.",
    )
    parser.add_argument(
        "--solver",
        default="ipopt",
        help="Solver name (default: ipopt).",
    )
    return parser.parse_args()


def _base_specs() -> Dict[str, object]:
    return {
        "mode": "Lag",
        "delP_bar": 10.0,
        "Temp_K": 298.15,
        "Am_cm2": 4.1,
        "rho_g_cm3": 1.0,
        "C_D_value": 50.0,
        "C_D_units": "mM",
        "M_F0_g": 200.0,
    }


def _resolve_run_doe(profile: str, cli_value: Optional[bool]) -> bool:
    """Resolve run_doe with profile-aware defaults and CLI override."""
    defaults = {
        "DATA1": False,
        "DATA2": False,
        "DATA3": True,
    }
    if cli_value is not None:
        return bool(cli_value)
    return bool(defaults[profile])


def _build_data1_config(file_path: Path, args: argparse.Namespace, run_doe: bool) -> UnifiedPipelineConfigV24:
    """Build DATA1 reproduction config for one MAT file."""
    model_options = ModelOptions(
        mode=ExperimentMode.DATA,
        run_mode=RunMode.ESTIMATION,
        b_form="single",
        nfe=int(args.nfe),
        use_sigma_logit_transform=False,
        use_multistart_for_mat_legacy=True,
        multistart_iterations=30,
    )
    return UnifiedPipelineConfigV24(
        file_path=str(file_path),
        selector=args.selector,
        specs=None,
        convert_to_concentration=bool(args.convert_to_concentration),
        plot=bool(args.plot),
        model_options=model_options,
        run_parmest=True,
        run_doe=run_doe,
        calc_cov=bool(args.calc_cov),
        cov_n=None,
        solver=args.solver,
        solver_options=None,
        tee=False,
    )


def _build_data2_config(file_path: Path, args: argparse.Namespace, run_doe: bool) -> UnifiedPipelineConfigV24:
    """Build canonical DATA2 reproduction config for one MAT file."""
    model_options = ModelOptions(
        mode=ExperimentMode.LAG,
        run_mode=RunMode.ESTIMATION,
        b_form="convection",
        nfe=int(args.nfe),
        use_sigma_logit_transform=False,
        fix_sigma_in_estimation=True,
    )
    return UnifiedPipelineConfigV24(
        file_path=str(file_path),
        selector=args.selector,
        specs=None,
        convert_to_concentration=bool(args.convert_to_concentration),
        plot=bool(args.plot),
        model_options=model_options,
        run_parmest=True,
        run_doe=run_doe,
        calc_cov=bool(args.calc_cov),
        cov_n=None,
        solver=args.solver,
        solver_options=None,
        tee=False,
    )


def _build_data3_config(args: argparse.Namespace, run_doe: bool) -> UnifiedPipelineConfigV24:
    """Build DATA3/new-experiment config (XLSX-centric path)."""
    file_path = Path(args.file_path).expanduser().resolve() if args.file_path else PRESET_DATA3_PATH
    selector = args.selector if args.selector is not None else PRESET_DATA3_SELECTOR
    model_options = ModelOptions(
        mode=ExperimentMode.DATA,
        run_mode=RunMode.ESTIMATION,
        b_form="single",
        nfe=int(args.nfe),
    )
    return UnifiedPipelineConfigV24(
        file_path=str(file_path),
        selector=selector,
        specs=_base_specs(),
        convert_to_concentration=bool(args.convert_to_concentration),
        plot=bool(args.plot),
        model_options=model_options,
        run_parmest=True,
        run_doe=run_doe,
        calc_cov=bool(args.calc_cov),
        cov_n=None,
        solver=args.solver,
        solver_options=None,
        tee=False,
    )


def _print_result(result: Dict[str, object]) -> None:
    print("Load OK:", result["load_ok"])
    for sev, msg in result["issues"]:
        print(f"{sev}: {msg}")

    print("Theta names:", result["theta_names"])
    print("Labeled outputs:", result["labeled_counts"]["outputs"])
    print("Labeled unknowns:", result["labeled_counts"]["unknowns"])
    print("Labeled design inputs:", result["labeled_counts"]["inputs"])

    if "parmest" in result:
        est = result["parmest"]
        if "warning" in est:
            print("WARNING:", est["warning"])
        if "covariance_warning" in est:
            print("WARNING:", est["covariance_warning"])
        print("ParmEst objective:", est.get("objective"))
        print("ParmEst theta:\n", est.get("theta"))
        if "covariance" in est:
            print("ParmEst covariance method:", est.get("covariance_method"))
            print("ParmEst covariance:\n", est.get("covariance"))

    if "doe" in result:
        doe = result["doe"]
        print("DoE FIM:\n", doe["fim"])
        print("DoE D-opt (log-det):", doe["d_opt_logdet"])


def main() -> None:
    args = _parse_args()
    run_doe = _resolve_run_doe(args.profile, args.run_doe)

    if args.profile == "DATA1":
        print("Profile: DATA1")
        print("Purpose: Reproduce canonical DATA1 baselines")
        if args.file_path:
            file_paths = [Path(args.file_path).expanduser().resolve()]
        else:
            file_paths = PRESET_DATA1_PATHS
        for i, file_path in enumerate(file_paths, start=1):
            config = _build_data1_config(file_path=file_path, args=args, run_doe=run_doe)
            print(f"\nDATA1 case [{i}/{len(file_paths)}]")
            print(f"Input file: {config.file_path}")
            result = run_unified_pipeline_v24(config)
            _print_result(result)
        return

    if args.profile == "DATA2":
        if args.file_path:
            file_paths = [Path(args.file_path).expanduser().resolve()]
        else:
            file_paths = PRESET_DATA2_PATHS
        print("Profile: DATA2")
        print("Purpose: Reproduce canonical DATA2 baselines")
        for i, file_path in enumerate(file_paths, start=1):
            config = _build_data2_config(file_path=file_path, args=args, run_doe=run_doe)
            print(f"\nDATA2 case [{i}/{len(file_paths)}]")
            print(f"Input file: {config.file_path}")
            result = run_unified_pipeline_v24(config)
            _print_result(result)
        return

    if args.profile == "DATA3":
        config = _build_data3_config(args, run_doe=run_doe)
        print("Profile: DATA3")
        print("Purpose: Analyze NEW/DATA3 experiments")
        print(f"Input file: {config.file_path}")
        if config.selector is not None:
            print(f"Selector: {config.selector}")
        result = run_unified_pipeline_v24(config)
        _print_result(result)
        return

    raise ValueError(f"Unknown profile: {args.profile}")


if __name__ == "__main__":
    main()
