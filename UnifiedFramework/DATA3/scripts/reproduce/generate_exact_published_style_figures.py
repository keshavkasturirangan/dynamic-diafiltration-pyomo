#!/usr/bin/env python3
"""Generate exact published-style DATA1/DATA2 comparison figures from the unified codebase."""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from UnifiedFramework.DATA3.scripts.reproduce.reproduce_data1_data2 import (
    CANONICAL_DATASETS,
    _build_guess_from_theta,
    build_model_options,
)
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library import RunMode
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.workflow.output_artifacts import extract_model_trajectories
from UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.workflow.unified_workflow import (
    DataLoader,
    DatasetRequest,
    DiafiltrationExperiment,
)


SIM_ROOT = REPO_ROOT / "UnifiedFramework" / "DATA3" / "docs" / "validation" / "simulation_validation_data_files"
DATA1_DATA = REPO_ROOT / "DATA1_matlab" / "data"
OUT_DIR = REPO_ROOT / "UnifiedFramework" / "DATA3" / "results" / "reproduction" / "20260408-codex-paper-check" / "figures_exact"


SIGMA_STYLES = {
    0.1: {"color": "#FF4D4D", "linestyle": "--"},
    0.5: {"color": "#1F4AFF", "linestyle": "-"},
    0.9: {"color": "#5CB85C", "linestyle": ":"},
}


def _style_axes(ax) -> None:
    ax.tick_params(direction="in", top=True, right=True, labelsize=11)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def _label_panel(ax, label: str) -> None:
    ax.text(-0.18, 1.05, label, transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")


def _render_data1_img004() -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    panel_specs = [
        ("A", SIM_ROOT / "data1_main" / "fig2" / "data1_main_fig2_a.csv", "Mass in Vial [g]", "mass_g"),
        ("B", SIM_ROOT / "data1_main" / "fig2" / "data1_main_fig2_b.csv", "Concentration [mM]", "concentration_mM"),
        ("C", SIM_ROOT / "data1_main" / "fig2" / "data1_main_fig2_c.csv", "Mass in Vial [g]", "mass_g"),
        ("D", SIM_ROOT / "data1_main" / "fig2" / "data1_main_fig2_d.csv", "Concentration [mM]", "concentration_mM"),
    ]

    for ax, (label, csv_path, ylabel, value_col) in zip(axes.flat, panel_specs):
        df = pd.read_csv(csv_path)
        if value_col == "mass_g":
            for series_id, series in df.groupby("series_id"):
                series = series.sort_values("time_min")
                if "measurement" in series_id:
                    ax.plot(series["time_min"], series[value_col], "o", color="red", markersize=2.5)
                else:
                    ax.plot(series["time_min"], series[value_col], color="#4C63FF", linewidth=2.0, alpha=0.9)
        else:
            for series_id, series in df.groupby("series_id"):
                series = series.sort_values("time_min")
                if "retentate_prediction" in series_id:
                    ax.plot(series["time_min"], series[value_col], color="green", linewidth=2.0)
                elif "permeate_prediction" in series_id:
                    ax.plot(series["time_min"], series[value_col], color="#FF5A5A", linewidth=2.0)
                elif "vial_prediction" in series_id:
                    ax.plot(series["time_min"], series[value_col], "^", color="#FF5A5A", markersize=5)
                elif "retentate_measurement" in series_id:
                    ax.plot(series["time_min"], series[value_col], "s", color="#CC33CC", markersize=4)
                elif "permeate_measurement" in series_id:
                    ax.plot(series["time_min"], series[value_col], "s", color="#22C1C3", markersize=4)

        ax.set_xlabel("Time [min]", fontsize=15, fontweight="bold")
        ax.set_ylabel(ylabel, fontsize=15, fontweight="bold")
        ax.set_ylim(bottom=0)
        _style_axes(ax)
        _label_panel(ax, label)

    fig.tight_layout()
    out = OUT_DIR / "data1_img004_unified_exact.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def _render_data1_img006() -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    panel_rows = [
        ("A", "501.1", "a"),
        ("B", "511.12", "b"),
    ]
    cols = [
        ("mass", "Mass [g]", "mass_g"),
        ("permeate", "Permeate Conc. [mM]", "permeate_concentration_mM"),
        ("retentate", "Retentate Conc. [mM]", "retentate_concentration_mM"),
    ]

    for row_idx, (panel_label, dataset, panel_key) in enumerate(panel_rows):
        for col_idx, (name, ylabel, value_col) in enumerate(cols):
            ax = axes[row_idx, col_idx]
            df = pd.read_csv(SIM_ROOT / "data1_main" / "fig4" / f"data1_main_fig4_{panel_key}_{name}.csv")
            for sigma in sorted(df["sigma"].unique()):
                series = df[df["sigma"] == sigma].sort_values(["series_id", "time_min"])
                style = SIGMA_STYLES[round(float(sigma), 1)]
                for _, subgroup in series.groupby("series_id"):
                    ax.plot(
                        subgroup["time_min"],
                        subgroup[value_col],
                        color=style["color"],
                        linestyle=style["linestyle"],
                        linewidth=2.0,
                    )
            if row_idx == 0:
                ax.set_title(["Mass", "Permeate", "Retentate"][col_idx], fontsize=17, fontweight="bold")
            if col_idx == 0:
                _label_panel(ax, panel_label)
            ax.set_xlabel("Time [min]", fontsize=14, fontweight="bold")
            ax.set_ylabel(ylabel, fontsize=14, fontweight="bold")
            ax.set_ylim(bottom=0)
            _style_axes(ax)

    handles = [
        plt.Line2D([0], [0], color=SIGMA_STYLES[s]["color"], linestyle=SIGMA_STYLES[s]["linestyle"], linewidth=2.0)
        for s in (0.1, 0.5, 0.9)
    ]
    labels = [r"$\sigma$ = 0.1", r"$\sigma$ = 0.5", r"$\sigma$ = 0.9"]
    axes[0, 2].legend(handles, labels, fontsize=13, loc="upper right")

    fig.tight_layout()
    out = OUT_DIR / "data1_img006_unified_exact.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def _plot_contour_grid(ax, df: pd.DataFrame, x_col: str, value_col: str, x_label: str, y_label: str) -> None:
    x = np.reshape(df[x_col].to_numpy(dtype=float), (50, 50))
    y = np.reshape(df["Lp"].to_numpy(dtype=float), (50, 50))
    z = np.reshape(df[value_col].to_numpy(dtype=float), (50, 50))
    levels = np.linspace(float(np.nanmin(z)), float(np.nanmax(z)), 10)
    cs = ax.contour(x, y, z, levels=levels, linewidths=1.6)
    ax.clabel(cs, cs.levels[::2], inline=True, fontsize=9, colors="k", fmt="%1.1f")
    idx = int(np.nanargmin(df[value_col].to_numpy(dtype=float)))
    ax.plot(
        float(df.iloc[idx][x_col]),
        float(df.iloc[idx]["Lp"]),
        "^",
        markersize=8,
        markeredgecolor="red",
        markerfacecolor=(1.0, 0.7, 0.7),
        clip_on=False,
    )
    ax.set_xlabel(x_label, fontsize=12, fontweight="bold")
    ax.set_ylabel(y_label, fontsize=12, fontweight="bold")
    _style_axes(ax)


def _render_data1_contour_exact(out_name: str, x_kind: str) -> Path:
    fig, axes = plt.subplots(3, 3, figsize=(13, 13))
    datasets = [
        # These paper panels come from the diafiltration family, not the filtration
        # 501.x cases. Using the diafiltration rows restores the published regions.
        ("A", "511.12 concpolar"),
        ("B", "511.11 concpolar"),
        ("C", "511.12"),
    ]
    cols = [
        ("Obj_mass", "Mass"),
        ("Obj_concentration", "Permeate"),
        ("Obj_retentate_concentration", "Retentate"),
    ]
    suffix = "x_sigma-y_Lp.csv" if x_kind == "sigma" else "x_B-y_Lp.csv"
    x_col = "sigma" if x_kind == "sigma" else "B"
    x_label = r"$\sigma$ [dimensionless]" if x_kind == "sigma" else r"B [$\mu$m $\cdot$ s$^{-1}$]"
    y_label = r"L$_p$ [L $\cdot$ m$^{-2}$ $\cdot$ h$^{-1}$ $\cdot$ bar$^{-1}$]"
    x_limits = (0.0, 1.0) if x_kind == "sigma" else (0.0, 2.0)
    y_limits = (0.5, 7.4)

    for row_idx, (panel_label, dataset_dir) in enumerate(datasets):
        df = pd.read_csv(DATA1_DATA / dataset_dir / f"contourdata-{suffix}")
        for col_idx, (value_col, title) in enumerate(cols):
            ax = axes[row_idx, col_idx]
            _plot_contour_grid(ax, df, x_col, value_col, x_label, y_label)
            ax.set_xlim(*x_limits)
            ax.set_ylim(*y_limits)
            if row_idx == 0:
                ax.set_title(title, fontsize=17, fontweight="bold")
            if col_idx == 0:
                _label_panel(ax, panel_label)

    fig.tight_layout()
    out = OUT_DIR / out_name
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def _render_data2_img007() -> Path | None:
    loader = DataLoader()
    datasets = ["DATA2_270511.123", "DATA2_270611.123"]
    sims: dict[str, tuple[object, dict[int, dict[str, np.ndarray]]]] = {}

    for dataset_id in datasets:
        spec = next(spec for spec in CANONICAL_DATASETS if spec.dataset_id == dataset_id)
        exp, _ = loader.load(DatasetRequest(file_path=str(REPO_ROOT / spec.data_file), dataset_id=dataset_id))
        theta = dict(spec.fixed_guess or spec.initial_guess or {})
        opts = build_model_options(spec, nfe=10, parameter_mode=RunMode.SIMULATION)
        guess = _build_guess_from_theta(theta, opts)
        try:
            model = DiafiltrationExperiment(exp, opts, guess=guess).build_simulation_model(
                theta_values=theta,
                solver_name="ipopt",
                solver_options={"max_iter": 1000, "tol": 1e-6, "acceptable_tol": 1e-6},
                tee=False,
            )
        except Exception:
            return None
        sims[dataset_id] = (exp, extract_model_trajectories(exp, model))

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    row_specs = [("A", "DATA2_270511.123"), ("D", "DATA2_270611.123")]
    row_specs += [("B", "DATA2_270511.123"), ("E", "DATA2_270611.123")]
    # The published layout is A/B/C and D/E/F; we map them below by panel position.
    panels = [
        ("A", "DATA2_270511.123", "mass"),
        ("B", "DATA2_270511.123", "conc"),
        ("C", "DATA2_270511.123", "stirred"),
        ("D", "DATA2_270611.123", "mass"),
        ("E", "DATA2_270611.123", "conc"),
        ("F", "DATA2_270611.123", "stirred"),
    ]
    axes_flat = axes.flat
    for ax, (panel_label, dataset_id, mode) in zip(axes_flat, panels):
        exp, sim = sims[dataset_id]
        t_delay = float(exp.vials[0].time_s[0])
        for i, vial in enumerate(exp.vials, start=1):
            t_meas = (np.asarray(vial.time_s, dtype=float) - t_delay) / 60.0
            t_pred = sim[i]["time_s"] / 60.0
            if mode == "mass":
                if vial.mass_g is not None:
                    ax.plot(t_meas, np.asarray(vial.mass_g, dtype=float), "o", color="red", markersize=2.5)
                ax.plot(t_pred, sim[i]["mV"], color="#4C63FF", linewidth=2.0)
                ax.set_ylabel("Mass in Vial [g]", fontsize=13, fontweight="bold")
            elif mode == "conc":
                cv_vals = np.asarray(np.atleast_1d(vial.cV_avg), dtype=float) if vial.cV_avg is not None else np.array([])
                cv_vals = cv_vals[np.isfinite(cv_vals)]
                if cv_vals.size == 1:
                    ax.plot([(float(vial.time_s[-1]) - t_delay) / 60.0], [float(cv_vals[0])], "s", color="#22C1C3", markersize=4)
                elif cv_vals.size > 1:
                    cv_times = np.asarray(np.atleast_1d(vial.time_s), dtype=float)[-cv_vals.size:]
                    ax.plot((cv_times - t_delay) / 60.0, cv_vals, "s", color="#22C1C3", markersize=4)
                if vial.retentate_signal is not None:
                    cf_meas = np.asarray(np.atleast_1d(vial.retentate_signal), dtype=float)
                    mask = np.isfinite(cf_meas)
                    ax.plot(t_meas[mask], cf_meas[mask], "s", color="#CC33CC", markersize=4)
                ax.plot(t_pred, sim[i]["cF"], color="green", linewidth=2.0)
                ax.plot(t_pred, sim[i]["cH"], color="#FF5A5A", linewidth=2.0)
                ax.plot([t_pred[-1]], [sim[i]["cV"][-1]], "^", color="#FF5A5A", markersize=5)
                ax.set_ylabel("Concentration [mM]", fontsize=13, fontweight="bold")
            else:
                if "mF" in sim[i]:
                    ax.plot(t_pred, sim[i]["mF"], color="#4C63FF", linewidth=2.0)
                ax.set_ylabel("Mass in Stirred Cell [g]", fontsize=13, fontweight="bold")
            ax.set_xlabel("Time [min]", fontsize=13, fontweight="bold")
            ax.set_ylim(bottom=0)
            _style_axes(ax)
        _label_panel(ax, panel_label)

    fig.tight_layout()
    out = OUT_DIR / "data2_img007_unified_exact.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    generated = [
        _render_data1_img004(),
        _render_data1_img006(),
        _render_data1_contour_exact("data1_img007_unified_exact.png", x_kind="sigma"),
        _render_data1_contour_exact("data1_img008_unified_exact.png", x_kind="B"),
    ]
    data2 = _render_data2_img007()
    if data2 is not None:
        generated.append(data2)
    for path in generated:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
