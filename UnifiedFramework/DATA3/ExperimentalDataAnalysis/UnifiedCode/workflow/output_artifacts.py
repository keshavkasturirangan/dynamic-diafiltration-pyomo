"""Workflow-owned plotting and validation artifact builders."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


def extract_model_trajectories(exp, model) -> Dict[int, Dict[str, np.ndarray]]:
    """Extract per-vial trajectories in physical time from a discretized model."""
    import pyomo.environ as pyo

    tau_vals = sorted(float(t) for t in list(model.tau))
    t_delay = float(exp.vials[0].time_s[0])
    out: Dict[int, Dict[str, np.ndarray]] = {}
    for n in range(1, len(exp.vials) + 1):
        vial = exp.vials[n - 1]
        ti = float(vial.time_s[0]) - t_delay
        tf = float(vial.time_s[-1]) - t_delay
        dur = tf - ti
        t_model = np.array([ti + dur * tau for tau in tau_vals], dtype=float)
        mV = np.array([float(pyo.value(model.mV[n, tau])) for tau in tau_vals], dtype=float)
        cF = np.array([float(pyo.value(model.cF[n, tau])) for tau in tau_vals], dtype=float)
        cH = np.array([float(pyo.value(model.cH[n, tau])) for tau in tau_vals], dtype=float)
        cV = np.array([float(pyo.value(model.cV[n, tau])) for tau in tau_vals], dtype=float)
        Jw = np.array([float(pyo.value(model.Jw[n, tau])) for tau in tau_vals], dtype=float)
        Js = np.array([float(pyo.value(model.Js[n, tau])) for tau in tau_vals], dtype=float)
        out[n] = {"time_s": t_model, "mV": mV, "cF": cF, "cH": cH, "cV": cV, "Jw": Jw, "Js": Js}
    return out


def legacy_style_objective(exp, sim: Dict[int, Dict[str, np.ndarray]], *, t_delay: float) -> float:
    """Approximate legacy utility.py objective on unified simulated trajectories."""
    n_vials = len(exp.vials)
    n_extra = int(getattr(exp, "n_extra", 0) or 0)
    n_v0 = int(getattr(exp, "n_v0", 1) or 1)
    collect_vial = max(1, n_vials - n_extra)

    obj_m = 0.0
    obj_cp = 0.0
    obj_cf0 = 0.0
    obj_cf = 0.0
    cnt_m = 0
    cnt_cp = 0
    cnt_cf = 0
    cnt_cf0 = 0

    for n, vial in enumerate(exp.vials, start=1):
        t_meas_s = np.asarray(vial.time_s, dtype=float)

        if vial.mass_g is not None:
            y_m = np.asarray(vial.mass_g, dtype=float)
            y_m_pred = np.interp(t_meas_s, sim[n]["time_s"] + t_delay, sim[n]["mV"])
            mask = np.isfinite(y_m)
            if mask.any() and n >= n_v0:
                res = (y_m_pred[mask] - y_m[mask]) / 0.01
                obj_m += float(np.sum(res**2))
                cnt_m += int(np.sum(mask))

        y_cv = float(vial.cV_avg) if isinstance(vial.cV_avg, (int, float, np.number)) else np.nan
        if np.isfinite(y_cv) and n > n_extra and abs(y_cv) > 1e-12:
            y_cv_pred = float(np.interp(float(vial.time_s[-1]), sim[n]["time_s"] + t_delay, sim[n]["cV"]))
            obj_cp += float(((y_cv_pred - y_cv) / (0.03 * y_cv)) ** 2)
            cnt_cp += 1

        if vial.retentate_signal is not None:
            y_cf = np.asarray(vial.retentate_signal, dtype=float)
            y_cf_pred = np.interp(t_meas_s, sim[n]["time_s"] + t_delay, sim[n]["cF"])
            mask = np.isfinite(y_cf) & (np.abs(y_cf) > 1e-12)
            if mask.any():
                res = (y_cf_pred[mask] - y_cf[mask]) / (0.003 * y_cf[mask])
                if n >= n_v0:
                    obj_cf += float(np.sum(res**2))
                    cnt_cf += int(np.sum(mask))
                else:
                    obj_cf0 += float(np.sum(res**2))
                    cnt_cf0 += int(np.sum(mask))

    term_m = obj_m / max(1, cnt_m)
    term_cp = obj_cp / max(1, cnt_cp if cnt_cp > 0 else collect_vial)
    term_cf = (obj_cf0 + obj_cf) / max(1, cnt_cf0 + cnt_cf)
    return float(1e4 * (term_m + term_cp + term_cf))


def build_data1_stage_a_outputs(
    *,
    out_dir: Path,
    exp,
    est: Dict[str, object],
    sim: Dict[int, Dict[str, np.ndarray]],
    tolerances,
    references: pd.DataFrame,
) -> Dict[str, object]:
    """Generate DATA1 Stage A concrete outputs."""
    import matplotlib.pyplot as plt

    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    t_delay = float(exp.vials[0].time_s[0])
    n_vials = len(exp.vials)
    theta = flatten_theta(est.get("theta"))
    lp = float(theta.get("Lp", np.nan))
    b = float(theta.get("B", np.nan))
    sigma = float(theta.get("sigma", np.nan))

    mass_true: List[float] = []
    mass_pred: List[float] = []
    cf_true: List[float] = []
    cf_pred: List[float] = []
    cv_true: List[float] = []
    cv_pred: List[float] = []

    fig2a = plt.figure(figsize=(4, 4))
    for i, v in enumerate(exp.vials, start=1):
        t = (np.asarray(v.time_s, dtype=float) - t_delay) / 60.0
        y = np.asarray(v.mass_g, dtype=float) if v.mass_g is not None else np.array([])
        if t.size and y.size:
            plt.plot(t, y, "r.", markersize=4, label="Measured" if i == 1 else None)
    plt.xlabel("Time [min]")
    plt.ylabel("Mass in Vial [g]")
    plt.ylim(bottom=0)
    plt.tick_params(direction="in")
    plt.legend(loc="best")
    fig2a_path = figures_dir / "data1_fig2_a_mass_measured.png"
    fig2a.savefig(fig2a_path, bbox_inches="tight")
    plt.close(fig2a)

    fig2b = plt.figure(figsize=(4, 4))
    for i, v in enumerate(exp.vials, start=1):
        t = (np.asarray(v.time_s, dtype=float) - t_delay) / 60.0
        y = np.asarray(v.retentate_signal, dtype=float) if v.retentate_signal is not None else np.array([])
        if t.size and y.size:
            plt.plot(t, y, "ms", markersize=4, label="Measured cF" if i == 1 else None)
    plt.xlabel("Time [min]")
    plt.ylabel("Retentate concentration [mM]")
    plt.ylim(bottom=0)
    plt.tick_params(direction="in", top=True, right=True)
    plt.legend(loc="best")
    fig2b_path = figures_dir / "data1_fig2_b_retentate_measured.png"
    fig2b.savefig(fig2b_path, bbox_inches="tight")
    plt.close(fig2b)

    fig2c = plt.figure(figsize=(4, 4))
    for i, v in enumerate(exp.vials, start=1):
        tf = (float(v.time_s[-1]) - t_delay) / 60.0
        yv = float(v.cV_avg) if isinstance(v.cV_avg, (int, float, np.number)) else np.nan
        if np.isfinite(yv):
            plt.plot([tf], [yv], "cs", markersize=6, label="Measured cV" if i == 1 else None)
    plt.xlabel("Time [min]")
    plt.ylabel("Vial concentration [mM]")
    plt.ylim(bottom=0)
    plt.tick_params(direction="in", top=True, right=True)
    plt.legend(loc="best")
    fig2c_path = figures_dir / "data1_fig2_c_vial_measured.png"
    fig2c.savefig(fig2c_path, bbox_inches="tight")
    plt.close(fig2c)

    fig2d = plt.figure(figsize=(4, 4))
    for i, v in enumerate(exp.vials, start=1):
        t = (np.asarray(v.time_s, dtype=float) - t_delay) / 60.0
        if v.mass_g is not None:
            plt.plot(t, np.asarray(v.mass_g, dtype=float), "r.", markersize=3, alpha=0.6, label="Mass meas." if i == 1 else None)
        if v.retentate_signal is not None:
            plt.plot(t, np.asarray(v.retentate_signal, dtype=float), "m-", linewidth=1.2, alpha=0.6, label="cF meas." if i == 1 else None)
        yv = float(v.cV_avg) if isinstance(v.cV_avg, (int, float, np.number)) else np.nan
        if np.isfinite(yv):
            plt.plot([(float(v.time_s[-1]) - t_delay) / 60.0], [yv], "c^", markersize=4, alpha=0.8, label="cV meas." if i == 1 else None)
    plt.xlabel("Time [min]")
    plt.ylabel("Measured responses")
    plt.tick_params(direction="in")
    plt.legend(loc="best")
    fig2d_path = figures_dir / "data1_fig2_d_combined_measured.png"
    fig2d.savefig(fig2d_path, bbox_inches="tight")
    plt.close(fig2d)

    fig3 = plt.figure(figsize=(4, 4))
    for n in range(1, n_vials + 1):
        v = exp.vials[n - 1]
        t_meas = (np.asarray(v.time_s, dtype=float) - t_delay) / 60.0
        t_sim = sim[n]["time_s"] / 60.0
        if v.mass_g is not None:
            y_mass = np.asarray(v.mass_g, dtype=float)
            plt.plot(t_meas, y_mass, "r.", markersize=3, alpha=0.6, label="Mass meas." if n == 1 else None)
            y_mass_pred = np.interp(t_meas * 60.0, sim[n]["time_s"], sim[n]["mV"])
            mass_true.extend(y_mass.tolist())
            mass_pred.extend(y_mass_pred.tolist())
        if v.retentate_signal is not None:
            y_cf = np.asarray(v.retentate_signal, dtype=float)
            plt.plot(t_meas, y_cf, "ms", markersize=3, alpha=0.6, label="cF meas." if n == 1 else None)
            plt.plot(t_sim, sim[n]["cF"], "g-", linewidth=1.4, alpha=0.8, label="cF pred." if n == 1 else None)
            y_cf_pred = np.interp(t_meas * 60.0, sim[n]["time_s"], sim[n]["cF"])
            cf_true.extend(y_cf.tolist())
            cf_pred.extend(y_cf_pred.tolist())
        y_cv = float(v.cV_avg) if isinstance(v.cV_avg, (int, float, np.number)) else np.nan
        if np.isfinite(y_cv):
            t_end = float(v.time_s[-1]) - t_delay
            y_cv_pred = float(np.interp(t_end, sim[n]["time_s"], sim[n]["cV"]))
            cv_true.append(y_cv)
            cv_pred.append(y_cv_pred)
            plt.plot([t_end / 60.0], [yv], "cs", markersize=5, alpha=0.8, label="cV meas." if n == 1 else None)
            plt.plot([sim[n]["time_s"][-1] / 60.0], [sim[n]["cV"][-1]], "r^", markersize=5, alpha=0.8, label="cV pred." if n == 1 else None)
    plt.xlabel("Time [min]")
    plt.ylabel("Concentration / mass responses")
    plt.tick_params(direction="in", top=True, right=True)
    plt.legend(loc="best")
    fig3_path = figures_dir / "data1_fig3_model_overlay.png"
    fig3.savefig(fig3_path, bbox_inches="tight")
    plt.close(fig3)

    objective = float(est.get("objective")) if est.get("objective") is not None else np.nan
    cov_method = est.get("covariance_method")
    cov_warn = est.get("covariance_warning")

    t1_metrics = pd.DataFrame(
        [
            {"dataset_id": "DATA1_511.12", "metric": "objective", "unified_value": objective, "target_id": "D1-M-T1"},
            {"dataset_id": "DATA1_511.12", "metric": "theta.Lp", "unified_value": lp, "target_id": "D1-M-T1"},
            {"dataset_id": "DATA1_511.12", "metric": "theta.B", "unified_value": b, "target_id": "D1-M-T1"},
            {"dataset_id": "DATA1_511.12", "metric": "theta.sigma", "unified_value": sigma, "target_id": "D1-M-T1"},
            {"dataset_id": "DATA1_511.12", "metric": "covariance.method", "unified_value": np.nan, "target_id": "D1-M-T1"},
        ]
    )
    t1_metrics.loc[t1_metrics["metric"] == "covariance.method", "unified_text"] = str(cov_method) if cov_method is not None else str(cov_warn)

    rmse_mass = _rmse(mass_true, mass_pred)
    rmse_cf = _rmse(cf_true, cf_pred)
    rmse_cv = _rmse(cv_true, cv_pred)
    nrmse_mass = nrmse_vs_range(np.asarray(mass_true, dtype=float), np.asarray(mass_pred, dtype=float))
    nrmse_cf = nrmse_vs_range(np.asarray(cf_true, dtype=float), np.asarray(cf_pred, dtype=float))
    nrmse_cv = nrmse_vs_range(np.asarray(cv_true, dtype=float), np.asarray(cv_pred, dtype=float))

    t2_metrics = pd.DataFrame(
        [
            {"dataset_id": "DATA1_511.12", "metric": "rmse.mass", "unified_value": rmse_mass, "target_id": "D1-M-T2"},
            {"dataset_id": "DATA1_511.12", "metric": "rmse.cF", "unified_value": rmse_cf, "target_id": "D1-M-T2"},
            {"dataset_id": "DATA1_511.12", "metric": "rmse.cV", "unified_value": rmse_cv, "target_id": "D1-M-T2"},
            {"dataset_id": "DATA1_511.12", "metric": "nrmse.mass", "unified_value": nrmse_mass, "target_id": "D1-M-T2"},
            {"dataset_id": "DATA1_511.12", "metric": "nrmse.cF", "unified_value": nrmse_cf, "target_id": "D1-M-T2"},
            {"dataset_id": "DATA1_511.12", "metric": "nrmse.cV", "unified_value": nrmse_cv, "target_id": "D1-M-T2"},
            {"dataset_id": "DATA1_511.12", "metric": "nrmse.pass_threshold", "unified_value": tolerances.curve_nrmse_range, "target_id": "D1-M-T2"},
        ]
    )
    unified_stage_metrics = pd.concat([t1_metrics, t2_metrics], ignore_index=True)
    unified_stage_metrics.to_csv(tables_dir / "data1_stageA_unified_metrics.csv", index=False)

    refs = references.copy()
    refs_t1 = refs[(refs["dataset_id"] == "DATA1_511.12") & (refs["target_id"] == "D1-M-T1")].copy()
    refs_t2 = refs[(refs["dataset_id"] == "DATA1_511.12") & (refs["target_id"] == "D1-M-T2")].copy()
    t1_side = refs_t1.merge(t1_metrics[["dataset_id", "metric", "unified_value"]], on=["dataset_id", "metric"], how="outer")
    t2_side = refs_t2.merge(t2_metrics[["dataset_id", "metric", "unified_value"]], on=["dataset_id", "metric"], how="outer")

    for df in (t1_side, t2_side):
        def _rel(row):
            pv = row.get("paper_value")
            uv = row.get("unified_value")
            if pd.isna(pv) or pd.isna(uv):
                return np.nan
            return rel_err(float(uv), float(pv))
        df["relative_error"] = df.apply(_rel, axis=1)

    t1_side["status"] = t1_side["relative_error"].map(
        lambda e: "MISSING_VALUE" if pd.isna(e) else ("PASS" if e <= tolerances.objective_rel else "FAIL")
    )
    t2_side["status"] = t2_side["relative_error"].map(
        lambda e: "MISSING_VALUE" if pd.isna(e) else ("PASS" if e <= tolerances.param_rel_primary else ("PASS_WITH_EXPLANATION" if e <= tolerances.param_rel_with_explanation else "FAIL"))
    )
    t1_path = tables_dir / "data1_table1_side_by_side.csv"
    t2_path = tables_dir / "data1_table2_side_by_side.csv"
    t1_side.to_csv(t1_path, index=False)
    t2_side.to_csv(t2_path, index=False)

    return {
        "figures": [str(fig2a_path), str(fig2b_path), str(fig2c_path), str(fig2d_path), str(fig3_path)],
        "tables": [str(t1_path), str(t2_path), str(tables_dir / "data1_stageA_unified_metrics.csv")],
        "summary": {
            "objective": objective,
            "theta": {"Lp": lp, "B": b, "sigma": sigma},
            "rmse": {"mass": rmse_mass, "cF": rmse_cf, "cV": rmse_cv},
            "nrmse": {"mass": nrmse_mass, "cF": nrmse_cf, "cV": nrmse_cv},
        },
    }


def build_data1_stage_b_outputs(*, out_dir: Path, exp, est: Dict[str, object], sim: Dict[int, Dict[str, np.ndarray]]) -> Dict[str, object]:
    """Generate DATA1 Stage B artifact bundle."""
    import matplotlib.pyplot as plt

    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    t_delay = float(exp.vials[0].time_s[0])
    n_vials = len(exp.vials)
    produced: List[tuple[str, str]] = []

    fig4a = plt.figure(figsize=(5, 4))
    for n in range(1, n_vials + 1):
        plt.plot(sim[n]["time_s"] / 60.0, sim[n]["Jw"], linewidth=1.2, alpha=0.7)
    plt.xlabel("Time [min]")
    plt.ylabel("Jw [cm/s]")
    plt.tick_params(direction="in", top=True, right=True)
    p = figures_dir / "data1_fig4_1_jw_trajectories.png"
    fig4a.savefig(p, bbox_inches="tight")
    plt.close(fig4a)
    produced.append(("D1-M-F4", str(p)))

    fig4b = plt.figure(figsize=(5, 4))
    for n in range(1, n_vials + 1):
        plt.plot(sim[n]["time_s"] / 60.0, sim[n]["Js"], linewidth=1.2, alpha=0.7)
    plt.xlabel("Time [min]")
    plt.ylabel("Js [umol/cm2/s]")
    plt.tick_params(direction="in", top=True, right=True)
    p = figures_dir / "data1_fig4_2_js_trajectories.png"
    fig4b.savefig(p, bbox_inches="tight")
    plt.close(fig4b)
    produced.append(("D1-M-F4", str(p)))

    mass_true: List[float] = []
    mass_pred: List[float] = []
    cf_true: List[float] = []
    cf_pred: List[float] = []
    cv_true: List[float] = []
    cv_pred: List[float] = []
    for n, v in enumerate(exp.vials, start=1):
        t_meas = np.asarray(v.time_s, dtype=float) - t_delay
        if v.mass_g is not None:
            y = np.asarray(v.mass_g, dtype=float)
            yp = np.interp(t_meas, sim[n]["time_s"], sim[n]["mV"])
            mass_true.extend(y.tolist())
            mass_pred.extend(yp.tolist())
        if v.retentate_signal is not None:
            y = np.asarray(v.retentate_signal, dtype=float)
            yp = np.interp(t_meas, sim[n]["time_s"], sim[n]["cF"])
            cf_true.extend(y.tolist())
            cf_pred.extend(yp.tolist())
        yv = float(v.cV_avg) if isinstance(v.cV_avg, (int, float, np.number)) else np.nan
        if np.isfinite(yv):
            yp = float(np.interp(float(v.time_s[-1] - t_delay), sim[n]["time_s"], sim[n]["cV"]))
            cv_true.append(yv)
            cv_pred.append(yp)

    def _parity_plot(y_t: List[float], y_p: List[float], label: str, path: Path) -> None:
        yt = np.asarray(y_t, dtype=float)
        yp = np.asarray(y_p, dtype=float)
        mask = np.isfinite(yt) & np.isfinite(yp)
        yt = yt[mask]
        yp = yp[mask]
        fig = plt.figure(figsize=(4, 4))
        if yt.size:
            plt.plot(yt, yp, "o", markersize=3, alpha=0.75)
            lo = float(min(np.min(yt), np.min(yp)))
            hi = float(max(np.max(yt), np.max(yp)))
            plt.plot([lo, hi], [lo, hi], "k--", linewidth=1.0)
        plt.xlabel(f"{label} measured")
        plt.ylabel(f"{label} predicted")
        plt.tick_params(direction="in", top=True, right=True)
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)

    p = figures_dir / "data1_fig5_1_mass_parity.png"
    _parity_plot(mass_true, mass_pred, "mass", p)
    produced.append(("D1-M-F5", str(p)))
    p = figures_dir / "data1_fig5_2_cF_parity.png"
    _parity_plot(cf_true, cf_pred, "cF", p)
    produced.append(("D1-M-F5", str(p)))
    p = figures_dir / "data1_fig5_3_cV_parity.png"
    _parity_plot(cv_true, cv_pred, "cV", p)
    produced.append(("D1-M-F5", str(p)))

    fig6 = plt.figure(figsize=(5, 4))
    for n, v in enumerate(exp.vials, start=1):
        if v.retentate_signal is None:
            continue
        t_meas = np.asarray(v.time_s, dtype=float) - t_delay
        y = np.asarray(v.retentate_signal, dtype=float)
        yp = np.interp(t_meas, sim[n]["time_s"], sim[n]["cF"])
        plt.plot(t_meas / 60.0, yp - y, "-", linewidth=1.0, alpha=0.7)
    plt.axhline(0.0, color="k", linestyle="--", linewidth=1.0)
    plt.xlabel("Time [min]")
    plt.ylabel("cF residual (pred - meas)")
    plt.tick_params(direction="in", top=True, right=True)
    p = figures_dir / "data1_fig6_1_cf_residuals.png"
    fig6.savefig(p, bbox_inches="tight")
    plt.close(fig6)
    produced.append(("D1-M-F6", str(p)))

    si_targets = ["D1-S-FS2", "D1-S-FS3", "D1-S-FS4", "D1-S-FS5", "D1-S-FS6", "D1-S-FS7"]
    group_size = max(1, int(np.ceil(n_vials / len(si_targets))))
    for i, target in enumerate(si_targets):
        start = i * group_size + 1
        end = min(n_vials, (i + 1) * group_size)
        fig = plt.figure(figsize=(5, 4))
        for n in range(start, end + 1):
            if n < 1 or n > n_vials:
                continue
            v = exp.vials[n - 1]
            t = (np.asarray(v.time_s, dtype=float) - t_delay) / 60.0
            if v.retentate_signal is not None:
                plt.plot(t, np.asarray(v.retentate_signal, dtype=float), ".", markersize=2, alpha=0.5)
                plt.plot(sim[n]["time_s"] / 60.0, sim[n]["cF"], "-", linewidth=1.0, alpha=0.7)
        plt.xlabel("Time [min]")
        plt.ylabel("Retentate concentration [mM]")
        plt.tick_params(direction="in", top=True, right=True)
        p = figures_dir / f"data1_si_figS{i+2}_1_overlay.png"
        fig.savefig(p, bbox_inches="tight")
        plt.close(fig)
        produced.append((target, str(p)))

    manifest = pd.DataFrame(produced, columns=["target_id", "artifact_path"])
    manifest_path = tables_dir / "data1_stageB_artifact_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    return {"manifest": str(manifest_path), "counts_by_target": manifest["target_id"].value_counts().to_dict()}


def build_data2_validation_artifacts(
    *,
    out_dir: Path,
    dataset_id: str,
    exp,
    curve_metrics: Dict[str, float],
    simulation_validation_root: Path,
) -> Dict[str, object]:
    """Generate DATA2 validation artifacts that are currently exercised in regression/nightly paths."""
    import matplotlib.pyplot as plt

    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    artifacts: Dict[str, object] = {"tables": [], "figures": []}

    if len(exp.vials) >= 4 and exp.vials[3].mass_g is not None:
        vial = exp.vials[3]
        measurement_trace = pd.DataFrame(
            {
                "time_min": np.asarray(vial.time_s, dtype=float) / 60.0,
                "mass_g": np.asarray(vial.mass_g, dtype=float),
            }
        ).sort_values("time_min")
        trace_csv = tables_dir / f"{dataset_id.lower()}_fig3_measurement_trace.csv"
        measurement_trace.to_csv(trace_csv, index=False)
        artifacts["tables"].append(str(trace_csv))

        fig = plt.figure(figsize=(4.5, 3.5))
        plt.plot(measurement_trace["time_min"], measurement_trace["mass_g"], "k.-", linewidth=1.2, markersize=3)
        plt.xlabel("Time [min]")
        plt.ylabel("Mass [g]")
        plt.tick_params(direction="in", top=True, right=True)
        fig_path = figures_dir / f"{dataset_id.lower()}_fig3_measurement_trace.png"
        fig.savefig(fig_path, bbox_inches="tight")
        plt.close(fig)
        artifacts["figures"].append(str(fig_path))

        baseline_path = simulation_validation_root / "data2_main" / "fig3" / "data2_main_fig3_mass_tc.csv"
        if dataset_id == "DATA2_270611.123" and baseline_path.exists():
            baseline = pd.read_csv(baseline_path)
            baseline_measure = baseline[baseline["series_id"] == "measurements"].copy().sort_values("time_min")
            common = baseline_measure.merge(measurement_trace, on="time_min", suffixes=("_baseline", "_unified"))
            err = common["mass_g_unified"].to_numpy(dtype=float) - common["mass_g_baseline"].to_numpy(dtype=float)
            comparison = pd.DataFrame(
                [
                    {
                        "target": "DATA2_Fig3_measurements",
                        "dataset_id": dataset_id,
                        "n_points": int(len(common)),
                        "mae": float(np.mean(np.abs(err))) if len(err) else np.nan,
                        "max_abs_err": float(np.max(np.abs(err))) if len(err) else np.nan,
                        "rmse": float(np.sqrt(np.mean(err**2))) if len(err) else np.nan,
                    }
                ]
            )
            comparison_csv = tables_dir / f"{dataset_id.lower()}_fig3_measurement_comparison.csv"
            comparison.to_csv(comparison_csv, index=False)
            artifacts["tables"].append(str(comparison_csv))

    curve_rows = []
    for metric_name in ("nrmse.mass", "nrmse.cF", "nrmse.cV"):
        curve_rows.append(
            {
                "dataset_id": dataset_id,
                "metric": metric_name,
                "unified_value": curve_metrics.get(metric_name, np.nan),
                "target_id": "D2-M-F6",
            }
        )
    curve_df = pd.DataFrame(curve_rows)
    curve_csv = tables_dir / f"{dataset_id.lower()}_fig6_curve_metrics.csv"
    curve_df.to_csv(curve_csv, index=False)
    artifacts["tables"].append(str(curve_csv))
    return artifacts


def flatten_theta(theta_obj: object) -> Dict[str, float]:
    """Normalize theta output into a simple float dict."""
    if isinstance(theta_obj, pd.Series):
        return {str(k): float(v) for k, v in theta_obj.items()}
    if isinstance(theta_obj, dict):
        out: Dict[str, float] = {}
        for k, v in theta_obj.items():
            try:
                out[str(k)] = float(v)
            except Exception:
                continue
        return out
    return {}


def rel_err(unified: float, reference: float) -> float:
    if reference == 0:
        return abs(unified - reference)
    return abs(unified - reference) / abs(reference)


def nrmse_vs_range(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if y_true.size == 0 or y_true.size != y_pred.size:
        return np.inf
    rng = float(np.nanmax(y_true) - np.nanmin(y_true))
    if rng == 0:
        return 0.0
    rmse = float(np.sqrt(np.nanmean((y_true - y_pred) ** 2)))
    return rmse / rng


def _rmse(y_true: List[float], y_pred: List[float]) -> float:
    yt = np.asarray(y_true, dtype=float)
    yp = np.asarray(y_pred, dtype=float)
    if yt.size == 0 or yt.size != yp.size:
        return np.nan
    return float(np.sqrt(np.mean((yt - yp) ** 2)))
