#!/usr/bin/env python3
"""Deterministic DATA1 figure regeneration from DiafiltrationPaperPlots logic."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import datetime
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use("Agg")


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA1_ROOT = REPO_ROOT / "DATA1_matlab"
DATA_DIR = DATA1_ROOT / "data"


def _prepare_output_dir(run_id: str, output_dir: str | None) -> Path:
    if output_dir:
        out = Path(output_dir).expanduser().resolve()
    else:
        out = REPO_ROOT / "results" / "reproduction" / run_id / "figures" / "data1_notebook_regen"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _copy_generated_pngs(tmp_dir: Path, out_dir: Path) -> list[str]:
    generated: list[str] = []
    for png in sorted(tmp_dir.glob("*.png")):
        dst = out_dir / png.name
        shutil.copy2(png, dst)
        generated.append(png.name)
    return generated


def _run_notebook_logic(tmp_dir: Path) -> None:
    import sys

    sys.path.insert(0, str(DATA1_ROOT))
    from diafiltration_plots import (  # pylint: disable=import-error
        calib_curve_cond,
        loadmat,
        plot_conc_range,
        plot_contour,
        plot_contour_sig_sen,
        plot_cr_measure,
        plot_sim,
        plot_sim_comparison,
        plot_sim_show,
    )

    def dpath(*parts: str) -> Path:
        return DATA_DIR.joinpath(*parts)

    def read_contour_csv(rel_path: str) -> pd.DataFrame:
        return pd.read_csv(dpath(*rel_path.split("/")))

    cwd_before = Path.cwd()
    try:
        # Functions save PNGs in current working directory.
        import os

        os.chdir(tmp_dir)

        # Filtration: M1/M2/M3/M4
        data_stru_f = loadmat(str(dpath("data_stru-dataset501.1.mat")))["data_stru"]
        fit_stru = loadmat(str(dpath("501.1 concpolar", "fit_stru.mat")))["fit_stru"]
        plot_sim_comparison(data_stru_f, fit_stru, plot_pred=True, cond=True, lg=False)
        plot_contour(read_contour_csv("501.1 concpolar/contourdata-x_B-y_Lp.csv"))
        plot_contour(read_contour_csv("501.1 concpolar/contourdata-x_sigma-y_Lp.csv"))

        data_stru_f = loadmat(str(dpath("data_stru-dataset501.11.mat")))["data_stru"]
        fit_stru = loadmat(str(dpath("501.11 concpolar", "fit_stru.mat")))["fit_stru"]
        plot_sim_comparison(data_stru_f, fit_stru, plot_pred=True, cond=True, lg=False)
        plot_contour(read_contour_csv("501.11 concpolar/contourdata-x_B-y_Lp.csv"))
        plot_contour(read_contour_csv("501.11 concpolar/contourdata-x_sigma-y_Lp.csv"))

        data_stru_f = loadmat(str(dpath("data_stru-dataset501.1.mat")))["data_stru"]
        fit_stru = loadmat(str(dpath("501.1", "fit_stru.mat")))["fit_stru"]
        plot_sim_comparison(data_stru_f, fit_stru, plot_pred=True, cond=True, lg=False)
        plot_contour(read_contour_csv("501.1/contourdata-x_B-y_Lp.csv"))
        plot_contour(read_contour_csv("501.1/contourdata-x_sigma-y_Lp.csv"))

        data_stru_f = loadmat(str(dpath("data_stru-dataset501.11.mat")))["data_stru"]
        fit_stru = loadmat(str(dpath("501.11", "fit_stru.mat")))["fit_stru"]
        plot_sim_comparison(data_stru_f, fit_stru, plot_pred=True, cond=False, lg=False)
        plot_contour(read_contour_csv("501.11/contourdata-x_B-y_Lp.csv"))
        plot_contour(read_contour_csv("501.11/contourdata-x_sigma-y_Lp.csv"))

        # Filtration sigma sensitivity.
        sigma = [0.1, 0.5, 0.9]
        colorstring = "rbg"
        for sig in sigma:
            sim_stru = loadmat(
                str(dpath("sigma sensitivity", f"sim_stru-dat501.1 C_Fin5.2843sig{sig:.1f}.mat"))
            )["sim_stru"]
            line = "dashed" if sig == 0.1 else "solid" if sig == 0.5 else "dotted"
            color = colorstring[sigma.index(sig)]
            plot_sim(sim_stru, color, line)
        plot_sim_show(sigma, colorstring)

        # Sigma sensitivity contour maps.
        data_stru_f = loadmat(str(dpath("data_stru-dataset501.1.mat")))["data_stru"]
        contour_sig_stru = loadmat(str(dpath("501.1 concpolar", "contour_sig_stru-dat501.1.mat")))[
            "contour_sig_stru"
        ]
        plot_contour_sig_sen(
            data_stru_f,
            contour_sig_stru,
            Vmin=0,
            Vmax=3,
            Level=[[1, 2, 3], [1, 2, 3], [1, 1.3, 1.6]],
            Colorbar_ticks=[0, 1, 2, 3],
            Manual_locations=[[(20, 40), (35, 60), (50, 80)], [(10, 20), (15, 30), (20, 40)], [(50, 70), (60, 60), (75, 57)]],
            name_append="",
            filled=True,
            bar=True,
        )

        data_stru_d = loadmat(str(dpath("data_stru-dataset511.12.mat")))["data_stru"]
        for endvial, level, ml in [
            (
                "endvial1",
                [[1, 1.5, 1.9], [1, 1.5, 2], [0.1, 0.2, 0.3]],
                [[(40, 60), (60, 100), (73.5, 110)], [(30, 40), (50, 50), (70, 60)], [(40, 100), (40, 50), (40, 30)]],
            ),
            (
                "endvial5",
                [[2, 5, 8], [2, 5, 8], [0.3, 0.4, 0.5]],
                [[(20, 40), (35, 60), (50, 80)], [(15, 40), (40, 50), (70, 60)], [(20, 30), (50, 20), (15, 20), (60, 80), (2, 5), (75, 90)]],
            ),
            (
                "endvial10",
                [[2, 5, 8], [2, 5, 8], [0.5, 0.6, 0.7]],
                [[(5, 10), (15, 30), (20, 45)], [(5, 50), (15, 60), (20, 70)], [(33, 40), (40, 58), (50, 80)]],
            ),
        ]:
            contour_sig_stru = loadmat(
                str(dpath("511.12 concpolar", f"contour_sig_stru-dat511.12-{endvial}.mat"))
            )["contour_sig_stru"]
            plot_contour_sig_sen(
                data_stru_d,
                contour_sig_stru,
                Vmin=0,
                Vmax=3,
                Level=level,
                Colorbar_ticks=[0, 1, 2, 3],
                Manual_locations=ml,
                name_append=f"-{endvial}",
                filled=True,
            )

        # Diafiltration: M1/M2/M3/M4
        data_stru_d = loadmat(str(dpath("data_stru-dataset511.12.mat")))["data_stru"]
        fit_stru = loadmat(str(dpath("511.12 concpolar", "fit_stru.mat")))["fit_stru"]
        plot_sim_comparison(data_stru_d, fit_stru, plot_pred=True, cond=True, lg=False)
        plot_contour(read_contour_csv("511.12 concpolar/contourdata-x_B-y_Lp.csv"))
        plot_contour(read_contour_csv("511.12 concpolar/contourdata-x_sigma-y_Lp.csv"))
        plot_sim_comparison(data_stru_d, fit_stru, plot_pred=True, cond=True, lg=False, preface=True)
        plot_contour(read_contour_csv("511.12 concpolar/contourdata-x_B-y_Lp.csv"), show_title=False, preface=True)
        plot_contour(read_contour_csv("511.12 concpolar/contourdata-x_sigma-y_Lp.csv"), show_title=False, preface=True)

        data_stru_d = loadmat(str(dpath("data_stru-dataset511.11.mat")))["data_stru"]
        fit_stru = loadmat(str(dpath("511.11 concpolar", "fit_stru.mat")))["fit_stru"]
        plot_sim_comparison(data_stru_d, fit_stru, plot_pred=True, cond=True, lg=False)
        plot_contour(read_contour_csv("511.11 concpolar/contourdata-x_B-y_Lp.csv"))
        plot_contour(read_contour_csv("511.11 concpolar/contourdata-x_sigma-y_Lp.csv"))

        data_stru_d = loadmat(str(dpath("data_stru-dataset511.12.mat")))["data_stru"]
        fit_stru = loadmat(str(dpath("511.12", "fit_stru.mat")))["fit_stru"]
        plot_sim_comparison(data_stru_d, fit_stru, plot_pred=True, cond=True, lg=False)
        plot_contour(read_contour_csv("511.12/contourdata-x_B-y_Lp.csv"))
        plot_contour(read_contour_csv("511.12/contourdata-x_sigma-y_Lp.csv"))

        data_stru_d = loadmat(str(dpath("data_stru-dataset511.11.mat")))["data_stru"]
        fit_stru = loadmat(str(dpath("511.11", "fit_stru.mat")))["fit_stru"]
        plot_sim_comparison(data_stru_d, fit_stru, plot_pred=True, cond=False, lg=False)
        plot_contour(read_contour_csv("511.11/contourdata-x_B-y_Lp.csv"))
        plot_contour(read_contour_csv("511.11/contourdata-x_sigma-y_Lp.csv"))

        # Diafiltration sigma sensitivity.
        sigma = [0.1, 0.5, 0.9]
        colorstring = "rbg"
        for sig in sigma:
            sim_stru = loadmat(
                str(dpath("sigma sensitivity", f"sim_stru-dat511.12 C_Fin15.2052sig{sig:.1f}.mat"))
            )["sim_stru"]
            line = "dashed" if sig == 0.1 else "solid" if sig == 0.5 else "dotted"
            color = colorstring[sigma.index(sig)]
            plot_sim(sim_stru, color, line)
        plot_sim_show(sigma, colorstring)

        # Experiment space + concentration-ratio measure + calibration.
        df_f = pd.read_csv(dpath("experiment space", "filtration.csv"), header=2)
        df_d = pd.read_csv(dpath("experiment space", "diafiltration.csv"), header=2)
        plot_conc_range(df_f, df_d)

        cr_pred = df_d.D3_cf
        data_stru_d = loadmat(str(dpath("data_stru-dataset511.12.mat")))["data_stru"]
        fit_stru = loadmat(str(dpath("511.12 concpolar", "fit_stru.mat")))["fit_stru"]
        plot_cr_measure(data_stru_d, fit_stru, cr_pred, ybottom=0, lg=True)

        calib_curve_data = pd.read_csv(
            dpath("experiment space", "conductivity_calibration.csv"), header=0, skiprows=[1]
        )
        calib_curve_cond(calib_curve_data)

        df = pd.read_csv(dpath("experiment space", "Classical_analysis-dat301.1.csv"), header=0)
        cr_pred = df.cf
        data_stru = loadmat(
            str(dpath("dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0", "data_stru-dataset301.1.mat"))
        )["data_stru"]
        fit_stru = loadmat(
            str(
                dpath(
                    "dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0",
                    "fit_stru-dat 301.1 oneCPNT holdup concpolar cvmv fixed ch0.mat",
                )
            )
        )["fit_stru"]
        plot_cr_measure(data_stru, fit_stru, cr_pred, ybottom=4, cond=False)

        df = pd.read_csv(dpath("experiment space", "Classical_analysis-dat501.1.csv"), header=0)
        cr_pred = df.cf
        data_stru_f = loadmat(str(dpath("data_stru-dataset501.1.mat")))["data_stru"]
        fit_stru = loadmat(str(dpath("501.1 concpolar", "fit_stru.mat")))["fit_stru"]
        plot_cr_measure(data_stru_f, fit_stru, cr_pred, ybottom=4)
    finally:
        plt.close("all")
        os.chdir(cwd_before)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", default=datetime.now().strftime("%Y%m%d-%H%M%S"))
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()

    out_dir = _prepare_output_dir(run_id=args.run_id, output_dir=args.output_dir)
    tmp_dir = out_dir / "_tmp_workdir"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    _run_notebook_logic(tmp_dir)
    generated = _copy_generated_pngs(tmp_dir, out_dir)
    shutil.rmtree(tmp_dir, ignore_errors=True)

    manifest = {
        "run_id": args.run_id,
        "output_dir": str(out_dir),
        "generated_png_count": len(generated),
        "generated_pngs": generated,
    }
    (out_dir / "manifest_data1_notebook_regen.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
