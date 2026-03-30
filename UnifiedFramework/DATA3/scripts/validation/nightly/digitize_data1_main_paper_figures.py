#!/usr/bin/env python3
"""Generate first-pass CSV digitizations for DATA1 main paper figures.

This script extracts colored raster traces from the embedded paper figures
already stored under the paper PDF extraction run. It intentionally produces
data-coordinate point clouds rather than claiming exact publication-quality
manual digitization.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import cv2
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[5]
DATA1_MAIN_IMAGE_ROOT = (
    REPO_ROOT
    / "UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/data1_main/images"
)
OUT_ROOT = (
    REPO_ROOT
    / "UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/paper/data1_main"
)


@dataclass(frozen=True)
class PanelSpec:
    figure_id: str
    panel_id: str
    image_name: str
    crop_box: tuple[int, int, int, int]
    x_min: float
    x_max: float
    y_min: float
    y_max: float
    notes: str = ""


PANELS: tuple[PanelSpec, ...] = (
    PanelSpec("fig2", "A", "img-004.jpg", (160, 56, 573, 615), 0.0, 40.0, 0.0, 0.8, "Mass vs time."),
    PanelSpec("fig2", "B", "img-004.jpg", (910, 56, 1490, 615), 0.0, 40.0, 0.0, 8.2, "Concentration vs time."),
    PanelSpec("fig2", "C", "img-004.jpg", (160, 775, 733, 1322), 0.0, 260.0, 0.0, 1.25, "Mass vs time."),
    PanelSpec("fig2", "D", "img-004.jpg", (910, 775, 1492, 1322), 0.0, 260.0, 0.0, 85.0, "Concentration vs time."),
    PanelSpec("fig3", "A", "img-005.jpg", (130, 23, 826, 786), 0.0, 90.0, 0.0, 26.0, "Permeate vs retentate scatter."),
    PanelSpec("fig4", "A_mass", "img-006.jpg", (168, 55, 589, 382), 0.0, 40.0, 0.0, 0.85, "Sigma sensitivity mass."),
    PanelSpec("fig4", "A_permeate", "img-006.jpg", (710, 55, 1087, 382), 0.0, 38.0, 0.62, 1.14, "Sigma sensitivity permeate."),
    PanelSpec("fig4", "A_retentate", "img-006.jpg", (1216, 55, 1570, 382), 0.0, 38.0, 5.2, 8.35, "Sigma sensitivity retentate."),
    PanelSpec("fig4", "B_mass", "img-006.jpg", (169, 526, 588, 923), 0.0, 260.0, 0.0, 4.0, "Sigma sensitivity mass, long horizon."),
    PanelSpec("fig4", "B_permeate", "img-006.jpg", (710, 526, 1087, 923), 0.0, 260.0, 0.8, 19.0, "Sigma sensitivity permeate, long horizon."),
    PanelSpec("fig4", "B_retentate", "img-006.jpg", (1216, 526, 1570, 923), 0.0, 260.0, 15.0, 165.0, "Sigma sensitivity retentate, long horizon."),
    PanelSpec("fig5", "A_mass", "img-007.jpg", (100, 61, 564, 504), 0.0, 1.0, 0.5, 7.4, "Contour plot over sigma and Lp."),
    PanelSpec("fig5", "A_permeate", "img-007.jpg", (637, 61, 1095, 504), 0.0, 1.0, 0.5, 7.4, "Contour plot over sigma and Lp."),
    PanelSpec("fig5", "A_retentate", "img-007.jpg", (1170, 61, 1560, 504), 0.0, 1.0, 0.5, 7.4, "Contour plot over sigma and Lp."),
    PanelSpec("fig5", "B_mass", "img-007.jpg", (100, 581, 564, 1024), 0.0, 1.0, 0.5, 7.4, "Contour plot over sigma and Lp."),
    PanelSpec("fig5", "B_permeate", "img-007.jpg", (637, 581, 1095, 1024), 0.0, 1.0, 0.5, 7.4, "Contour plot over sigma and Lp."),
    PanelSpec("fig5", "B_retentate", "img-007.jpg", (1170, 581, 1560, 1024), 0.0, 1.0, 0.5, 7.4, "Contour plot over sigma and Lp."),
    PanelSpec("fig5", "C_mass", "img-007.jpg", (100, 1104, 564, 1545), 0.0, 1.0, 0.5, 7.4, "Contour plot over sigma and Lp."),
    PanelSpec("fig5", "C_permeate", "img-007.jpg", (637, 1104, 1095, 1545), 0.0, 1.0, 0.5, 7.4, "Contour plot over sigma and Lp."),
    PanelSpec("fig5", "C_retentate", "img-007.jpg", (1170, 1104, 1560, 1545), 0.0, 1.0, 0.5, 7.4, "Contour plot over sigma and Lp."),
    PanelSpec("fig6", "A_mass", "img-008.jpg", (100, 61, 570, 503), 0.0, 2.0, 0.5, 7.4, "Contour plot over B and Lp."),
    PanelSpec("fig6", "A_permeate", "img-008.jpg", (636, 61, 1096, 503), 0.0, 2.0, 0.5, 7.4, "Contour plot over B and Lp."),
    PanelSpec("fig6", "A_retentate", "img-008.jpg", (1170, 61, 1560, 503), 0.0, 2.0, 0.5, 7.4, "Contour plot over B and Lp."),
    PanelSpec("fig6", "B_mass", "img-008.jpg", (100, 581, 570, 1022), 0.0, 2.0, 0.5, 7.4, "Contour plot over B and Lp."),
    PanelSpec("fig6", "B_permeate", "img-008.jpg", (636, 581, 1096, 1022), 0.0, 2.0, 0.5, 7.4, "Contour plot over B and Lp."),
    PanelSpec("fig6", "B_retentate", "img-008.jpg", (1170, 581, 1560, 1022), 0.0, 2.0, 0.5, 7.4, "Contour plot over B and Lp."),
    PanelSpec("fig6", "C_mass", "img-008.jpg", (100, 1104, 570, 1545), 0.0, 2.0, 0.5, 7.4, "Contour plot over B and Lp."),
    PanelSpec("fig6", "C_permeate", "img-008.jpg", (636, 1104, 1096, 1545), 0.0, 2.0, 0.5, 7.4, "Contour plot over B and Lp."),
    PanelSpec("fig6", "C_retentate", "img-008.jpg", (1170, 1104, 1560, 1545), 0.0, 2.0, 0.5, 7.4, "Contour plot over B and Lp."),
)


COLOR_RANGES = {
    "red": [((0, 70, 70), (12, 255, 255)), ((170, 70, 70), (180, 255, 255))],
    "orange": [((10, 70, 70), (28, 255, 255))],
    "yellow": [((28, 60, 90), (40, 255, 255))],
    "green": [((40, 40, 50), (90, 255, 255))],
    "cyan": [((80, 40, 50), (105, 255, 255))],
    "blue": [((100, 40, 50), (136, 255, 255))],
    "magenta": [((136, 40, 50), (170, 255, 255))],
    "gray": [((0, 0, 60), (180, 45, 230))],
}


def _resolve_points(mask: np.ndarray, min_area: int = 6) -> np.ndarray:
    mask_u8 = (mask.astype(np.uint8)) * 255
    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask_u8, connectivity=8)
    keep = np.zeros_like(mask_u8)
    for label in range(1, num_labels):
        area = int(stats[label, cv2.CC_STAT_AREA])
        if area >= min_area:
            keep[labels == label] = 255
    ys, xs = np.where(keep > 0)
    if xs.size == 0:
        return np.empty((0, 2), dtype=float)
    return np.column_stack([xs, ys]).astype(float)


def _mask_for_family(hsv: np.ndarray, family: str) -> np.ndarray:
    mask = np.zeros(hsv.shape[:2], dtype=bool)
    for lower, upper in COLOR_RANGES[family]:
        current = cv2.inRange(hsv, np.array(lower, dtype=np.uint8), np.array(upper, dtype=np.uint8)) > 0
        mask |= current
    return mask


def _pixel_to_data(points: np.ndarray, panel: PanelSpec) -> np.ndarray:
    x0, y0, x1, y1 = panel.crop_box
    width = max(x1 - x0 - 1, 1)
    height = max(y1 - y0 - 1, 1)
    x_data = panel.x_min + (points[:, 0] / width) * (panel.x_max - panel.x_min)
    y_data = panel.y_max - (points[:, 1] / height) * (panel.y_max - panel.y_min)
    return np.column_stack([x_data, y_data])


def _write_points_csv(panel: PanelSpec, rows: list[dict[str, object]]) -> Path:
    out_dir = OUT_ROOT / panel.figure_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"data1_main_{panel.figure_id}_{panel.panel_id.lower()}_points.csv"
    fieldnames = [
        "figure_id",
        "panel_id",
        "series_family",
        "x_data",
        "y_data",
        "x_pixel",
        "y_pixel",
        "source_image",
        "notes",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return out_path


def _write_metadata(rows: list[dict[str, object]]) -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    out_path = OUT_ROOT / "data1_main_panel_metadata.csv"
    fieldnames = [
        "figure_id",
        "panel_id",
        "source_image",
        "crop_x0",
        "crop_y0",
        "crop_x1",
        "crop_y1",
        "x_min",
        "x_max",
        "y_min",
        "y_max",
        "notes",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_readme() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    readme_path = OUT_ROOT / "README.md"
    readme_path.write_text(
        "# DATA1 Main Paper Digitized CSVs\n\n"
        "These CSVs are first-pass raster digitizations generated from the extracted\n"
        "embedded figures in the DATA1 main paper. They are intended as reproducible\n"
        "paper-side baselines for follow-up cleanup, not as a claim of perfect manual\n"
        "WebPlotDigitizer-quality extraction.\n\n"
        "Conventions:\n"
        "- one CSV per panel,\n"
        "- `series_family` groups colored traces/markers,\n"
        "- `x_data` and `y_data` are derived from manual panel calibration,\n"
        "- contour figures export trace clouds by color family rather than labeled\n"
        "  contour values,\n"
        "- red triangle reference markers are included under the `red` family where present.\n",
        encoding="utf-8",
    )


def main() -> None:
    metadata_rows: list[dict[str, object]] = []
    _write_readme()
    for panel in PANELS:
        source_path = DATA1_MAIN_IMAGE_ROOT / panel.image_name
        image = cv2.imread(str(source_path))
        if image is None:
            raise FileNotFoundError(f"Missing source image: {source_path}")

        x0, y0, x1, y1 = panel.crop_box
        crop = image[y0:y1, x0:x1]
        hsv = cv2.cvtColor(crop, cv2.COLOR_BGR2HSV)

        rows: list[dict[str, object]] = []
        for family in COLOR_RANGES:
            points_px = _resolve_points(_mask_for_family(hsv, family))
            if points_px.size == 0:
                continue
            points_data = _pixel_to_data(points_px, panel)
            for (xp, yp), (xd, yd) in zip(points_px, points_data, strict=True):
                rows.append(
                    {
                        "figure_id": panel.figure_id,
                        "panel_id": panel.panel_id,
                        "series_family": family,
                        "x_data": round(float(xd), 6),
                        "y_data": round(float(yd), 6),
                        "x_pixel": int(xp),
                        "y_pixel": int(yp),
                        "source_image": panel.image_name,
                        "notes": panel.notes,
                    }
                )
        _write_points_csv(panel, rows)
        metadata_rows.append(
            {
                "figure_id": panel.figure_id,
                "panel_id": panel.panel_id,
                "source_image": panel.image_name,
                "crop_x0": x0,
                "crop_y0": y0,
                "crop_x1": x1,
                "crop_y1": y1,
                "x_min": panel.x_min,
                "x_max": panel.x_max,
                "y_min": panel.y_min,
                "y_max": panel.y_max,
                "notes": panel.notes,
            }
        )
    _write_metadata(metadata_rows)
    print(f"[digitize-data1-main] wrote panel CSVs under: {OUT_ROOT}")


if __name__ == "__main__":
    main()
