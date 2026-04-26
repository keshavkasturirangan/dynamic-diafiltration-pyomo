#!/usr/bin/env python3
"""Generate cleaner series-level CSVs for DATA1 main paper Figures 2 and 3."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import cv2
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[5]
IMAGE_ROOT = (
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


PANELS = {
    "fig2A": PanelSpec("fig2", "A", "img-004.jpg", (160, 56, 573, 615), 0.0, 40.0, 0.0, 0.8),
    "fig2B": PanelSpec("fig2", "B", "img-004.jpg", (910, 56, 1490, 615), 0.0, 40.0, 0.0, 8.2),
    "fig2C": PanelSpec("fig2", "C", "img-004.jpg", (160, 775, 733, 1322), 0.0, 260.0, 0.0, 1.25),
    "fig2D": PanelSpec("fig2", "D", "img-004.jpg", (910, 775, 1492, 1322), 0.0, 260.0, 0.0, 85.0),
    "fig3A": PanelSpec("fig3", "A", "img-005.jpg", (130, 23, 826, 786), 0.0, 90.0, 0.0, 26.0),
}


def _family_mask(hsv: np.ndarray, family: str) -> np.ndarray:
    ranges = {
        "red": [((0, 70, 70), (12, 255, 255)), ((170, 70, 70), (180, 255, 255))],
        "orange": [((10, 70, 70), (28, 255, 255))],
        "green": [((40, 40, 50), (90, 255, 255))],
        "cyan": [((80, 40, 50), (105, 255, 255))],
        "blue": [((100, 40, 50), (136, 255, 255))],
        "magenta": [((136, 40, 50), (170, 255, 255))],
    }
    mask = np.zeros(hsv.shape[:2], dtype=np.uint8)
    for lower, upper in ranges[family]:
        mask |= cv2.inRange(hsv, np.array(lower, dtype=np.uint8), np.array(upper, dtype=np.uint8))
    return mask


def _pixel_to_data(panel: PanelSpec, x_px: np.ndarray, y_px: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x0, y0, x1, y1 = panel.crop_box
    width = max(x1 - x0 - 1, 1)
    height = max(y1 - y0 - 1, 1)
    x_data = panel.x_min + (x_px / width) * (panel.x_max - panel.x_min)
    y_data = panel.y_max - (y_px / height) * (panel.y_max - panel.y_min)
    return x_data, y_data


def _components(mask: np.ndarray, *, min_area: int, max_area: int | None = None) -> list[dict[str, float]]:
    num, labels, stats, centroids = cv2.connectedComponentsWithStats((mask > 0).astype("uint8"), 8)
    rows: list[dict[str, float]] = []
    for idx in range(1, num):
        area = int(stats[idx, cv2.CC_STAT_AREA])
        if area < min_area:
            continue
        if max_area is not None and area > max_area:
            continue
        x, y, w, h = [int(v) for v in stats[idx, :4]]
        cx, cy = [float(v) for v in centroids[idx]]
        rows.append({"area": area, "x": cx, "y": cy, "bbox_x": x, "bbox_y": y, "bbox_w": w, "bbox_h": h})
    return rows


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _sample_curve(mask: np.ndarray, panel: PanelSpec, series_name: str) -> list[dict[str, object]]:
    ys, xs = np.where(mask > 0)
    rows: list[dict[str, object]] = []
    if xs.size == 0:
        return rows
    for x in sorted(np.unique(xs)):
        y_values = ys[xs == x]
        if y_values.size < 2:
            continue
        y = float(np.median(y_values))
        x_data, y_data = _pixel_to_data(panel, np.array([x], dtype=float), np.array([y], dtype=float))
        rows.append(
            {
                "figure_id": panel.figure_id,
                "panel_id": panel.panel_id,
                "series_name": series_name,
                "point_index": len(rows),
                "x_data": round(float(x_data[0]), 6),
                "y_data": round(float(y_data[0]), 6),
                "x_pixel": int(x),
                "y_pixel": round(float(y), 3),
            }
        )
    return rows


def _cluster(values: list[float], gap: float) -> list[list[float]]:
    if not values:
        return []
    clusters: list[list[float]] = [[values[0]]]
    for value in values[1:]:
        if value - clusters[-1][-1] <= gap:
            clusters[-1].append(value)
        else:
            clusters.append([value])
    return clusters


def _blue_series_rows(panel: PanelSpec, hsv: np.ndarray, gap: float) -> list[dict[str, object]]:
    blue = _family_mask(hsv, "blue")
    comps = _components(blue, min_area=120)
    filtered = [c for c in comps if c["bbox_h"] >= 18 and c["bbox_w"] <= 80]
    x_centers = sorted(c["x"] for c in filtered)
    clusters = _cluster(x_centers, gap=gap)
    rows: list[dict[str, object]] = []
    ys, xs = np.where(blue > 0)
    if xs.size == 0:
        return rows
    for series_idx, cluster in enumerate(clusters, start=1):
        lo = min(cluster) - gap / 2
        hi = max(cluster) + gap / 2
        keep = (xs >= lo) & (xs <= hi)
        x_keep = xs[keep]
        y_keep = ys[keep]
        if x_keep.size < 20:
            continue
        coeff = np.polyfit(y_keep.astype(float), x_keep.astype(float), deg=1)
        y_min = float(np.min(y_keep))
        y_max = float(np.max(y_keep))
        y_samples = np.linspace(y_min, y_max, 80)
        x_samples = coeff[0] * y_samples + coeff[1]
        x_data, y_data = _pixel_to_data(panel, x_samples, y_samples)
        for point_idx, (xp, yp, xd, yd) in enumerate(zip(x_samples, y_samples, x_data, y_data, strict=True)):
            rows.append(
                {
                    "figure_id": panel.figure_id,
                    "panel_id": panel.panel_id,
                    "series_name": f"blue_line_{series_idx}",
                    "point_index": point_idx,
                    "x_data": round(float(xd), 6),
                    "y_data": round(float(yd), 6),
                    "x_pixel": round(float(xp), 3),
                    "y_pixel": round(float(yp), 3),
                }
            )
    return rows


def _marker_rows(panel: PanelSpec, hsv: np.ndarray, family: str, *, min_area: int, max_area: int | None = None, min_w: int = 0, min_h: int = 0, opened_kernel: int | None = None) -> list[dict[str, object]]:
    mask = _family_mask(hsv, family)
    if opened_kernel is not None:
        kernel = np.ones((opened_kernel, opened_kernel), np.uint8)
        mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
    comps = _components(mask, min_area=min_area, max_area=max_area)
    rows: list[dict[str, object]] = []
    kept = 0
    for comp in sorted(comps, key=lambda item: (item["x"], item["y"])):
        if comp["bbox_w"] < min_w or comp["bbox_h"] < min_h:
            continue
        x_data, y_data = _pixel_to_data(panel, np.array([comp["x"]]), np.array([comp["y"]]))
        rows.append(
            {
                "figure_id": panel.figure_id,
                "panel_id": panel.panel_id,
                "series_name": f"{family}_marker",
                "point_index": kept,
                "x_data": round(float(x_data[0]), 6),
                "y_data": round(float(y_data[0]), 6),
                "x_pixel": round(float(comp["x"]), 3),
                "y_pixel": round(float(comp["y"]), 3),
                "area_px": int(comp["area"]),
                "bbox_w": int(comp["bbox_w"]),
                "bbox_h": int(comp["bbox_h"]),
            }
        )
        kept += 1
    return rows


def _write_fig2_refined() -> None:
    source = cv2.imread(str(IMAGE_ROOT / "img-004.jpg"))
    for key in ("fig2A", "fig2B", "fig2C", "fig2D"):
        panel = PANELS[key]
        x0, y0, x1, y1 = panel.crop_box
        crop = source[y0:y1, x0:x1]
        hsv = cv2.cvtColor(crop, cv2.COLOR_BGR2HSV)

        trace_rows: list[dict[str, object]] = []
        marker_rows: list[dict[str, object]] = []

        if panel.panel_id in {"A", "C"}:
            gap = 38.0 if panel.panel_id == "A" else 48.0
            trace_rows.extend(_blue_series_rows(panel, hsv, gap=gap))
            max_area = 260 if panel.panel_id == "A" else 320
            marker_rows.extend(_marker_rows(panel, hsv, "red", min_area=8, max_area=max_area))
        else:
            trace_rows.extend(_sample_curve(_family_mask(hsv, "green"), panel, "green_trace"))
            trace_rows.extend(_sample_curve(_family_mask(hsv, "red"), panel, "red_trace"))
            marker_rows.extend(_marker_rows(panel, hsv, "magenta", min_area=120, min_w=14, min_h=10))
            marker_rows.extend(_marker_rows(panel, hsv, "cyan", min_area=35, min_w=5, min_h=5))
            marker_rows.extend(_marker_rows(panel, hsv, "red", min_area=150, min_w=12, min_h=12, opened_kernel=9))

        fig_dir = OUT_ROOT / panel.figure_id
        _write_csv(
            fig_dir / f"data1_main_{panel.figure_id}_{panel.panel_id.lower()}_traces_refined.csv",
            ["figure_id", "panel_id", "series_name", "point_index", "x_data", "y_data", "x_pixel", "y_pixel"],
            trace_rows,
        )
        _write_csv(
            fig_dir / f"data1_main_{panel.figure_id}_{panel.panel_id.lower()}_markers_refined.csv",
            ["figure_id", "panel_id", "series_name", "point_index", "x_data", "y_data", "x_pixel", "y_pixel", "area_px", "bbox_w", "bbox_h"],
            marker_rows,
        )


def _write_fig3_refined() -> None:
    panel = PANELS["fig3A"]
    source = cv2.imread(str(IMAGE_ROOT / panel.image_name))
    x0, y0, x1, y1 = panel.crop_box
    crop = source[y0:y1, x0:x1]
    hsv = cv2.cvtColor(crop, cv2.COLOR_BGR2HSV)

    rows: list[dict[str, object]] = []
    family_specs = {
        "green": {"min_area": 150},
        "blue": {"min_area": 150},
        "orange": {"min_area": 150},
        "cyan": {"min_area": 90},
        "magenta": {"min_area": 400},
        "red": {"min_area": 250},
    }
    for family, params in family_specs.items():
        for row in _marker_rows(panel, hsv, family, **params, min_w=8, min_h=8):
            row["series_name"] = family
            rows.append(row)

    _write_csv(
        OUT_ROOT / panel.figure_id / "data1_main_fig3_a_markers_refined.csv",
        ["figure_id", "panel_id", "series_name", "point_index", "x_data", "y_data", "x_pixel", "y_pixel", "area_px", "bbox_w", "bbox_h"],
        rows,
    )


def main() -> None:
    _write_fig2_refined()
    _write_fig3_refined()
    print(f"[refine-data1-main] wrote refined Figure 2/3 CSVs under: {OUT_ROOT}")


if __name__ == "__main__":
    main()
