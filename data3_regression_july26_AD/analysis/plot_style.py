"""
Shared publication-quality plotting style for all DATA3 figures.

Follows the guidance in
https://ndcbe.github.io/data-and-computing/notebooks/01/Publication-Quality-Figures.html
Every figure-producing script in this project should call
`apply_style()` once (near the top, after `matplotlib.use("Agg")`) and use
`save_fig()` to write each figure, so every plot in the report shares one
visual language: font sizes, line widths, inward ticks on all four sides,
bold axis labels (with math bolded via `\boldsymbol{}`, see below) and
bracketed units, boxed legends, a colorblind-safe palette, and PNG+PDF
output.

This module only affects appearance -- it has no dependence on, and no
effect on, any numerical results.

Figure *physical* size matters as much as font size: LaTeX's
\\includegraphics[width=...] rescales the whole figure (fonts included) by
(display_width / saved_figure_width). If a figure is saved narrower than
its eventual display width, its fonts render LARGER on the page than the
rcParams point sizes suggest, and vice versa. To keep on-page figure text
close to the report's 11pt body text, every figure should be *saved* at
(approximately) its true on-page display width -- see `fig_width()` below,
measured directly from docs/reports/main.tex (11pt article, 1in margins):
    \\the\\linewidth        = 469.755 pt = 6.50 in   (full-width figures)
    \\the\\linewidth * 0.49 = 230.183 pt = 3.19 in   (0.49\\textwidth minipage
                                                       two-up figures)
When `\\includegraphics` uses some OTHER width (e.g. `0.6\\linewidth` for a
tall/narrow figure), pass that same fraction into the figure's *saved* width
too (`fig_width("full") * 0.6`) rather than leaving it at a mismatched
default -- see `docs/prompts/figure_compliance_audit.md` for the full
checklist this module implements.

Legends are boxed by default (`apply_style()` sets `legend.frameon=True`
with a light gray edge and an opaque face) -- this SUPERSEDES the project's
earlier frameon=False convention (see git history / prior prompts), adopted
because unboxed legends were hard to distinguish from data at a glance in
multi-panel figures. Every legend in every figure should render with a
visible box; there is no per-call opt-out needed since it is a global
rcParams default -- if a script explicitly passes `frameon=False` to a
`legend()`/`fig.legend()` call, remove that override.

Bold axis labels/titles do NOT automatically bold embedded mathtext:
matplotlib's bold rcParams (`axes.labelweight`, or manual `fontweight="bold"`)
only bolds the plain-text parts of a label -- `r"$|\\chi|$ [mM]"` next to a
bold `[mM]` still renders `|\\chi|` in normal weight. Two additional facts
that trip this up:
  - `\\mathbf{}` does NOT bold lowercase Greek letters (chi, delta, ...) in
    matplotlib's mathtext engine -- it silently no-ops on them.
  - `\\boldsymbol{}` DOES bold both Latin and Greek mathtext, and is the
    correct choice everywhere in this project.
So every axis label / title containing inline math should wrap the math in
`\\boldsymbol{}`, e.g. `r"$\\boldsymbol{|\\chi|}$ [mM]"`,
`r"$\\boldsymbol{\\delta^*}$ [--]"`, `r"$\\Delta \\boldsymbol{c_m}$ [mM]"` --
NOT `r"$\\mathbf{|\\chi|}$ [mM]"` (silently non-bold for the Greek/symbol
parts) and NOT a plain `r"$|\\chi|$ [mM]"` (math left un-bolded next to a
bold unit). Plain non-math annotations (e.g. plot titles with no formulas)
already bold correctly via rcParams/`fontweight` and need no `\\boldsymbol{}`.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

# Measured from docs/reports/main.tex (11pt article, geometry margin=1in):
# \the\linewidth = 469.755pt; the two-up figures sit in 0.49\textwidth
# minipages = 0.49 * 469.755pt = 230.183pt. 72.27 pt/in (TeX point).
FULL_WIDTH_IN = 469.755 / 72.27   # ~6.50 in
HALF_WIDTH_IN = 230.183 / 72.27   # ~3.19 in

# Okabe--Ito colorblind-safe palette (Okabe & Ito, 2008), the standard
# palette recommended for scientific figures.
OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
}

# Ordered palette used for assigning one color per dataset/series, in a
# fixed, reproducible order (first color is reserved for "primary"/single
# -dataset plots elsewhere in the project, so multi-dataset figures start
# from the second color to stay visually distinct from those).
#
# 12 entries: the 8 Okabe-Ito colors plus 4 extra hand-picked, visually
# distinct colors, so figures with up to 12 categories (e.g. the 9 NaCl
# datasets) never wrap around and collide (an 8-color palette repeating at
# index 8 caused two datasets to render in the same blue -- see
# docs/prompts/nacl_report_refinement2.md item 3).
PALETTE = [
    OKABE_ITO["blue"],
    OKABE_ITO["vermillion"],
    OKABE_ITO["bluish_green"],
    OKABE_ITO["orange"],
    OKABE_ITO["reddish_purple"],
    OKABE_ITO["sky_blue"],
    OKABE_ITO["black"],
    OKABE_ITO["yellow"],
    "#8B4513",  # saddle brown
    "#999999",  # medium gray
    "#4B0082",  # indigo
    "#2E8B57",  # sea green
]


def apply_style() -> None:
    """Set matplotlib rcParams for consistent, publication-quality figures.

    Call once per process, before creating any figures.
    """
    matplotlib.rcParams.update(
        {
            "figure.figsize": (FULL_WIDTH_IN, FULL_WIDTH_IN),
            "figure.dpi": 100,  # on-screen/interactive; savefig.dpi controls file dpi
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            # Sized to match the ~11pt report body text once each figure is
            # saved at its true display width via fig_width() (see module
            # docstring) and included at width=\linewidth in main.tex.
            "font.size": 9,
            "axes.labelsize": 11,
            "axes.labelweight": "bold",
            "axes.titlesize": 10,
            "axes.titleweight": "bold",
            "figure.titleweight": "bold",
            "axes.linewidth": 1.0,
            "lines.linewidth": 2,
            "lines.markersize": 5,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            "xtick.minor.visible": False,
            "ytick.minor.visible": False,
            "legend.fontsize": 9,
            "legend.title_fontsize": 9,
            # Boxed legends (supersedes the earlier frameon=False convention,
            # see module docstring): a light, visible frame around every
            # legend so it reads clearly against data in busy multi-panel
            # figures.
            "legend.frameon": True,
            "legend.framealpha": 1.0,
            "legend.edgecolor": "#888888",
            "legend.borderaxespad": 0.5,
        }
    )


def fig_width(kind: str = "full") -> float:
    """Return the figure width (in inches) to *save* a figure at, so that
    when included in main.tex at width=\\linewidth (full) or inside a
    0.49\\textwidth minipage (half), the display scale is ~1 and figure
    text renders at approximately its rcParams point size (~body-text
    size), rather than being magnified or shrunk by LaTeX's rescaling.
    """
    widths = {"full": FULL_WIDTH_IN, "half": HALF_WIDTH_IN}
    if kind not in widths:
        raise ValueError(f"kind must be one of {list(widths)}, got {kind!r}")
    return widths[kind]


def get_dataset_colors(labels) -> dict:
    """Map an ordered sequence of dataset labels to a fixed, reproducible
    color from PALETTE, so the same dataset gets the same color across
    every figure in a script (and across scripts, if the label order is
    kept consistent)."""
    return {label: PALETTE[i % len(PALETTE)] for i, label in enumerate(labels)}


def save_fig(fig, png_path) -> None:
    """Save `fig` as both PNG (dpi from rcParams, i.e. 300) and a PDF
    sibling with the same stem, both with a tight bounding box."""
    png_path = Path(png_path)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, bbox_inches="tight")
    fig.savefig(png_path.with_suffix(".pdf"), bbox_inches="tight")
