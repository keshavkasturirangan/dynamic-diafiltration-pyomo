#!/usr/bin/env python
"""Build a clean combined DATA3 collaborator deck (2026-06-11).

Combines:
  - parameter-bounds confirmation (from contour-diagnostics v6 slide 2)
  - sigma x Lp contour maps of all parameters / identifiability 3-stage story (v5/v6)
  - mass + concentration fits w/ contour-consistent theta (B_form=1 -> curved cF)
    from meeting_prep_2026-06-10/consistent_initial_guess

No model re-run: existing PNGs are embedded.  Clean legend + corrected
concentration-axis caption are added natively on each fit slide.
"""
import csv
import re
from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from pptx.oxml.ns import qn

# ----------------------------------------------------------------------------
REPO = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo")
NF270 = REPO / "UnifiedFramework/DATA3/results/paper_artifacts/nf270"
FIT_DIR = NF270 / "meeting_prep_2026-06-10/consistent_initial_guess"
FITS_FIXED = NF270 / "meeting_prep_2026-06-10/fits_fixed"  # startup-vial-corrected, DATA2-style
CONTOUR_DIR = NF270 / "matlab_ports_5ch_full"
ANIM_DIR = NF270 / "animations3d"
ALLPAIRS = NF270 / "matlab_ports_5ch_allpairs/MC2.05.07.24_NaCl"
LOGO = REPO / "refactored_codes_v1/notre_dame_logo.png"
OUT = REPO / "refactored_codes_v1/DATA3_followup_analysis_2026-06-11.pptx"

# palette --------------------------------------------------------------------
NAVY = RGBColor(0x0C, 0x23, 0x40)
NAVY2 = RGBColor(0x12, 0x33, 0x5B)
GOLD = RGBColor(0xC9, 0x97, 0x00)
INK = RGBColor(0x1B, 0x1B, 0x1B)
GRAY = RGBColor(0x5B, 0x66, 0x70)
LGRAY = RGBColor(0x8A, 0x93, 0x9B)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
PANEL = RGBColor(0xF4, 0xF6, 0xF8)
BOXBG = RGBColor(0xFB, 0xF7, 0xEA)   # soft cream callout
S1 = RGBColor(0x2E, 0x7D, 0x32)      # well-identified / GOOD (green)
S2 = RGBColor(0xD9, 0x8A, 0x00)      # inconsistent / OKAY (amber)
S3 = RGBColor(0xB0, 0x41, 0x3E)      # not estimable / POOR (red)

HFONT = "Calibri"
BFONT = "Calibri"

EMU_IN = 914400
prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
SW, SH = 13.333, 7.5
BLANK = prs.slide_layouts[6]


# ---- helpers ---------------------------------------------------------------
def slide():
    return prs.slides.add_slide(BLANK)


def bg(s, color):
    s.background.fill.solid()
    s.background.fill.fore_color.rgb = color


def _set_font(run, size, color, bold=False, italic=False, font=BFONT):
    run.font.size = Pt(size)
    run.font.color.rgb = color
    run.font.bold = bold
    run.font.italic = italic
    run.font.name = font


def textbox(s, x, y, w, h, lines, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP,
            wrap=True, space_after=2):
    """lines: list of list-of-run-dicts (each run: text,size,color,bold,italic,font).
    A plain string line == one default run."""
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = wrap
    tf.vertical_anchor = anchor
    for m in ("margin_left", "margin_right", "margin_top", "margin_bottom"):
        setattr(tf, m, Inches(0.04))
    for i, line in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(space_after)
        p.space_before = Pt(0)
        if isinstance(line, str):
            line = [{"text": line}]
        for rd in line:
            r = p.add_run()
            r.text = rd.get("text", "")
            _set_font(r, rd.get("size", 14), rd.get("color", INK),
                      rd.get("bold", False), rd.get("italic", False),
                      rd.get("font", BFONT))
    return tb


def rect(s, x, y, w, h, fill, line=None, line_w=None, round_=False, shadow=False):
    shp = s.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE if round_ else MSO_SHAPE.RECTANGLE,
        Inches(x), Inches(y), Inches(w), Inches(h))
    if fill is None:
        shp.fill.background()
    else:
        shp.fill.solid()
        shp.fill.fore_color.rgb = fill
    if line is None:
        shp.line.fill.background()
    else:
        shp.line.color.rgb = line
        shp.line.width = Pt(line_w or 1)
    shp.shadow.inherit = False
    if round_:
        try:
            shp.adjustments[0] = 0.10
        except Exception:
            pass
    return shp


def img_fit(path, x, y, w):
    """place image scaled to width w, return (w,h) used."""
    from PIL import Image
    iw, ih = Image.open(path).size
    h = w * ih / iw
    return x, y, w, h


def add_image(s, path, x, y, w):
    x, y, w, h = img_fit(path, x, y, w)
    s.shapes.add_picture(str(path), Inches(x), Inches(y), Inches(w), Inches(h))
    return h


# ---- LaTeX-style equation rendering (matplotlib mathtext -> transparent PNG) -
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as _plt

EQ_DIR = REPO / "refactored_codes_v1/_deck_eq"
EQ_DIR.mkdir(exist_ok=True)
_eq_n = [0]


def render_eq(tex, color=INK, fontsize=26, dpi=300):
    from PIL import Image
    path = EQ_DIR / f"eq_{_eq_n[0]}.png"
    _eq_n[0] += 1
    fig = _plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, tex, fontsize=fontsize, color="#" + str(color))
    fig.savefig(path, dpi=dpi, transparent=True, bbox_inches="tight", pad_inches=0.02)
    _plt.close(fig)
    return path, Image.open(path).size


# ---- native OMML equations (editable in PowerPoint) ------------------------
import sys as _sys
_sys.path.insert(0, str(REPO / "refactored_codes_v1"))
from _eq2omml import latex_to_omath
from pptx.oxml import parse_xml
from xml.sax.saxutils import escape as _xesc


def place_omath(s, latex, x, y, w, sz=20, color=INK, fallback="", align=PP_ALIGN.LEFT, h=0.6):
    """Insert a NATIVE, editable PowerPoint equation (OMML) in a textbox.
    `fallback` is clean unicode shown by non-PowerPoint renderers (e.g. LibreOffice)."""
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    for m in ("margin_left", "margin_right", "margin_top", "margin_bottom"):
        setattr(tf, m, Inches(0.02))
    p = tf.paragraphs[0]
    p.alignment = align
    omath = latex_to_omath(latex)
    hexc = str(color)
    sz100 = int(round(sz * 100))
    algn = {PP_ALIGN.LEFT: "l", PP_ALIGN.CENTER: "ctr", PP_ALIGN.RIGHT: "r"}.get(align, "l")
    ppr = (f'<a:pPr xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" algn="{algn}">'
           f'<a:defRPr sz="{sz100}"><a:solidFill><a:srgbClr val="{hexc}"/></a:solidFill></a:defRPr></a:pPr>')
    ac = (
        '<mc:AlternateContent xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006">'
        '<mc:Choice xmlns:a14="http://schemas.microsoft.com/office/drawing/2010/main" Requires="a14">'
        '<a14:m xmlns:m="http://schemas.openxmlformats.org/officeDocument/2006/math">'
        f'<m:oMathPara>{omath}</m:oMathPara></a14:m></mc:Choice>'
        '<mc:Fallback xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
        f'<a:r><a:rPr lang="en-US" sz="{sz100}" dirty="0"><a:solidFill><a:srgbClr val="{hexc}"/></a:solidFill>'
        f'<a:latin typeface="Cambria Math"/></a:rPr><a:t>{_xesc(fallback)}</a:t></a:r></mc:Fallback>'
        '</mc:AlternateContent>')
    p._p.insert(0, parse_xml(ppr))
    p._p.append(parse_xml(ac))
    return tb


def place_eq(s, tex, x, y, height, color=INK, fontsize=26, align="left", box_w=None):
    """Render `tex` (mathtext, no surrounding $) and place at (x,y) scaled to `height` inches.
    If align != left and box_w given, x is the left of a box of width box_w to align within."""
    path, (iw, ih) = render_eq(tex, color, fontsize)
    w = height * iw / ih
    if box_w is not None and w > box_w:          # clamp to box width, shrink height
        height = height * box_w / w
        w = box_w
    if box_w is not None:
        if align == "center":
            x = x + (box_w - w) / 2
        elif align == "right":
            x = x + (box_w - w)
    s.shapes.add_picture(str(path), Inches(x), Inches(y), height=Inches(height))
    return w


def footer(s, idx, total):
    textbox(s, 0.55, 7.12, 8.5, 0.3,
            [[{"text": "DATA3 single-salt · identifiability & fits · collaborator meeting · June 11 2026",
               "size": 9, "color": LGRAY}]])
    textbox(s, 11.4, 7.12, 1.4, 0.3,
            [[{"text": f"{idx} / {total}", "size": 9, "color": LGRAY}]],
            align=PP_ALIGN.RIGHT)


def logo(s, dark=False):
    # The only logo asset is a wide dark wordmark (aspect 4.3) that overflows the
    # right margin and collides with stage chips, so it is intentionally omitted
    # on content slides — the footer carries the identity line. Drop in a square /
    # white ND mark here later if desired.
    return


def stage_chip(s, label, color, x=10.9, y=0.34, w=1.9, h=0.42):
    r = rect(s, x, y, w, h, color, round_=True)
    tf = r.text_frame
    tf.word_wrap = False
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    run = p.add_run()
    run.text = label
    _set_font(run, 11.5, WHITE, bold=True, font=HFONT)
    return r


def title_block(s, kicker, title, subtitle):
    if kicker:
        textbox(s, 0.55, 0.30, 9.5, 0.3,
                [[{"text": kicker, "size": 12, "color": GOLD, "bold": True,
                   "font": HFONT}]])
    textbox(s, 0.55, 0.58, 11.0, 0.7,
            [[{"text": title, "size": 26, "color": NAVY, "bold": True,
               "font": HFONT}]])
    if subtitle:
        textbox(s, 0.55, 1.18, 12.2, 0.45,
                [[{"text": subtitle, "size": 13.5, "color": GRAY,
                   "italic": True}]])


# ---- load theta ------------------------------------------------------------
THETA = {}
with open(FIT_DIR / "summary_contour_consistent_theta.csv") as f:
    for row in csv.DictReader(f):
        THETA[row["sheet"]] = row


def theta_line(sheet):
    r = THETA[sheet]
    Lp = float(r["Lp_contour"]); B = float(r["B_warm"])
    sig = float(r["sigma_contour"]); w3 = float(r["WSSE_3ch"])
    sflag = " ↑" if sig >= 0.999 else (" ↓" if sig <= 0.001 else "")
    bflag = " ↑" if B >= 29.99 else ""
    return Lp, B, sig, w3, sflag, bflag


def theta_math(sheet, wsse=True):
    """Native-math LaTeX for the contour-consistent θ readout."""
    Lp, B, sig, w3, sflag, bflag = theta_line(sheet)
    parts = [rf"L_p = {Lp:.2f}", rf"B = {B:.2f}\ \mu\mathrm{{m/s}}", rf"\sigma = {sig:.2f}"]
    eq = r",\quad ".join(parts)
    if wsse:
        eq += rf",\quad \mathrm{{WSSE}}_{{3ch}} = {w3:,.0f}".replace(",", "{,}")
    return eq


def theta_fallback(sheet, wsse=True):
    Lp, B, sig, w3, sflag, bflag = theta_line(sheet)
    t = f"Lp = {Lp:.2f},  B = {B:.2f} µm/s,  σ = {sig:.2f}{sflag}"
    if wsse:
        t += f",  WSSE_3ch = {w3:,.0f}"
    return t


def short(sheet):
    # MC3.07.22.24_SNaCl -> ("MC3", "07.22.24", "SNaCl")
    mc, rest = sheet.split(".", 1)
    date, salt = rest.rsplit("_", 1)
    # date is like 07.22.24 ; sometimes mc already MC3
    return mc, date, salt


def pretty(sheet):
    mc, date, salt = short(sheet)
    salt = (salt.replace("CaCl2", "CaCl₂").replace("LaCl3", "LaCl₃"))
    return f"{mc} · {date} · {salt}"


# ===========================================================================
# slide registry (built first so we can number with a known total)
# ===========================================================================
slides_built = []
TOTAL = 0  # patched at end


# legend (native) for fit slides --------------------------------------------
def fit_legend(s, x, y, w):
    """clean shared legend for the mass + concentration panels."""
    rect(s, x, y, w, 1.02, PANEL, line=RGBColor(0xD7, 0xDC, 0xE1), line_w=0.75,
         round_=True)
    textbox(s, x + 0.15, y + 0.06, w - 0.3, 0.3,
            [[{"text": "How to read the panels", "size": 11.5, "color": NAVY,
               "bold": True, "font": HFONT}]])
    rows = [
        [{"text": "● ", "size": 11, "color": INK, "bold": True},
         {"text": "experiment   ", "size": 11, "color": INK},
         {"text": "—— ", "size": 11, "color": INK, "bold": True},
         {"text": "model fit", "size": 11, "color": INK},
         {"text": "      colour = vial index (sequential)", "size": 11, "color": GRAY}],
        [{"text": "mass panel: ", "size": 11, "color": NAVY, "bold": True},
         {"text": "cumulative permeate mass per vial vs time", "size": 11, "color": INK}],
        [{"text": "conc. panel: ", "size": 11, "color": NAVY, "bold": True},
         {"text": "retentate c", "size": 11, "color": INK},
         {"text": "f", "size": 8, "color": INK},
         {"text": " (● exp / — model)   permeate c", "size": 11, "color": INK},
         {"text": "h", "size": 8, "color": INK},
         {"text": " (▲ exp / –– model, dashed)", "size": 11, "color": INK}],
    ]
    textbox(s, x + 0.15, y + 0.33, w - 0.3, 0.66, rows, space_after=1)


def contour_legend(s, x, y, w, h=0.44):
    """Slim single-strip legend for the contour / scrub panels."""
    rect(s, x, y, w, h, PANEL, line=RGBColor(0xD7, 0xDC, 0xE1), line_w=0.75, round_=True)
    textbox(s, x + 0.15, y + 0.02, w - 0.3, h - 0.02, [[
        {"text": "Legend:   ", "size": 10.5, "color": NAVY, "bold": True, "font": HFONT},
        {"text": "▲", "size": 10.5, "color": S3, "bold": True},
        {"text": " = per-panel minimum (argmin)    ", "size": 10.5, "color": INK},
        {"text": "contour lines", "size": 10.5, "color": INK, "bold": True},
        {"text": " = log₁₀ weighted SSR (cooler → better fit)    ", "size": 10.5, "color": INK},
        {"text": "each panel = one measurement channel    panels agree → identifiable", "size": 10.5, "color": INK},
    ]], anchor=MSO_ANCHOR.MIDDLE)


def verdict_chip(s, kind, x, y, w=2.3):
    cmap = {"GOOD": (S1, "Identifiability: GOOD"),
            "OKAY": (S2, "Identifiability: OKAY"),
            "POOR": (S3, "Identifiability: POOR")}
    c, lab = cmap[kind]
    r = rect(s, x, y, w, 0.4, c, round_=True)
    tf = r.text_frame; tf.word_wrap = False
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    run = p.add_run(); run.text = lab
    _set_font(run, 11, WHITE, bold=True, font=HFONT)


# === BUILD ==================================================================

# 1 — TITLE -----------------------------------------------------------------
def s_title():
    s = slide(); bg(s, NAVY)
    rect(s, 0, 6.95, SW, 0.55, GOLD)
    textbox(s, 0.9, 2.05, 11.5, 0.4,
            [[{"text": "DATA3  ·  NF270 SINGLE-SALT CAMPAIGN", "size": 15,
               "color": GOLD, "bold": True, "font": HFONT}]])
    textbox(s, 0.9, 2.5, 11.5, 1.5,
            [[{"text": "Identifiability & Analysis — Follow-up",
               "size": 40, "color": WHITE, "bold": True, "font": HFONT}]])
    textbox(s, 0.9, 3.55, 11.5, 0.6,
            [[{"text": "Paper-justified measurement error  ·  σ-identifiability sensitivity  ·  "
                       "axis convention  ·  B(ionic strength) & per-band θ", "size": 16,
               "color": RGBColor(0xCC, 0xD6, 0xE2)}]])
    textbox(s, 0.9, 5.6, 11.5, 0.9,
            [[{"text": "Keshav Kasturi Rangan", "size": 16, "color": WHITE, "bold": True}],
             [{"text": "Post-meeting analysis  ·  June 11, 2026", "size": 13,
               "color": RGBColor(0xAE, 0xBC, 0xCE)}]])
    return s


# 2 — AGENDA / WHAT'S NEW ----------------------------------------------------
def s_agenda():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "WHERE WE ARE",
                "What's new since the last collaborator deck",
                "Four threads, building on the single-salt fits you saw last time (v9)")
    items = [
        ("01", "Tighten the parameter bounds", S2,
         "Per-salt B upper bound, halved by cation charge. Lₚ and σ bounds unchanged from the DATA1/DATA2 KCl baseline."),
        ("02", "Contour-map every parameter", NAVY,
         "Forward-simulate the full Spiegler–Kedem DAE on a σ×Lₚ grid; weighted SSR per channel. No fitting — a pure identifiability diagnostic across 5 measurement channels."),
        ("03", "Diagnose identifiability: good / okay / poor", GOLD,
         "Three stages across all 11 sheets — well-identified, inconsistent, not-estimable — and how each reflects in the mass & concentration fit."),
        ("04", "Consistent contour-derived fits", S1,
         "Use the joint contour argmin (mass + perm-conc + ret-conc) as a consistent initial guess. Variable B (B_form = 1) → curved cₕ(t), not a straight line."),
    ]
    y = 1.95
    for num, head, col, body in items:
        rect(s, 0.7, y + 0.03, 0.62, 0.62, col, round_=True)
        textbox(s, 0.7, y + 0.12, 0.62, 0.45,
                [[{"text": num, "size": 17, "color": WHITE, "bold": True, "font": HFONT}]],
                align=PP_ALIGN.CENTER)
        textbox(s, 1.5, y - 0.02, 11.2, 0.4,
                [[{"text": head, "size": 16.5, "color": NAVY, "bold": True, "font": HFONT}]])
        textbox(s, 1.5, y + 0.36, 11.2, 0.7,
                [[{"text": body, "size": 12.5, "color": GRAY}]])
        y += 1.18
    return s


# 3 — PARAMETER BOUNDS -------------------------------------------------------
def s_bounds():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "01 · CONFIRM FIRST",
                "The parameter bounds I'm changing",
                "Lₚ and σ inherited from the DATA1/DATA2 KCl baseline — the per-salt B upper bound is the only new knob")
    # left: inherited baseline card
    rect(s, 0.6, 1.95, 5.4, 1.55, PANEL, line=RGBColor(0xD7, 0xDC, 0xE1), line_w=0.75, round_=True)
    textbox(s, 0.8, 2.06, 5.0, 0.35,
            [[{"text": "Inherited baseline (KCl · DATA1/DATA2)", "size": 13.5,
               "color": NAVY, "bold": True, "font": HFONT}]])
    textbox(s, 0.8, 2.45, 5.0, 1.0,
            [[{"text": "Lₚ  ∈  [0.5, 50]", "size": 14, "color": INK, "bold": True},
              {"text": "   L/(m²·hr·bar)", "size": 12, "color": GRAY}],
             [{"text": "σ   ∈  [0, 1]", "size": 14, "color": INK, "bold": True},
              {"text": "   (dimensionless reflection coeff.)", "size": 12, "color": GRAY}],
             [{"text": "B   ∈  [0, 30]", "size": 14, "color": INK, "bold": True},
              {"text": "   µm/s   — unchanged for monovalent", "size": 12, "color": GRAY}]],
            space_after=4)
    # right: per-salt B table
    rect(s, 6.35, 1.95, 6.35, 1.55, BOXBG, line=GOLD, line_w=1, round_=True)
    textbox(s, 6.55, 2.06, 6.0, 0.35,
            [[{"text": "New: per-salt B upper bound — halve per +1 cation charge",
               "size": 13.5, "color": NAVY, "bold": True, "font": HFONT}]])
    tbl = [("NaCl", "[0, 30]", "monovalent — same as KCl"),
           ("CaCl₂", "[0, 15]", "divalent — halved"),
           ("LaCl₃", "[0, 10]", "trivalent — a third")]
    yy = 2.45
    for salt, b, basis in tbl:
        textbox(s, 6.55, yy, 1.5, 0.3, [[{"text": salt, "size": 13, "color": INK, "bold": True}]])
        textbox(s, 8.0, yy, 1.7, 0.3, [[{"text": b + " µm/s", "size": 13, "color": S2, "bold": True}]])
        textbox(s, 9.7, yy, 3.0, 0.3, [[{"text": basis, "size": 11.5, "color": GRAY}]])
        yy += 0.34
    # how contours are built
    rect(s, 0.6, 3.8, 12.13, 2.75, RGBColor(0xF8, 0xFA, 0xFC),
         line=RGBColor(0xD7, 0xDC, 0xE1), line_w=0.75, round_=True)
    textbox(s, 0.85, 3.92, 11.6, 0.35,
            [[{"text": "How the contour panels are built", "size": 15, "color": NAVY,
               "bold": True, "font": HFONT}]])
    textbox(s, 0.85, 4.3, 11.7, 2.1,
            [[{"text": "At every (σ, Lₚ) grid point we forward-simulate the full Spiegler–Kedem DAE",
               "size": 13.5, "color": INK},
              {"text": "  (B fixed at warm-start, S₀/S at warm-start)", "size": 12.5, "color": GRAY}],
             [{"text": "then evaluate the weighted sum-of-squared-residuals for each channel.",
               "size": 13.5, "color": INK}],
             [{"text": "  ", "size": 6, "color": INK}],
             [{"text": "5 channels per panel:  ", "size": 13.5, "color": NAVY, "bold": True},
              {"text": "mass  ·  permeate conc.  ·  retentate conc.  ·  permeate conductivity  ·  retentate conductivity",
               "size": 13, "color": INK}],
             [{"text": "  ", "size": 6, "color": INK}],
             [{"text": "This is a diagnostic, not a fit — ", "size": 13.5, "color": S1, "bold": True},
              {"text": "pure forward-simulation. Where the 5 channels agree on a minimum, the data identify the parameters; where they split, they don't.",
               "size": 13, "color": INK}]],
            space_after=3)
    return s


# 4 — 3-STAGE FRAMEWORK ------------------------------------------------------
def s_framework():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "02 · IDENTIFIABILITY",
                "Three stages of identifiability across all 11 sheets",
                "Agreement across the 5 measurement channels = trust. What we trust most → what we trust least.")
    cards = [
        (S1, "STAGE 1", "Well-identified", "GOOD",
         "2 sheets · NaCl",
         "Both are diluting NaCl runs that start at high cf (~95 mM). High concentration → large Δπ → σ is elucidated (DATA1 §4.3.2), and Lₚ is pinned by mass. MC3 SNaCl keeps an interior σ; MC5 S2NaCl's σ drifts but Lₚ holds."),
        (S2, "STAGE 2", "Inconsistent", "OKAY",
         "7 sheets · 4 NaCl + 3 CaCl₂",
         "Lₚ tight, but σ rail-pinned at a wall. Channels each find a minimum, but at different (σ, Lₚ) walls — a weighting split, not a model-form failure."),
        (S3, "STAGE 3", "Not estimable", "POOR",
         "2 sheets · LaCl₃",
         "Same wrong-model bias on both sheets; WSSE ~20× the NaCl gold standard. Constant-(σ,B) Spiegler–Kedem can't represent La³⁺ ion-pairing."),
    ]
    x = 0.6; w = 3.94; gap = 0.18
    for col, tag, name, verdict, count, body in cards:
        y = 2.0
        rect(s, x, y, w, 4.5, WHITE, line=col, line_w=1.5, round_=True)
        rect(s, x, y, w, 0.9, col, round_=True)
        rect(s, x, y + 0.55, w, 0.35, col)  # square off bottom of header band
        textbox(s, x + 0.2, y + 0.12, w - 0.4, 0.4,
                [[{"text": tag, "size": 13, "color": WHITE, "bold": True, "font": HFONT}]])
        textbox(s, x + 0.2, y + 0.46, w - 0.4, 0.4,
                [[{"text": name, "size": 19, "color": WHITE, "bold": True, "font": HFONT}]])
        textbox(s, x + 0.2, y + 1.05, w - 0.4, 0.35,
                [[{"text": count, "size": 12.5, "color": col, "bold": True}]])
        textbox(s, x + 0.2, y + 1.45, w - 0.4, 2.3,
                [[{"text": body, "size": 13, "color": INK}]])
        # verdict ribbon
        rect(s, x + 0.2, y + 3.95, w - 0.4, 0.42, col, round_=True)
        textbox(s, x + 0.2, y + 4.0, w - 0.4, 0.35,
                [[{"text": "fit quality: " + verdict, "size": 12.5, "color": WHITE, "bold": True}]],
                align=PP_ALIGN.CENTER)
        x += w + gap
    return s


# 5 — THETA SUMMARY TABLE ----------------------------------------------------
def s_table():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "02 · IDENTIFIABILITY",
                "Where the optimizer landed on every sheet",
                "Warm-start θ across the campaign, colour-coded by stage. ↑/↓ = parameter pinned at a bound.")
    order = [
        ("MC3.07.22.24_SNaCl", S1, "Stage 1A · interior σ — gold standard"),
        ("MC5.07.23.24_S2NaCl", S1, "Stage 1B · diluting, σ harmless"),
        ("MC4.07.11.24_SNaCl", S2, "Stage 1B-ish · diluting BUT poor WSSE ⚠"),
        ("MC2.05.07.24_NaCl", S2, "Stage 2 · σ wall, concentrating"),
        ("MC5.07.23.24_NaCl", S2, "Stage 2 · σ wall, concentrating"),
        ("MC5.07.23.24_SNaCl", S2, "Stage 2 · σ wall, concentrating"),
        ("MC2.05.07.24_CaCl2", S2, "Stage 2 · channels split by side"),
        ("MC3.07.11.24_SCaCl2", S2, "Stage 2 · σ at wall"),
        ("MC3.07.12.24_S2CaCl2", S2, "Stage 2 · σ AND B at walls"),
        ("MC2.05.21.24_LaCl3", S3, "Stage 3 · wrong model — all channels disagree"),
        ("MC4.07.11.24_SLaCl3", S3, "Stage 3 · wrong model — replicate"),
    ]
    cols = [("Sheet", 0.6, 3.0, PP_ALIGN.LEFT),
            ("Lₚ", 3.6, 0.95, PP_ALIGN.CENTER),
            ("B [µm/s]", 4.55, 1.15, PP_ALIGN.CENTER),
            ("σ", 5.7, 1.05, PP_ALIGN.CENTER),
            ("WSSE₃ch", 6.75, 1.25, PP_ALIGN.CENTER),
            ("Stage assignment", 8.0, 4.75, PP_ALIGN.LEFT)]
    yh = 1.95
    rect(s, 0.55, yh, 12.2, 0.42, NAVY)
    for name, cx, cw, al in cols:
        textbox(s, cx, yh + 0.06, cw, 0.32,
                [[{"text": name, "size": 12.5, "color": WHITE, "bold": True, "font": HFONT}]],
                align=al)
    y = yh + 0.46
    rh = 0.405
    for i, (sheet, col, note) in enumerate(order):
        Lp, B, sig, w3, sflag, bflag = theta_line(sheet)
        if i % 2 == 1:
            rect(s, 0.55, y, 12.2, rh, PANEL)
        rect(s, 0.55, y, 0.09, rh, col)  # stage tab
        textbox(s, 0.6, y + 0.05, 3.0, 0.32, [[{"text": pretty(sheet), "size": 11.5, "color": INK, "bold": True}]])
        textbox(s, 3.6, y + 0.05, 0.95, 0.32, [[{"text": f"{Lp:.2f}", "size": 11.5, "color": INK}]], align=PP_ALIGN.CENTER)
        textbox(s, 4.55, y + 0.05, 1.15, 0.32, [[{"text": f"{B:.2f}{bflag}", "size": 11.5, "color": INK}]], align=PP_ALIGN.CENTER)
        textbox(s, 5.7, y + 0.05, 1.05, 0.32, [[{"text": f"{sig:.2f}{sflag}", "size": 11.5, "color": INK, "bold": bool(sflag)}]], align=PP_ALIGN.CENTER)
        textbox(s, 6.75, y + 0.05, 1.25, 0.32, [[{"text": f"{w3:,.0f}", "size": 11.5, "color": INK}]], align=PP_ALIGN.CENTER)
        textbox(s, 8.0, y + 0.05, 4.75, 0.32, [[{"text": note, "size": 11, "color": col, "bold": True}]])
        y += rh
    textbox(s, 0.6, y + 0.05, 12.0, 0.3,
            [[{"text": "Only MC3 SNaCl landed at an interior σ — that is the reference 'well-identified' case.",
               "size": 11.5, "color": GRAY, "italic": True}]])
    return s


# contour example slide ------------------------------------------------------
def s_contour(kicker, sheet, stage_col, verdict, headline, takeaway):
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, kicker, headline, pretty(sheet) + "   ·   σ × Lₚ, 5 channels")
    verdict_chip(s, verdict, 10.45, 0.36)
    h = add_image(s, CONTOUR_DIR / sheet / "objcontour-x_sigma-y_Lp.png", 0.55, 1.7, 12.2)
    yL = 1.7 + h + 0.06
    contour_legend(s, 0.55, yL, 12.2)
    yb = yL + 0.54
    Lp, B, sig, w3, sflag, bflag = theta_line(sheet)
    rect(s, 0.55, yb, 12.2, 7.05 - yb, BOXBG, line=stage_col, line_w=1, round_=True)
    textbox(s, 0.78, yb + 0.09, 1.3, 0.34, [[{"text": "θ here:", "size": 13, "color": NAVY, "bold": True, "font": HFONT}]])
    place_omath(s, theta_math(sheet), 1.95, yb + 0.06, 10.6, sz=14, color=INK,
                fallback=theta_fallback(sheet), align=PP_ALIGN.LEFT, h=0.4)
    textbox(s, 0.78, yb + 0.52, 12.0, 7.05 - yb - 0.56,
            [[{"text": takeaway, "size": 12.5, "color": INK}]])
    return s


# fit slide ------------------------------------------------------------------
def _fixed_plot(sheet, kind):
    hits = sorted((FITS_FIXED / sheet).glob(f"{kind}-*.png"))
    return hits[0] if hits else None


def s_fit(sheet, stage_col, stage_tag, verdict, note):
    s = slide(); bg(s, WHITE); logo(s)
    stage_chip(s, stage_tag, stage_col, x=9.55, y=0.34, w=3.25)
    title_block(s, "04 · CONTOUR-CONSISTENT FIT", pretty(sheet),
                "Mass & concentration per vial · DATA2-style · prediction shown only where permeate is collected")
    mp, cp = _fixed_plot(sheet, "mass"), _fixed_plot(sheet, "concentration")
    if mp:
        add_image(s, mp, 0.5, 1.95, 3.45)            # square (~3.3 h)
    if cp:
        add_image(s, cp, 4.15, 2.12, 6.25)           # wide 2:1 (~3.05 h); legend built in
    # right rail: theta + verdict
    rx = 10.6
    Lp, B, sig, w3, sflag, bflag = theta_line(sheet)
    rect(s, rx, 1.95, 2.25, 2.5, BOXBG, line=stage_col, line_w=1.25, round_=True)
    textbox(s, rx + 0.12, 2.04, 2.0, 0.3,
            [[{"text": "θ (contour-consistent)", "size": 11.5, "color": NAVY, "bold": True, "font": HFONT}]])
    sarr = r"\,\uparrow" if "↑" in sflag else (r"\,\downarrow" if "↓" in sflag else "")
    barr = r"\,\uparrow" if "↑" in bflag else ""
    place_omath(s, rf"L_p = {Lp:.2f}", rx + 0.14, 2.42, 2.0, sz=13, color=INK, fallback=f"Lp = {Lp:.2f}", h=0.32)
    place_omath(s, rf"B = {B:.2f}{barr}\ \mu\mathrm{{m/s}}", rx + 0.14, 2.84, 2.0, sz=12, color=INK, fallback=f"B = {B:.2f}{bflag} µm/s", h=0.3)
    place_omath(s, rf"\sigma = {sig:.2f}{sarr}", rx + 0.14, 3.24, 2.0, sz=13, color=INK, fallback=f"σ = {sig:.2f}{sflag}", h=0.32)
    textbox(s, rx + 0.14, 3.7, 2.0, 0.3, [[{"text": "WSSE₃ch", "size": 11, "color": GRAY}]])
    place_omath(s, rf"{w3:,.0f}".replace(",", "{,}"), rx + 0.14, 3.94, 2.0, sz=15, color=stage_col, fallback=f"{w3:,.0f}", h=0.34)
    verdict_chip(s, verdict, rx, 4.55, w=2.25)
    rect(s, rx, 5.08, 2.25, 1.55, PANEL, line=RGBColor(0xD7, 0xDC, 0xE1), line_w=0.75, round_=True)
    textbox(s, rx + 0.12, 5.16, 2.02, 1.45, [[{"text": note, "size": 11, "color": INK}]])
    # caption: startup-vial convention
    textbox(s, 0.5, 5.55, 9.8, 0.9,
            [[{"text": "Leading hold-up / startup vial excluded ", "size": 11.5, "color": NAVY, "bold": True},
              {"text": "(n", "size": 11.5, "color": NAVY, "bold": True},
              {"text": "v0", "size": 8, "color": NAVY, "bold": True},
              {"text": " = 2): the model integrates J", "size": 11.5, "color": INK},
              {"text": "w", "size": 8, "color": INK},
              {"text": " across the hold-up to carry the cell state, but mass is plotted only where permeate is collected in the vial — DATA1/DATA2 convention.", "size": 11.5, "color": INK}],
             [{"text": "Concentration legend is built into the panel (retentate / vial / permeate, measurement vs prediction).", "size": 11, "color": GRAY, "italic": True}]],
            space_after=2)
    return s


# section divider ------------------------------------------------------------
def s_divider(tag, title, sub, col):
    s = slide(); bg(s, NAVY)
    rect(s, 0.9, 2.55, 0.16, 1.45, col, round_=True)  # left accent bar (motif)
    textbox(s, 1.25, 2.55, 11.3, 0.4, [[{"text": tag, "size": 15, "color": col, "bold": True, "font": HFONT}]])
    textbox(s, 1.25, 2.9, 11.3, 0.9, [[{"text": title, "size": 34, "color": WHITE, "bold": True, "font": HFONT}]])
    textbox(s, 1.25, 3.62, 11.3, 0.6, [[{"text": sub, "size": 15, "color": RGBColor(0xBB, 0xC6, 0xD4)}]])
    return s


# all-11 small multiples -----------------------------------------------------
def s_smallmult():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "02 · IDENTIFIABILITY",
                "Every sheet's σ × Lₚ contour, colour-tabbed by stage",
                "Basin consistent (Stage 1) → split (Stage 2) → absent (Stage 3).  Each tile: 5 channels, contour = log₁₀ WSSE, ▲ = argmin.")
    # stage-colour key
    kx = 8.15
    for lab, col in [("Stage 1 well-id.", S1), ("Stage 2 inconsist.", S2), ("Stage 3 not-est.", S3)]:
        rect(s, kx, 1.55, 0.16, 0.16, col, round_=True)
        textbox(s, kx + 0.2, 1.47, 1.55, 0.3, [[{"text": lab, "size": 9.5, "color": col, "bold": True}]])
        kx += 1.55
    grid = [
        ("MC3.07.22.24_SNaCl", S1), ("MC5.07.23.24_S2NaCl", S1), ("MC4.07.11.24_SNaCl", S2),
        ("MC2.05.07.24_NaCl", S2), ("MC5.07.23.24_NaCl", S2), ("MC5.07.23.24_SNaCl", S2),
        ("MC2.05.07.24_CaCl2", S2), ("MC3.07.11.24_SCaCl2", S2), ("MC3.07.12.24_S2CaCl2", S2),
        ("MC2.05.21.24_LaCl3", S3), ("MC4.07.11.24_SLaCl3", S3),
    ]
    cols = 3; cw = 3.6; gx = 0.22; gy = 0.10
    x0 = 1.05; y0 = 1.82
    from PIL import Image
    for i, (sheet, col) in enumerate(grid):
        r, c = divmod(i, cols)
        x = x0 + c * (cw + gx)
        iw, ih = Image.open(CONTOUR_DIR / sheet / "objcontour-x_sigma-y_Lp.png").size
        h = cw * ih / iw
        rowh = h + 0.30
        y = y0 + r * (rowh + gy)
        rect(s, x - 0.03, y - 0.03, cw + 0.06, rowh + 0.06, WHITE, line=col, line_w=1.25, round_=True)
        rect(s, x - 0.03, y - 0.03, 0.09, rowh + 0.06, col)
        textbox(s, x + 0.12, y + 0.01, cw - 0.2, 0.26, [[{"text": pretty(sheet), "size": 10, "color": col, "bold": True}]])
        s.shapes.add_picture(str(CONTOUR_DIR / sheet / "objcontour-x_sigma-y_Lp.png"),
                             Inches(x + 0.05), Inches(y + 0.28), Inches(cw - 0.1), Inches((cw - 0.1) * ih / iw))
    return s


# animated identifiability ---------------------------------------------------
SCRUB = {
    "sigma": ("scrub-sigma_LpvsB.gif", "scrubbing σ  ·  basin shown in Lₚ × B"),
    "Lp":    ("scrub-Lp_Bvssigma.gif", "scrubbing Lₚ  ·  basin shown in B × σ"),
    "B":     ("scrub-B_Lpvssigma.gif", "scrubbing B  ·  basin shown in Lₚ × σ"),
}


def s_anim_intro():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "03 · IDENTIFIABILITY IN MOTION",
                "Watch the objective basin move",
                "Same 5-channel diagnostic, animated — hold one parameter fixed and sweep it frame by frame")
    rect(s, 0.6, 1.95, 12.13, 1.7, RGBColor(0xF8, 0xFA, 0xFC),
         line=RGBColor(0xD7, 0xDC, 0xE1), line_w=0.75, round_=True)
    textbox(s, 0.85, 2.08, 11.6, 0.35,
            [[{"text": "How to read the scrub animations", "size": 15, "color": NAVY, "bold": True, "font": HFONT}]])
    textbox(s, 0.85, 2.45, 11.7, 1.1,
            [[{"text": "Each frame fixes one parameter and shows the objective basin for the other two, across 3 channels ",
               "size": 13.5, "color": INK},
              {"text": "(mass · permeate-conc · retentate-conc)", "size": 12.5, "color": GRAY},
              {"text": ".  The red ▲ marks the per-frame minimum.", "size": 13.5, "color": INK}],
             [{"text": " ", "size": 6, "color": INK}],
             [{"text": "Read the motion as identifiability:  ", "size": 13.5, "color": NAVY, "bold": True},
              {"text": "a basin that stays put and the 3 channels agree → ", "size": 13, "color": INK},
              {"text": "identifiable", "size": 13, "color": S1, "bold": True},
              {"text": ".   A minimum that slides along a wall, or channels that disagree → ", "size": 13, "color": INK},
              {"text": "not", "size": 13, "color": S3, "bold": True},
              {"text": ".", "size": 13, "color": INK}]],
            space_after=2)
    # three verdict legend chips
    chips = [("GOOD — basin stays, channels agree", S1),
             ("OKAY — minimum slides to a wall", S2),
             ("POOR — no coherent basin", S3)]
    x = 0.6
    for lab, col in chips:
        rect(s, x, 4.0, 3.95, 0.55, col, round_=True)
        textbox(s, x + 0.15, 4.13, 3.7, 0.35, [[{"text": lab, "size": 12, "color": WHITE, "bold": True}]], align=PP_ALIGN.CENTER)
        x += 4.04
    textbox(s, 0.6, 4.95, 12.1, 0.9,
            [[{"text": "▶  These slides contain animated GIFs — they play automatically in PowerPoint slideshow mode. "
                       "The σ-sweep (B × Lₚ basin moving through σ) is the one to watch: it is the animated twin of the static σ × Lₚ panels.",
               "size": 12.5, "color": GRAY, "italic": True}]])
    return s


def s_anim(sheet, scrub, stage_col, stage_tag, verdict, headline, takeaway):
    s = slide(); bg(s, WHITE); logo(s)
    fname, sub = SCRUB[scrub]
    gif = ANIM_DIR / sheet / fname
    title_block(s, stage_tag, headline, pretty(sheet) + "   ·   " + sub)
    verdict_chip(s, verdict, 10.45, 0.36)
    if gif.exists():
        h = add_image(s, gif, 0.9, 1.78, 11.5)
        badge = "▶ animated — plays in slideshow"
    else:
        # fallback to the static 5-channel contour until the 3D grid finishes
        h = add_image(s, CONTOUR_DIR / sheet / "objcontour-x_sigma-y_Lp.png", 0.55, 1.78, 12.2)
        badge = "static σ × Lₚ (3D animation pending)"
    yb = 1.78 + h + 0.16
    Lp, B, sig, w3, sflag, bflag = theta_line(sheet)
    rect(s, 0.6, yb, 12.13, 7.05 - yb, BOXBG, line=stage_col, line_w=1, round_=True)
    textbox(s, 0.82, yb + 0.09, 11.8, 0.32,
            [[{"text": badge + "      ", "size": 11.5, "color": stage_col, "bold": True},
              {"text": f"warm-start θ:  Lₚ = {Lp:.2f}   B = {B:.2f} µm/s   σ = {sig:.2f}{sflag}   WSSE₃ch = {w3:,.0f}",
               "size": 12.5, "color": INK, "bold": True, "font": "Consolas"}]])
    textbox(s, 0.82, yb + 0.43, 11.8, 7.05 - yb - 0.45,
            [[{"text": takeaway, "size": 12.5, "color": INK}]])
    return s


# ---- narrative / math slides ------------------------------------------------
def s_thesis():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "THE GOAL",
                "Different experiments pin different parameters — use that",
                "Each run identifies some of θ = (Lₚ, B, σ) well and others not at all. Two questions follow:")
    # two big questions
    rect(s, 0.6, 1.95, 6.0, 1.5, RGBColor(0xEF, 0xF4, 0xEE), line=S1, line_w=1.25, round_=True)
    textbox(s, 0.8, 2.08, 5.6, 0.4, [[{"text": "Q1 · Combine", "size": 14, "color": S1, "bold": True, "font": HFONT}]])
    textbox(s, 0.8, 2.45, 5.6, 0.95,
            [[{"text": "Can we combine experiments into a better ", "size": 14, "color": INK},
              {"text": "initial guess", "size": 14, "color": INK, "bold": True},
              {"text": " for the joint fit — taking each parameter from the run that actually constrains it?", "size": 14, "color": INK}]])
    rect(s, 6.75, 1.95, 5.98, 1.5, RGBColor(0xFB, 0xF1, 0xE6), line=S2, line_w=1.25, round_=True)
    textbox(s, 6.95, 2.08, 5.6, 0.4, [[{"text": "Q2 · Design", "size": 14, "color": S2, "bold": True, "font": HFONT}]])
    textbox(s, 6.95, 2.45, 5.6, 0.95,
            [[{"text": "Which ", "size": 14, "color": INK},
              {"text": "new experiments", "size": 14, "color": INK, "bold": True},
              {"text": " would break the correlations that leave σ and B unidentified for the harder salts?", "size": 14, "color": INK}]])
    # logic flow
    textbox(s, 0.6, 3.75, 12.0, 0.35, [[{"text": "The logic flow of this talk", "size": 14, "color": NAVY, "bold": True, "font": HFONT}]])
    steps = [("1", "Model", "θ = (Lₚ, B, σ)\nin Spiegler–Kedem"),
             ("2", "Diagnose", "WSSE contour maps\n+ animations"),
             ("3", "Read where\nθ is pinned", "good / okay / poor\nper experiment"),
             ("4", "Combine →\nseed", "better initial\nguess per salt"),
             ("5", "Design new\nexperiments", "break σ–B,\npin Lₚ")]
    x = 0.6; w = 2.27; gap = 0.18
    for i, (n, head, body) in enumerate(steps):
        col = NAVY if i % 2 == 0 else GOLD
        rect(s, x, 4.2, w, 1.95, WHITE, line=col, line_w=1.25, round_=True)
        rect(s, x, 4.2, w, 0.5, col, round_=True)
        textbox(s, x, 4.27, w, 0.4, [[{"text": "STEP " + n, "size": 11.5, "color": WHITE, "bold": True, "font": HFONT}]], align=PP_ALIGN.CENTER)
        textbox(s, x + 0.12, 4.78, w - 0.24, 0.7, [[{"text": head, "size": 13.5, "color": NAVY, "bold": True, "font": HFONT}]], align=PP_ALIGN.CENTER)
        textbox(s, x + 0.12, 5.5, w - 0.24, 0.6, [[{"text": body, "size": 11, "color": GRAY}]], align=PP_ALIGN.CENTER)
        if i < len(steps) - 1:
            textbox(s, x + w - 0.02, 4.95, 0.22, 0.4, [[{"text": "→", "size": 20, "color": col, "bold": True}]], align=PP_ALIGN.CENTER)
        x += w + gap
    textbox(s, 0.6, 6.4, 12.1, 0.5,
            [[{"text": "Takeaway up front:  ", "size": 13.5, "color": S1, "bold": True},
              {"text": "the campaign already contains everything we need to seed every salt's fit — and tells us exactly which two experiments would close the gaps.",
               "size": 13.5, "color": INK}]])
    return s


def s_model_math():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "01 · THE MODEL",
                "What we simulate and what we estimate",
                "Spiegler–Kedem transport through the NF membrane — three parameters per experiment: θ = (Lₚ, B, σ)")
    # paper notation (Ouimet et al., J. Membrane Sci. 641, 2022 — DATA1; DATA2)
    eqs = [
        ("Water flux  ·  Eq. (1)",
         r"J_w = L_p\left(\Delta P - \sigma\,\Delta\pi\right)",
         "Jw = Lp (ΔP − σ Δπ)", 22,
         "σ scales how much osmotic back-pressure opposes the applied ΔP"),
        ("Osmotic pressure  ·  van 't Hoff",
         r"\Delta\pi = n R T\left(c_{in} - c_h\right)",
         "Δπ = n R T (c_in − c_h)", 22,
         "n = dissolved species per formula unit (NaCl 2 · CaCl₂ 3 · LaCl₃ 4)"),
        ("Concentration polarization  ·  Eq. (3)",
         r"c_{in} = \left(c_f - c_h\right)\exp\!\left(\frac{J_w}{k}\right) + c_h",
         "c_in = (c_f − c_h) · exp(Jw / k) + c_h", 20,
         "wall enrichment from the boundary layer (mass-transfer coefficient k)"),
        ("Solute flux  ·  Eq. (2)",
         r"J_s = B\left(c_{in} - c_h\right)",
         "Js = B (c_in − c_h)", 22,
         "diffusion-dominant; B = solute permeability coefficient"),
    ]
    y = 1.85
    for name, latex, fb, sz, note in eqs:
        rect(s, 0.6, y, 12.13, 1.0, RGBColor(0xF8, 0xFA, 0xFC), line=RGBColor(0xDD, 0xE2, 0xE7), line_w=0.75, round_=True)
        textbox(s, 0.8, y + 0.08, 4.7, 0.35, [[{"text": name, "size": 13, "color": NAVY, "bold": True, "font": HFONT}]])
        textbox(s, 0.8, y + 0.48, 5.2, 0.45, [[{"text": note, "size": 11.5, "color": GRAY, "italic": True}]])
        place_omath(s, latex, 6.1, y + 0.16, 6.4, sz=sz, color=INK, fallback=fb, align=PP_ALIGN.CENTER, h=0.66)
        y += 1.1
    rect(s, 0.6, 6.42, 12.13, 0.66, BOXBG, line=GOLD, line_w=1, round_=True)
    textbox(s, 0.8, 6.5, 12.0, 0.55,
            [[{"text": "Same notation as the DATA1 / DATA2 papers. ", "size": 12.5, "color": NAVY, "bold": True},
              {"text": "Concentration fits curve because the state ODEs are nonlinear; DATA3 also allows B = B(c_f) ", "size": 12.5, "color": INK},
              {"text": "(B_form = 1)", "size": 12.5, "color": INK, "bold": True, "font": "Consolas"},
              {"text": " when a constant B underfits.", "size": 12.5, "color": INK}]])
    return s


def s_diag_math():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "02 · THE DIAGNOSTIC",
                "Contour maps = forward-simulated WSSE, no fitting",
                "We sweep a parameter grid, simulate the full DAE at each node, and score the weighted residual per channel")
    rect(s, 0.6, 1.95, 12.13, 1.05, RGBColor(0xF8, 0xFA, 0xFC), line=RGBColor(0xDD, 0xE2, 0xE7), line_w=0.75, round_=True)
    textbox(s, 0.8, 2.04, 11.0, 0.35, [[{"text": "Weighted sum of squared residuals, summed over the 3 fitted channels (mass · permeate cᵥ · retentate cᵣ)", "size": 12.5, "color": GRAY, "italic": True}]])
    place_omath(s, r"\mathrm{WSSE}(\theta)=\sum_{ch\,\in\,\{m,c_v,c_r\}}\frac{1}{N_{ch}}\sum_i\left(y_i^{\,model}(\theta)-y_i^{\,exp}\right)^2",
                0.8, 2.4, 11.8, sz=17, color=INK,
                fallback="WSSE(θ) = Σ ch∈{m,cv,cr}  (1/N_ch) Σ_i ( y_i model − y_i exp )²",
                align=PP_ALIGN.CENTER, h=0.6)
    rect(s, 0.6, 3.2, 5.95, 1.7, WHITE, line=NAVY, line_w=1.25, round_=True)
    textbox(s, 0.8, 3.3, 5.6, 0.35, [[{"text": "Best fit", "size": 13, "color": NAVY, "bold": True, "font": HFONT}]])
    place_omath(s, r"\theta^{*} = \arg\min_{\theta}\,\mathrm{WSSE}(\theta)", 0.8, 3.74, 5.55, sz=19, color=INK,
                fallback="θ* = argmin_θ  WSSE(θ)", align=PP_ALIGN.CENTER, h=0.5)
    textbox(s, 0.8, 4.35, 5.6, 0.5, [[{"text": "θ = (Lₚ, B, σ).  If the minimum is a sharp interior point, θ is identified; if it lies on a wall or in a flat valley, it is not.", "size": 11.5, "color": GRAY}]])
    rect(s, 6.78, 3.2, 5.95, 1.7, WHITE, line=GOLD, line_w=1.25, round_=True)
    textbox(s, 6.98, 3.3, 5.6, 0.35, [[{"text": "Contour-consistent seed", "size": 13, "color": NAVY, "bold": True, "font": HFONT}]])
    place_omath(s, r"\theta_{init}=\arg\min_{(\sigma,\,L_p)}\left[\mathrm{WSSE}_m+\mathrm{WSSE}_{c_v}+\mathrm{WSSE}_{c_r}\right]",
                6.98, 3.78, 5.55, sz=14, color=INK,
                fallback="θ_init = argmin (σ,Lp) [ WSSE_m + WSSE_cv + WSSE_cr ]", align=PP_ALIGN.CENTER, h=0.5)
    textbox(s, 6.98, 4.35, 5.6, 0.5, [[{"text": "Read the joint 2-D minimum straight off the contour (B, S₀, S held at warm-start) and use it to seed the full fit.", "size": 11.5, "color": GRAY}]])
    rect(s, 0.6, 5.15, 12.13, 1.9, RGBColor(0xF4, 0xF6, 0xF8), line=RGBColor(0xDD, 0xE2, 0xE7), line_w=0.75, round_=True)
    textbox(s, 0.8, 5.25, 12.0, 0.35, [[{"text": "Two views of the same objective, both shown in this deck", "size": 13, "color": NAVY, "bold": True, "font": HFONT}]])
    textbox(s, 0.8, 5.62, 11.9, 1.3,
            [[{"text": "•  Static σ × Lₚ panels — ", "size": 12.5, "color": INK, "bold": True},
              {"text": "all 5 measurement channels at once (mass, perm/ret conc, perm/ret conductivity). Where channels agree → trust.", "size": 12.5, "color": INK}],
             [{"text": "•  Animated scrubs — ", "size": 12.5, "color": INK, "bold": True},
              {"text": "hold one parameter and sweep it; watch whether the 2-D basin for the other two stays put. This is how to read a good initial guess for each experiment.", "size": 12.5, "color": INK}],
             [{"text": "   ① B × Lₚ basin as σ : 0 → 1     ② σ × Lₚ basin as B : 0 → B_max", "size": 12, "color": GRAY, "font": "Consolas"}]],
            space_after=4)
    return s


def s_identif_math():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "02 · WHY EXPERIMENTS DIFFER",
                "σ leaves a fingerprint only when osmotic pressure is large",
                "The sensitivity of the data to σ is set by Δπ — so the operating regime decides what each run can identify")
    rect(s, 0.6, 1.9, 12.13, 1.0, BOXBG, line=GOLD, line_w=1, round_=True)
    textbox(s, 0.8, 1.98, 4.5, 0.5, [[{"text": "σ enters the model only through σ·Δπ, so", "size": 12.5, "color": GRAY, "italic": True}]])
    place_omath(s, r"\frac{\partial J_w}{\partial \sigma}=-L_p\,\Delta\pi=-L_p\,n R T\left(c_{in}-c_h\right)",
                5.1, 2.0, 7.5, sz=21, color=NAVY,
                fallback="∂Jw/∂σ = −Lp Δπ = −Lp · n R T (c_in − c_h)", align=PP_ALIGN.CENTER, h=0.8)
    cases = [
        (S1, "High c_f  →  σ visible", r"\Delta\pi = n R T\,c_f\quad(c_h\approx 0)",
         "Δπ = n R T cf   (c_h ≈ 0)",
         "Large cf → large Δπ → ∂Jw/∂σ large → σ leaves a fingerprint (DATA1 §4.3.2). Both gold sheets are diluting runs starting at cf ≈ 95 mM."),
        (NAVY, "Low c_f  →  Lₚ only", r"c_f\to 0\;\Rightarrow\;J_w = L_p\,\Delta P",
         "cf → 0  ⟹  Jw = Lp ΔP",
         "Δπ → 0: Lₚ reads directly from mass, independent of σ and B. This is a pure-water run, or the end of any diluting run."),
        (S2, "σ–B coupling", r"\frac{dc_f}{dt}=\frac{A_m}{m_f}\left(J_w c_d - J_s\right)",
         "dc_f/dt = (A_m / m_f)(Jw c_d − Js)",
         "σ acts only through Jw, B only through Js — but both drive c_f(t). Vary ΔP (hence Jw) and their relative pull changes, separating them."),
    ]
    x = 0.6; w = 3.94
    for col, head, latex, fb, note in cases:
        rect(s, x, 3.1, w, 3.55, WHITE, line=col, line_w=1.25, round_=True)
        rect(s, x, 3.1, w, 0.55, col, round_=True)
        textbox(s, x + 0.15, 3.18, w - 0.3, 0.4, [[{"text": head, "size": 14, "color": WHITE, "bold": True, "font": HFONT}]], align=PP_ALIGN.CENTER)
        rect(s, x + 0.25, 3.9, w - 0.5, 0.85, RGBColor(0xF6, 0xF8, 0xFA), round_=True)
        place_omath(s, latex, x + 0.35, 4.0, w - 0.7, sz=17, color=NAVY, fallback=fb, align=PP_ALIGN.CENTER, h=0.66)
        textbox(s, x + 0.25, 4.92, w - 0.5, 1.65, [[{"text": note, "size": 12.5, "color": INK}]])
        x += w + 0.18
    textbox(s, 0.6, 6.78, 12.1, 0.4,
            [[{"text": "So: ", "size": 13, "color": NAVY, "bold": True},
              {"text": "a diluting run from high cf pins Lₚ AND σ; a concentrating run from low cf pins only Lₚ; a pressure sweep separates σ from B.",
               "size": 13, "color": INK}]])
    return s


def s_anim2(sheet, stage_col, stage_tag, verdict, headline, readoff):
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, stage_tag, headline, pretty(sheet) + "   ·   read the joint minimum off both sweeps")
    verdict_chip(s, verdict, 10.45, 0.36)
    gw = 7.3
    g1 = ANIM_DIR / sheet / "scrub-sigma_LpvsB.gif"
    g2 = ANIM_DIR / sheet / "scrub-B_Lpvssigma.gif"
    textbox(s, 0.5, 1.46, 7.3, 0.26, [[{"text": "①  Lₚ × B  basin   as   σ : 0 → 1   (Lₚ on Y)", "size": 12, "color": stage_col, "bold": True, "font": HFONT}]])
    if g1.exists():
        add_image(s, g1, 0.5, 1.7, gw)
    else:
        add_image(s, CONTOUR_DIR / sheet / "objcontour-x_sigma-y_Lp.png", 0.5, 1.72, gw)
    textbox(s, 0.5, 4.02, 7.3, 0.26, [[{"text": "②  σ × Lₚ  basin   as   B : 0 → B_max", "size": 12, "color": stage_col, "bold": True, "font": HFONT}]])
    if g2.exists():
        add_image(s, g2, 0.5, 4.26, gw)
    else:
        textbox(s, 0.5, 4.8, 7.2, 0.5, [[{"text": "(σ × Lₚ sweep renders once the 3-D grid finishes)", "size": 12, "color": GRAY, "italic": True}]])
    contour_legend(s, 0.5, 6.62, 7.3, h=0.42)
    # right rail: read-off seed (native math)
    rx = 8.45; rw = 4.35
    Lp, B, sig, w3, sflag, bflag = theta_line(sheet)
    rect(s, rx, 1.7, rw, 2.62, RGBColor(0xEF, 0xF4, 0xEE) if verdict == "GOOD" else PANEL,
         line=stage_col, line_w=1.5, round_=True)
    textbox(s, rx + 0.18, 1.8, rw - 0.36, 0.32, [[{"text": "Initial guess to read off", "size": 13.5, "color": NAVY, "bold": True, "font": HFONT}]])
    place_omath(s, rf"L_p \approx {Lp:.1f}", rx + 0.22, 2.2, rw - 0.5, sz=17, color=INK, fallback=f"Lp ≈ {Lp:.1f}", h=0.4)
    place_omath(s, rf"B \approx {B:.1f}\ \mu\mathrm{{m/s}}", rx + 0.22, 2.68, rw - 0.5, sz=17, color=INK, fallback=f"B ≈ {B:.1f} µm/s", h=0.4)
    place_omath(s, rf"\sigma \approx {sig:.2f}", rx + 0.22, 3.16, rw - 0.5, sz=17, color=INK, fallback=f"σ ≈ {sig:.2f}", h=0.4)
    textbox(s, rx + 0.18, 3.78, rw - 0.36, 0.5, [[{"text": "joint argmin of mass + cᵥ + cᵣ", "size": 10.5, "color": GRAY, "italic": True}]])
    rect(s, rx, 4.46, rw, 2.42, BOXBG, line=stage_col, line_w=1, round_=True)
    textbox(s, rx + 0.18, 4.56, rw - 0.36, 2.2, [[{"text": readoff, "size": 12.5, "color": INK}]])
    return s


def s_synthesis():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "Q1 · COMBINE",
                "Take each parameter from the run that pins it",
                "Same θ = (Lₚ, B, σ), but each component seeded from the experiment whose data actually constrain it")
    rows = [
        (S1, "Lₚ", r"J_w = L_p\,\Delta P", "Jw = Lp ΔP", "Pure-water run, or any diluting run's mass channel",
         "Δπ ≈ 0 → Lₚ falls straight out of the mass data, untouched by σ or B."),
        (GOLD, "σ", r"\arg\min_{\sigma}\,\mathrm{WSSE}", "argmin_σ WSSE", "A high-cf (diluting) run — MC3 SNaCl, σ ≈ 0.79",
         "Only where Δπ/ΔP is appreciable does σ leave an interior minimum."),
        (NAVY, "B", r"J_s = B\left(c_{in}-c_h\right)", "Js = B (c_in − c_h)", "The low-residual solute channel where Jw spans a range",
         "Tighten with the per-salt bound (halve per +1 cation charge)."),
    ]
    y = 1.98
    for col, sym, latex, fb, src, note in rows:
        rect(s, 0.6, y, 12.13, 1.32, WHITE, line=col, line_w=1.25, round_=True)
        rect(s, 0.6, y, 1.0, 1.32, col, round_=True)
        textbox(s, 0.6, y + 0.42, 1.0, 0.5, [[{"text": sym, "size": 26, "color": WHITE, "bold": True, "font": HFONT}]], align=PP_ALIGN.CENTER)
        place_omath(s, latex, 1.8, y + 0.22, 4.3, sz=18, color=INK, fallback=fb, align=PP_ALIGN.LEFT, h=0.5)
        textbox(s, 1.8, y + 0.8, 4.3, 0.45, [[{"text": "from:  " + src, "size": 11.5, "color": NAVY, "bold": True}]])
        textbox(s, 6.4, y + 0.3, 6.1, 0.8, [[{"text": note, "size": 12.5, "color": INK}]])
        y += 1.42
    rect(s, 0.6, 6.32, 12.13, 0.8, NAVY, round_=True)
    textbox(s, 0.85, 6.4, 12.0, 0.7,
            [[{"text": "θ_init per salt = ( Lₚ : water / diluting   ·   σ : σ-probing run   ·   B : clean solute channel )",
               "size": 12.5, "color": WHITE, "bold": True, "font": "Consolas"}],
             [{"text": "Transfer across salts of the same valence, then run the joint multistart from this seed instead of a blind LHS.",
               "size": 11.5, "color": RGBColor(0xCC, 0xD6, 0xE2)}]], space_after=2)
    return s


# appendix B contours --------------------------------------------------------
def s_appendix_B():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "APPENDIX",
                "MC2 NaCl — B-direction contours behind the σ=1 wall",
                "The σ=1 wall on the Stage-2 NaCl sheets hides a B–Lₚ trade-off")
    h1 = add_image(s, ALLPAIRS / "objcontour-x_B-y_Lp.png", 0.55, 1.95, 6.0)
    # axis-convention: prefer the B-on-Y file (objcontour-x_sigma-y_B.png); fall
    # back to the legacy sigma-on-Y file until the regen lands.
    bsig_new = ALLPAIRS / "objcontour-x_sigma-y_B.png"
    bsig = bsig_new if bsig_new.exists() else (ALLPAIRS / "objcontour-x_B-y_sigma.png")
    add_image(s, bsig, 6.75, 1.95, 6.0)
    bsig_lbl = ("σ × B  (Lₚ fixed) — B on Y per axis convention"
                if bsig_new.exists() else "B × σ  (Lₚ fixed) — mass argmins at interior σ ≈ 0.32")
    textbox(s, 0.55, 1.95 + h1 + 0.1, 6.0, 0.4, [[{"text": "B × Lₚ  (σ fixed at warm-start = 1) — Lₚ on Y", "size": 12, "color": GRAY, "italic": True}]], align=PP_ALIGN.CENTER)
    textbox(s, 6.75, 1.95 + h1 + 0.1, 6.0, 0.4, [[{"text": bsig_lbl, "size": 12, "color": GRAY, "italic": True}]], align=PP_ALIGN.CENTER)
    return s


# closing: suggested experiments (math-tagged) -------------------------------
def s_data():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "Q2 · DESIGN",
                "New experiments, chosen to fix exactly what the math says is missing",
                "Each one targets a specific identifiability gap from the sensitivity analysis")
    items = [
        ("PRIORITY 1", S3, "Pure-water Lₚ before & after every run", "→ pins Lₚ",
         r"J_w=L_p\,\Delta P", "Jw = Lp ΔP",
         "Δπ = 0 isolates Lₚ so σ and B can't compensate for it. Doubles as a fouling check (pre/post drift)."),
        ("PRIORITY 2", S2, "Multi-ΔP sweep within a run  +  dilute runs", "→ separates σ–B",
         r"\frac{\partial J_w}{\partial\sigma}=-L_p\,\Delta\pi", "∂Jw/∂σ = −Lp Δπ",
         "Stepping ΔP changes Jw, shifting how σ (via Jw) and B (via Js) each act on c_f — breaking their correlation. Dilute runs raise Δπ/ΔP so σ identifies for CaCl₂/LaCl₃."),
        ("PRIORITY 3", NAVY, "More LaCl₃ + retentate conductivity + ICP speciation", "→ tests model form",
         r"n_{eff}=?\;\;\left(\mathrm{La}^{3+}\rightleftharpoons\mathrm{LaCl}^{2+}\right)", "n_eff = ?  (La³⁺ ⇌ LaCl²⁺)",
         "n = 2 is an anecdote (need ≥5). Retentate-side conductivity kills the CaCl₂ side-split; speciation tests the ion-pairing hypothesis directly."),
    ]
    y = 1.98
    for tag, col, head, tagline, latex, fb, body in items:
        rect(s, 0.6, y, 12.13, 1.5, WHITE, line=col, line_w=1.25, round_=True)
        rect(s, 0.6, y, 1.85, 1.5, col, round_=True)
        rect(s, 1.1, y, 1.35, 1.5, col)
        textbox(s, 0.65, y + 0.32, 1.75, 0.4, [[{"text": tag, "size": 13, "color": WHITE, "bold": True, "font": HFONT}]], align=PP_ALIGN.CENTER)
        textbox(s, 0.65, y + 0.78, 1.75, 0.4, [[{"text": tagline, "size": 12, "color": WHITE, "bold": True}]], align=PP_ALIGN.CENTER)
        textbox(s, 2.65, y + 0.14, 9.9, 0.4, [[{"text": head, "size": 15, "color": NAVY, "bold": True, "font": HFONT}]])
        place_omath(s, latex, 2.65, y + 0.56, 4.2, sz=16, color=col, fallback=fb, align=PP_ALIGN.LEFT, h=0.42)
        textbox(s, 2.65, y + 1.08, 9.9, 0.4, [[{"text": body, "size": 11.5, "color": GRAY}]])
        y += 1.62
    return s


# closing: questions ---------------------------------------------------------
def s_questions():
    s = slide(); bg(s, NAVY)
    textbox(s, 0.9, 0.7, 11.5, 0.4, [[{"text": "QUESTIONS FOR DISCUSSION", "size": 15, "color": GOLD, "bold": True, "font": HFONT}]])
    textbox(s, 0.9, 1.1, 11.5, 0.7, [[{"text": "What I'd like your read on before the re-run", "size": 26, "color": WHITE, "bold": True, "font": HFONT}]])
    qs = [
        ("1", "Does the 3-stage diagnosis match your physical intuition? Increasing cation valence → identifiability degrades. Is the σ=1 wall a regime problem, or the cᶠ weighting forcing it there?"),
        ("2", "Two paths to well-identified — MC3 SNaCl (interior σ = 0.79) vs MC5 S2NaCl (diluting, σ irrelevant). What in the run-sheet metadata differs from the other NaCl sheets?"),
        ("3", "MC4 SNaCl is diluting but fits ~10× worse (WSSE ≈ 19.5k). Likely a data-quality issue (fouling / leak / calibration) rather than model — worth a closer look?"),
        ("4", "LaCl₃ fix: parameterize σ(cᶠ), add an ion-pairing equilibrium block, or fit only the diluting portion? Either way we need 3–5 more sheets."),
        ("5", "Per-salt B bounds — NaCl [0,30], CaCl₂ [0,15], LaCl₃ [0,10]. Anything you'd revise before re-running the full campaign?"),
    ]
    y = 2.1
    for n, q in qs:
        rect(s, 0.9, y + 0.02, 0.5, 0.5, GOLD, round_=True)
        textbox(s, 0.9, y + 0.08, 0.5, 0.4, [[{"text": n, "size": 16, "color": NAVY, "bold": True, "font": HFONT}]], align=PP_ALIGN.CENTER)
        textbox(s, 1.6, y, 11.0, 0.9, [[{"text": q, "size": 13.5, "color": RGBColor(0xE6, 0xEC, 0xF3)}]])
        y += 1.0
    return s


# ---- follow-up results slides ----------------------------------------------
def s_meas_error():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "FOLLOW-UP · MEASUREMENT ERROR",
                "The cF weight was physically unjustified — and σ depends on it",
                "Paper: Lilonfe et al., 'Soft Sensors…', ChemRxiv 2026 (doi:10.26434/chemrxiv.15002541/v1)")
    # left: the weighting fix
    rect(s, 0.6, 1.95, 6.0, 2.3, PANEL, line=RGBColor(0xD7, 0xDC, 0xE1), line_w=0.75, round_=True)
    textbox(s, 0.8, 2.05, 5.6, 0.35, [[{"text": "Per-channel measurement error (DATA3)", "size": 13.5, "color": NAVY, "bold": True, "font": HFONT}]])
    rows = [("Retentate cF", "0.3% → 2%", "conductivity soft-sensor MAPE ⊕ cell-const", S2),
            ("Permeate cV", "3% (kept)", "ICP-OES reference — already physical", NAVY),
            ("Mass m", "0.01 g (kept)", "balance resolution", GRAY)]
    yy = 2.5
    for ch, val, why, col in rows:
        textbox(s, 0.8, yy, 1.7, 0.3, [[{"text": ch, "size": 12.5, "color": INK, "bold": True}]])
        textbox(s, 2.5, yy, 1.5, 0.3, [[{"text": val, "size": 12.5, "color": col, "bold": True}]])
        textbox(s, 0.8, yy + 0.3, 5.6, 0.3, [[{"text": why, "size": 10.5, "color": GRAY, "italic": True}]])
        yy += 0.6
    # right: the sigma-fragility finding
    rect(s, 6.75, 1.95, 5.98, 2.3, BOXBG, line=S3, line_w=1.25, round_=True)
    textbox(s, 6.95, 2.05, 5.6, 0.35, [[{"text": "σ flips with the weight (multistart-verified)", "size": 13.5, "color": NAVY, "bold": True, "font": HFONT}]])
    textbox(s, 6.95, 2.45, 5.7, 1.7,
            [[{"text": "MC3 SNaCl:  ", "size": 12.5, "color": INK, "bold": True},
              {"text": "σ 0.42 (interior, 0.3%)  →  1.0 (wall, 2%)", "size": 12.5, "color": INK}],
             [{"text": "MC5 S2NaCl: ", "size": 12.5, "color": INK, "bold": True},
              {"text": "σ 0.0 (wall, 0.3%)  →  0.35 (interior, 2%)", "size": 12.5, "color": INK}],
             [{"text": " ", "size": 6, "color": INK}],
             [{"text": "σ's global minimizer is weight-dependent → σ is poorly identified.",
               "size": 12, "color": S3, "bold": True}]], space_after=3)
    # bottom takeaway
    rect(s, 0.6, 4.45, 12.13, 2.1, NAVY, round_=True)
    textbox(s, 0.85, 4.6, 12.0, 0.4, [[{"text": "Takeaway", "size": 14, "color": GOLD, "bold": True, "font": HFONT}]])
    textbox(s, 0.85, 5.0, 12.0, 1.5,
            [[{"text": "The interior σ we reported was partly an ", "size": 14, "color": WHITE},
              {"text": "artifact of the over-tight 0.3% cF weight", "size": 14, "color": WHITE, "bold": True},
              {"text": ".  With the honest 2% error, σ is revealed to be ", "size": 14, "color": WHITE},
              {"text": "not robustly identifiable", "size": 14, "color": GOLD, "bold": True},
              {"text": " from single-salt runs.", "size": 14, "color": WHITE}],
             [{"text": " ", "size": 6, "color": WHITE}],
             [{"text": "→  This is a sensitivity result, not a fit fix.  σ needs better EXPERIMENTS (high-cF, multi-ΔP DoE), not reweighting. "
                       "The weight now feeds objective + Fisher-information + parmest consistently; DATA1/DATA2 byte-identical.",
               "size": 12.5, "color": RGBColor(0xCC, 0xD6, 0xE2)}]], space_after=3)
    return s


def s_directions():
    s = slide(); bg(s, WHITE); logo(s)
    title_block(s, "FOLLOW-UP · DIRECTIONS",
                "Where this points: B(ionic strength) and per-band θ",
                "Concentration-dependent transport, made salt-transferable and band-resolved")
    cards = [
        (NAVY, "β(c) → β(I)", [
            ("B = J_w[β₀ + β₁·c + β₂·c²]", True),
            ("plot β₀/β₁ vs Lₚ or σ (not lumped B)", False),
            ("I = k_I·c   (NaCl 1 · CaCl₂ 3 · LaCl₃ 6)", True),
            ("B(I) ≡ B(c) per salt (β rescaled), but β become transferable across salts → multisalt", False)]),
        (S2, "Per-band / per-vial θ", [
            ("group vials by cF band → one θ per band", False),
            ("mask the objective to a band's vials", False),
            ("B_form='pervial' already gives per-vial B", True),
            ("correlate band θ with contour basins → σ-drift with concentration", False)]),
        (S1, "Axis convention", [
            ("Y-priority  Lₚ > B > σ", True),
            ("animations + contours regenerated", False),
            ("σ-scrub now Lₚ × B (Lₚ on Y)", False),
            ("rule of thumb going forward", False)]),
    ]
    x = 0.6; w = 3.94; gap = 0.18
    for col, head, items in cards:
        rect(s, x, 1.95, w, 4.55, WHITE, line=col, line_w=1.25, round_=True)
        rect(s, x, 1.95, w, 0.62, col, round_=True)
        textbox(s, x + 0.18, 2.06, w - 0.36, 0.4, [[{"text": head, "size": 15, "color": WHITE, "bold": True, "font": HFONT}]], align=PP_ALIGN.CENTER)
        yy = 2.82
        for txt, mono in items:
            rect(s, x + 0.2, yy + 0.05, 0.1, 0.1, col)
            textbox(s, x + 0.42, yy - 0.02, w - 0.62, 0.85,
                    [[{"text": txt, "size": 11.5, "color": INK, "font": ("Consolas" if mono else BFONT)}]])
            yy += 0.9
        x += w + gap
    return s


# ===========================================================================
# assemble order
# ===========================================================================
builders = []
builders.append(s_title)
# ---- post-meeting follow-up headlines (new) ----
builders.append(s_meas_error)
builders.append(s_directions)
builders.append(s_thesis)
# section 01 — the model
builders.append(s_model_math)
builders.append(s_bounds)
# section 02 — diagnose
builders.append(lambda: s_divider("02 · DIAGNOSE", "Map every parameter, then read what's identified",
                                   "Forward-simulated WSSE contours + animations — good · okay · poor, per experiment", GOLD))
builders.append(s_diag_math)
builders.append(s_framework)
builders.append(s_table)
builders.append(lambda: s_contour("STAGE 1 · WELL-IDENTIFIED", "MC3.07.22.24_SNaCl", S1, "GOOD",
    "GOOD — interior σ from a high-cf diluting run (gold standard)",
    "A diluting NaCl run starting at cf ≈ 95 mM, so Δπ is large early and σ leaves a real fingerprint (DATA1 §4.3.2). All 5 channels at the warm-start B; mass and permeate-concentration both argmin at an interior σ near Lₚ ≈ 8.6. This is the only NaCl sheet where σ is genuinely identified, and the cleanest fit in the campaign."))
builders.append(s_smallmult)
builders.append(s_identif_math)
# section 03 — in motion (both scrubs, read the initial guess)
builders.append(lambda: s_divider("03 · IN MOTION", "Read the initial guess off the moving basin",
                                   "Two sweeps per experiment: ① B×Lₚ as σ:0→1   ② σ×Lₚ as B:0→Bₘₐₓ", GOLD))
builders.append(s_anim_intro)
builders.append(lambda: s_anim2("MC3.07.22.24_SNaCl", S1, "STAGE 1 · WELL-IDENTIFIED", "GOOD",
    "GOOD — the basin holds under both sweeps",
    "Across both sweeps the minimum stays anchored (Lₚ ≈ 8.6, interior σ) and the channels track together. A sharp, stationary basin = a trustworthy seed. This is the gold-standard read: take all three parameters from this run."))
builders.append(lambda: s_anim2("MC2.05.07.24_NaCl", S2, "STAGE 2 · INCONSISTENT (NaCl)", "OKAY",
    "OKAY — Lₚ holds, σ slides to the wall",
    "Sweep ① : the B×Lₚ minimum keeps the same Lₚ but drifts to the σ=1 wall. Sweep ② : as B varies the σ×Lₚ minimum barely moves in Lₚ — Lₚ is solid. Take Lₚ from here; do NOT trust σ — get σ from the gold-standard run instead."))
builders.append(lambda: s_anim2("MC2.05.07.24_CaCl2", S2, "STAGE 2 · INCONSISTENT (CaCl₂)", "OKAY",
    "OKAY — Lₚ usable, σ and B trade off",
    "The two sweeps show the σ–B trade-off directly: moving B shifts where σ wants to sit. Lₚ is still the reliable read; σ and B are correlated and need a pressure sweep to separate. Seed Lₚ from here, σ from the σ-probing run, B from its bound."))
builders.append(lambda: s_anim2("MC2.05.21.24_LaCl3", S3, "STAGE 3 · NOT ESTIMABLE", "POOR",
    "POOR — no stationary basin under either sweep",
    "Neither sweep produces a consistent minimum and the channels never agree — the constant-(σ,B) model can't represent La³⁺. Don't read a seed from this run; it needs a model fix (ion-pairing) and more replicates before any θ is meaningful."))
# section 04 — fits
builders.append(lambda: s_divider("04 · FITS", "Mass & concentration fits from the contour seed",
                                   "θ from the joint contour argmin  ·  variable B (B_form = 1) → curved cₕ(t)", S1))
NOTE_B = "Variable B (B_form = 1): the concentration curve bends — a constant-B model would draw a straight line here."
fit_order = [
    ("MC3.07.22.24_SNaCl", S1, "STAGE 1A · well-identified", "GOOD",
     "Gold standard. Interior σ → both mass and cᶠ tracked. " + NOTE_B),
    ("MC5.07.23.24_S2NaCl", S1, "STAGE 1B · diluting", "GOOD",
     "Diluting run → cᶠ falls. Lₚ pinned by mass; σ irrelevant. " + NOTE_B),
    ("MC4.07.11.24_SNaCl", S2, "STAGE 1B · diluting ⚠", "OKAY",
     "Diluting but WSSE ~10× worse — suspected data-quality issue. " + NOTE_B),
    ("MC2.05.07.24_NaCl", S2, "STAGE 2 · NaCl, σ=1 wall", "OKAY",
     "Mass fit good; concentration biased by σ pinned at 1. " + NOTE_B),
    ("MC5.07.23.24_NaCl", S2, "STAGE 2 · NaCl, σ=1 wall", "OKAY",
     "Concentrating; σ at wall. " + NOTE_B),
    ("MC5.07.23.24_SNaCl", S2, "STAGE 2 · NaCl, σ=1 wall", "OKAY",
     "Concentrating; σ at wall. " + NOTE_B),
    ("MC2.05.07.24_CaCl2", S2, "STAGE 2 · CaCl₂", "OKAY",
     "Channel side-split; high retentate-conc WSSE. " + NOTE_B),
    ("MC3.07.11.24_SCaCl2", S2, "STAGE 2 · CaCl₂", "OKAY",
     "σ at wall. " + NOTE_B),
    ("MC3.07.12.24_S2CaCl2", S2, "STAGE 2 · CaCl₂", "OKAY",
     "σ AND B pinned at walls (B=30). " + NOTE_B),
    ("MC2.05.21.24_LaCl3", S3, "STAGE 3 · LaCl₃", "POOR",
     "Wrong-model bias — model can't track cᶠ. " + NOTE_B),
    ("MC4.07.11.24_SLaCl3", S3, "STAGE 3 · LaCl₃", "POOR",
     "Replicate wrong-model bias. " + NOTE_B),
]
for sheet, col, tag, verdict, note in fit_order:
    builders.append(lambda sh=sheet, c=col, t=tag, v=verdict, n=note: s_fit(sh, c, t, v, n))
# section 05 — the takeaway
builders.append(lambda: s_divider("05 · THE TAKEAWAY", "Combine what we have, then design what we're missing",
                                   "Q1 — a better initial guess from the existing runs   ·   Q2 — the experiments that close the gaps", S1))
builders.append(s_synthesis)   # Q1 — combine
builders.append(s_data)        # Q2 — design / suggested experiments
builders.append(s_questions)
builders.append(s_appendix_B)  # appendix at the very end

TOTAL = len(builders)
for i, b in enumerate(builders, 1):
    s = b()
    # footer on non-title/non-divider light slides handled generically:
    # add footer to every slide except pure title (slide 1) and dividers (navy)
    # detect navy bg
    try:
        is_navy = (s.background.fill.type is not None and
                   s.background.fill.fore_color.rgb == NAVY)
    except Exception:
        is_navy = False
    if i != 1 and not is_navy:
        footer(s, i, TOTAL)

prs.save(str(OUT))
print(f"saved {OUT}  ({TOTAL} slides)")
