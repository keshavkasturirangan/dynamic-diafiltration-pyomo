"""DATA3_contour_diagnostics_2026-06-04.pptx — advisor-discussion deck.

NARRATIVE ARC (three stages of parameter identifiability):

  Stage 1 — WELL-IDENTIFIED         : MC2 NaCl   σ × Lp  (4 of 5 channels converge)
  Stage 2 — INCONSISTENT BUT MIN'D  : MC2 CaCl₂  σ × Lp  (channels split by side)
  Stage 3 — NOT ESTIMABLE            : MC2 LaCl₃  σ × Lp  (no consistent basin)

The deck builds from well-behaved to pathological, showing the progression
of parameter identifiability with cation valence on NF270.

Supporting slides after the main narrative provide the deeper view:
  - cross-sheet argmin summary
  - NaCl in all 3 parameter pairs (the σ × Lp 'well-identified' case has
    subtle wrinkles when you also sweep B × Lp and B × σ)
  - 3D L_p × B × σ slice composites for NaCl + CaCl₂
"""

from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE


REPO = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo")
DST  = REPO / "refactored_codes_v1" / "DATA3_contour_diagnostics_2026-06-04.pptx"
ART  = REPO / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"


INK     = RGBColor(0x1A, 0x1A, 0x1A)
GREY    = RGBColor(0x55, 0x55, 0x55)
LIGHT   = RGBColor(0x8A, 0x8A, 0x8A)
ACCENT  = RGBColor(0xC0, 0x39, 0x2B)
SOFT_BG = RGBColor(0xF6, 0xF8, 0xFB)
SOFT_B  = RGBColor(0xCB, 0xD5, 0xE1)
HEAD_BG = RGBColor(0xEC, 0xEC, 0xEC)
ALT_BG  = RGBColor(0xF7, 0xF7, 0xF7)

# Stage badge colors
STAGE1_FILL = RGBColor(0xE6, 0xF4, 0xEA)   # soft green
STAGE1_TEXT = RGBColor(0x1E, 0x7C, 0x3F)   # forest green
STAGE2_FILL = RGBColor(0xFD, 0xF3, 0xD7)   # soft amber
STAGE2_TEXT = RGBColor(0xA8, 0x70, 0x00)   # amber
STAGE3_FILL = RGBColor(0xFC, 0xE3, 0xE3)   # soft red
STAGE3_TEXT = RGBColor(0xB0, 0x21, 0x21)   # red


# ---------- helpers ----------
def add_text(slide, x, y, w, h, text, *, size=11, bold=False, color=INK,
             align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, italic=False, font="Calibri"):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.margin_left = tf.margin_right = Emu(36000)
    tf.margin_top = tf.margin_bottom = Emu(18000)
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run(); r.text = line
        r.font.size = Pt(size); r.font.bold = bold; r.font.italic = italic
        r.font.name = font; r.font.color.rgb = color
    return tb


def set_cell(cell, text, *, size=10, bold=False, color=INK, fill=None,
             align=PP_ALIGN.LEFT, italic=False):
    cell.margin_left = cell.margin_right = Emu(50000)
    cell.margin_top = cell.margin_bottom = Emu(30000)
    cell.vertical_anchor = MSO_ANCHOR.MIDDLE
    if fill is not None:
        cell.fill.solid(); cell.fill.fore_color.rgb = fill
    tf = cell.text_frame; tf.word_wrap = True; tf.text = ""
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run(); r.text = line
        r.font.size = Pt(size); r.font.bold = bold; r.font.italic = italic
        r.font.name = "Calibri"; r.font.color.rgb = color


def add_title_block(slide, label, title, subtitle=""):
    add_text(slide, 0.45, 0.30, 0.85, 0.55, label,
             size=24, bold=True, color=LIGHT)
    add_text(slide, 1.35, 0.30, 11.5, 0.55, title,
             size=22, bold=True, color=INK)
    if subtitle:
        add_text(slide, 1.35, 0.90, 11.5, 0.40, subtitle,
                 size=13, italic=True, color=GREY)


def add_stage_badge(slide, stage_num, label_text, x=11.5, y=0.30):
    """Add a colored stage badge to top-right of slide."""
    if stage_num == 1:
        fill, text_color = STAGE1_FILL, STAGE1_TEXT
    elif stage_num == 2:
        fill, text_color = STAGE2_FILL, STAGE2_TEXT
    else:
        fill, text_color = STAGE3_FILL, STAGE3_TEXT
    badge = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                   Inches(x), Inches(y), Inches(1.55), Inches(0.50))
    badge.fill.solid(); badge.fill.fore_color.rgb = fill
    badge.line.color.rgb = text_color; badge.line.width = Pt(1.0)
    btf = badge.text_frame
    btf.margin_left = btf.margin_right = Emu(36000)
    btf.margin_top = btf.margin_bottom = Emu(18000)
    btf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = btf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = f"STAGE {stage_num} · {label_text}"
    r.font.size = Pt(10); r.font.bold = True; r.font.color.rgb = text_color


def add_image_centered(slide, png_path, *, y_top=1.50, max_height=5.20, max_width=12.40):
    from PIL import Image
    img = Image.open(png_path)
    iw, ih = img.size
    aspect = iw / ih
    if max_width / aspect <= max_height:
        w = max_width; h = max_width / aspect
    else:
        h = max_height; w = max_height * aspect
    x = (13.333 - w) / 2.0
    slide.shapes.add_picture(str(png_path), Inches(x), Inches(y_top),
                             width=Inches(w), height=Inches(h))
    return (x, y_top + h)


def add_takeaway_box(slide, y, text, *, height=0.95, color=ACCENT):
    box = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                  Inches(0.45), Inches(y),
                                  Inches(12.4), Inches(height))
    box.fill.solid(); box.fill.fore_color.rgb = SOFT_BG
    box.line.color.rgb = color; box.line.width = Pt(1.0)
    tf = box.text_frame
    tf.margin_left = tf.margin_right = Emu(72000)
    tf.margin_top = tf.margin_bottom = Emu(36000)
    tf.word_wrap = True
    tf.vertical_anchor = MSO_ANCHOR.TOP
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.LEFT
    for txt, bold, c in [("Takeaway:  ", True, color),
                          (text,        False, INK)]:
        r = p.add_run(); r.text = txt
        r.font.size = Pt(11); r.font.bold = bold; r.font.color.rgb = c


# ---------- build ----------
prs = Presentation()
prs.slide_width  = Inches(13.333333)
prs.slide_height = Inches(7.5)
blank = prs.slide_layouts[6]


# ======================== Slide 1 — title ============================
s = prs.slides.add_slide(blank)
add_text(s, 0.5, 1.6, 12.3, 1.0,
         "DATA3 — Contour Diagnostics",
         size=44, bold=True, color=INK, align=PP_ALIGN.CENTER)
add_text(s, 0.5, 2.85, 12.3, 0.6,
         "Three stages of parameter identifiability on NF270 single-salt fits",
         size=18, italic=True, color=GREY, align=PP_ALIGN.CENTER)

# Three-stage roadmap teaser
roadmap_y = 4.10
roadmap_h = 1.50
for i, (stage_num, label, sheet, descr) in enumerate([
    (1, "Well-identified", "MC2 NaCl",
     "4 of 5 objective channels\nconverge at the warm-start θ"),
    (2, "Inconsistent",     "MC2 CaCl₂",
     "Channels each find a minimum,\nbut at opposite (σ, Lp) walls"),
    (3, "Not estimable",    "MC2 LaCl₃",
     "Every channel argmins at a\ndifferent Lp — no shared basin"),
]):
    x = 0.5 + i * 4.3
    if stage_num == 1:
        fill, text_c = STAGE1_FILL, STAGE1_TEXT
    elif stage_num == 2:
        fill, text_c = STAGE2_FILL, STAGE2_TEXT
    else:
        fill, text_c = STAGE3_FILL, STAGE3_TEXT
    box = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                              Inches(x), Inches(roadmap_y),
                              Inches(4.10), Inches(roadmap_h))
    box.fill.solid(); box.fill.fore_color.rgb = fill
    box.line.color.rgb = text_c; box.line.width = Pt(1.0)
    add_text(s, x, roadmap_y + 0.10, 4.10, 0.35,
             f"STAGE {stage_num}",
             size=11, bold=True, color=text_c, align=PP_ALIGN.CENTER)
    add_text(s, x, roadmap_y + 0.42, 4.10, 0.40,
             label,
             size=15, bold=True, color=INK, align=PP_ALIGN.CENTER)
    add_text(s, x, roadmap_y + 0.82, 4.10, 0.30,
             sheet,
             size=11, italic=True, color=GREY, align=PP_ALIGN.CENTER)
    add_text(s, x, roadmap_y + 1.12, 4.10, 0.40,
             descr,
             size=10, color=INK, align=PP_ALIGN.CENTER)

add_text(s, 0.5, 6.6, 12.3, 0.4,
         "Discussion · 2026-06-04",
         size=12, italic=True, color=LIGHT, align=PP_ALIGN.CENTER)


# ======================== Slide 2 — context ==========================
s = prs.slides.add_slide(blank)
add_title_block(s, "01", "01 · CONTEXT — parameter bounds",
                "L_p ∈ [0.5, 50] and σ ∈ [0, 1] inherited from DATA1/DATA2 KCl;  per-salt B bounds")

kcl_y = 1.55
kcl_box = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                              Inches(0.45), Inches(kcl_y), Inches(12.4), Inches(0.75))
kcl_box.fill.solid(); kcl_box.fill.fore_color.rgb = SOFT_BG
kcl_box.line.color.rgb = SOFT_B
kcl_tf = kcl_box.text_frame
kcl_tf.margin_left = kcl_tf.margin_right = Emu(72000)
kcl_tf.vertical_anchor = MSO_ANCHOR.MIDDLE
p = kcl_tf.paragraphs[0]; p.alignment = PP_ALIGN.LEFT
for txt, bold in [("KCl baseline (inherited):  ", True),
                  ("L_p ∈ [0.5, 50] L/(m²·hr·bar),    B ∈ [0, 30] µm/s,    σ ∈ [0, 1]",
                   False)]:
    r = p.add_run(); r.text = txt
    r.font.size = Pt(12); r.font.bold = bold; r.font.color.rgb = INK

add_text(s, 0.45, 2.55, 12.4, 0.30,
         "Per-salt B upper bound (physics-informed — halve per +1 valence):",
         size=12, bold=True, color=ACCENT)
tbl = s.shapes.add_table(4, 3,
                         Inches(0.45), Inches(2.90),
                         Inches(12.4), Inches(2.10)).table
for i, w in enumerate([2.5, 2.0, 7.9]):
    tbl.columns[i].width = Inches(w)
tbl.rows[0].height = Inches(0.45)
for i in range(1, 4):
    tbl.rows[i].height = Inches(0.55)
for j, h in enumerate(["Salt", "B bound [µm/s]", "Basis"]):
    set_cell(tbl.cell(0, j), h, size=11.5, bold=True, fill=HEAD_BG,
             align=PP_ALIGN.CENTER)
rows = [
    ("NaCl",  "[0, 30]", "monovalent — same as KCl"),
    ("CaCl₂", "[0, 15]", "divalent counter-ion — halved"),
    ("LaCl₃", "[0, 10]", "trivalent — third"),
]
for i, (a, b, c) in enumerate(rows, start=1):
    fill = ALT_BG if i % 2 == 1 else None
    set_cell(tbl.cell(i, 0), a, size=11.5, bold=True, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(tbl.cell(i, 1), b, size=11.5, bold=True, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(tbl.cell(i, 2), c, size=10.5, fill=fill)

add_text(s, 0.45, 5.30, 12.4, 1.80,
         "How the contour panels below are built:  for every (σ, Lp) grid point we forward-simulate "
         "the full Spiegler-Kedem DAE (B and the state initialization fixed at warm-start values), "
         "then compute the weighted sum-of-squared-residuals per channel.  5 channels per panel:  "
         "mass + permeate concentration + retentate concentration + permeate conductivity + retentate "
         "conductivity.  No fitting — pure forward-simulation diagnostics.",
         size=11, color=INK)


# ====================== STAGE 1 — well-identified =====================
s = prs.slides.add_slide(blank)
add_title_block(s, "02",
                "02 · MC2 NaCl — Well-identified case",
                "B fixed at warm-start (9.73 µm/s) · σ × Lp swept · all 5 channel SSRs evaluated at every cell")
add_stage_badge(s, 1, "Well-identified")
add_image_centered(s, ART / "matlab_ports_5ch_allpairs" / "MC2.05.07.24_NaCl" / "objcontour-x_sigma-y_Lp.png",
                   y_top=1.55, max_height=4.40, max_width=12.6)
add_takeaway_box(s, 6.10,
    "4 of 5 channels converge at the warm-start θ = (σ=1, Lp=9.95):  mass, permeate concentration, "
    "retentate concentration, and retentate conductivity all argmin at the same point.  The five "
    "channels collectively pin down (σ, Lp) — the parameters are well-identified by the data.  "
    "(Permeate conductivity wants Lp ≈ 15.9 — a subtle wrinkle revisited later.)",
    color=STAGE1_TEXT)


# ====================== STAGE 2 — inconsistent ========================
s = prs.slides.add_slide(blank)
add_title_block(s, "03",
                "03 · MC2 CaCl₂ — Inconsistent across channels",
                "B fixed at warm-start (0.51 µm/s) · σ × Lp swept · channels split by side of the membrane")
add_stage_badge(s, 2, "Inconsistent")
add_image_centered(s, ART / "matlab_ports_5ch" / "MC2.05.07.24_CaCl2" / "objcontour-x_sigma-y_Lp.png",
                   y_top=1.55, max_height=4.40, max_width=12.6)
add_takeaway_box(s, 6.10,
    "Each channel finds its own minimum — the objective is minimized, but at inconsistent (σ, Lp).  "
    "Both permeate channels (concentration + conductivity) argmin at σ=0 wall, Lp=8.  "
    "Both retentate channels argmin at σ=1 wall, Lp=4.  Mass disagrees with all — interior σ ≈ 0.53.  "
    "The channels can be fit individually but no shared θ satisfies all of them.",
    color=STAGE2_TEXT)


# ====================== STAGE 3 — not estimable =======================
s = prs.slides.add_slide(blank)
add_title_block(s, "04",
                "04 · MC2 LaCl₃ — Parameters not estimable",
                "B fixed at warm-start (0.31 µm/s) · σ × Lp swept · every channel argmins at a different Lp")
add_stage_badge(s, 3, "Not estimable")
add_image_centered(s, ART / "matlab_ports_5ch" / "MC2.05.21.24_LaCl3" / "objcontour-x_sigma-y_Lp.png",
                   y_top=1.55, max_height=4.40, max_width=12.6)
add_takeaway_box(s, 6.10,
    "Every channel argmins at a different Lp:  mass=4.98, perm conc=0.75, ret conc=2.49, perm cond=1.25, "
    "ret cond=2.99.  The objective is bumpy across the entire (σ, Lp) plane — no consistent basin.  "
    "Single-point ParmEst will pick whichever channel dominates the weighted sum.  "
    "Path forward:  ion-pairing physics (LaCl²⁺/LaCl⁺ equilibria) + σ(cF) sub-model — slide 12 of v10 deck.",
    color=STAGE3_TEXT)


# ====================== Cross-sheet summary ===========================
s = prs.slides.add_slide(blank)
add_title_block(s, "05",
                "05 · CROSS-SHEET SUMMARY — argmin per channel",
                "valence progression — the same σ × Lp sweep applied to three sheets")

add_text(s, 0.45, 1.55, 12.4, 0.30,
         "Argmin L_p [L/(m²·hr·bar)] per channel",
         size=12, bold=True, color=ACCENT)
tbl1 = s.shapes.add_table(4, 7,
                          Inches(0.45), Inches(1.92),
                          Inches(12.4), Inches(2.00)).table
for i, w in enumerate([1.50, 1.40, 1.70, 1.80, 1.80, 1.80, 2.40]):
    tbl1.columns[i].width = Inches(w)
tbl1.rows[0].height = Inches(0.50)
for i in range(1, 4):
    tbl1.rows[i].height = Inches(0.50)
heads = ["Sheet", "Mass", "Perm conc", "Ret conc", "Perm cond", "Ret cond", "Warm-start Lp"]
for j, h in enumerate(heads):
    set_cell(tbl1.cell(0, j), h, size=11, bold=True, fill=HEAD_BG, align=PP_ALIGN.CENTER)
rows1 = [
    ("NaCl ✓",  "9.95", "9.95", "9.95", "15.93 ↑",      "9.95", "9.95"),
    ("CaCl₂ ⚠", "8.00", "8.00", "4.00",  "8.00",        "4.00", "4.00"),
    ("LaCl₃ ✗", "4.98", "0.75", "2.49",  "1.25",        "2.99", "2.49"),
]
for i, row in enumerate(rows1, start=1):
    fill = ALT_BG if i % 2 == 1 else None
    for j, v in enumerate(row):
        bold = j == 0
        color = ACCENT if "↑" in v else INK
        if j == 0:
            if "✓" in v: color = STAGE1_TEXT
            elif "⚠" in v: color = STAGE2_TEXT
            elif "✗" in v: color = STAGE3_TEXT
        set_cell(tbl1.cell(i, j), v, size=11, bold=bold, color=color, fill=fill,
                 align=PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER)

add_text(s, 0.45, 4.10, 12.4, 0.30,
         "Argmin σ per channel",
         size=12, bold=True, color=ACCENT)
tbl2 = s.shapes.add_table(4, 6,
                          Inches(0.45), Inches(4.47),
                          Inches(12.4), Inches(2.00)).table
for i, w in enumerate([1.50, 1.55, 2.10, 2.10, 2.50, 2.65]):
    tbl2.columns[i].width = Inches(w)
tbl2.rows[0].height = Inches(0.50)
for i in range(1, 4):
    tbl2.rows[i].height = Inches(0.50)
heads2 = ["Sheet", "Mass", "Perm conc", "Ret conc", "Perm cond", "Ret cond"]
for j, h in enumerate(heads2):
    set_cell(tbl2.cell(0, j), h, size=11, bold=True, fill=HEAD_BG, align=PP_ALIGN.CENTER)
rows2 = [
    ("NaCl ✓",  "1.00",         "1.00",         "1.00",         "1.00",         "1.00"),
    ("CaCl₂ ⚠", "0.53 interior","0.00 wall",    "1.00 wall",    "0.00 wall",    "1.00 wall"),
    ("LaCl₃ ✗", "0.16 interior","1.00 wall",    "1.00 wall",    "0.00 wall",    "1.00 wall"),
]
for i, row in enumerate(rows2, start=1):
    fill = ALT_BG if i % 2 == 1 else None
    for j, v in enumerate(row):
        bold = j == 0
        color = ACCENT if "wall" in v else INK
        if j == 0:
            if "✓" in v: color = STAGE1_TEXT
            elif "⚠" in v: color = STAGE2_TEXT
            elif "✗" in v: color = STAGE3_TEXT
        set_cell(tbl2.cell(i, j), v, size=10.5, bold=bold, color=color, fill=fill,
                 align=PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER)

add_text(s, 0.45, 6.70, 12.4, 0.65,
         "Identifiability degrades monotonically with cation valence:  NaCl (well-identified) → "
         "CaCl₂ (inconsistent, σ pinned at opposite walls per channel) → LaCl₃ (not estimable).",
         size=11, italic=True, color=GREY)


# ====================== Even the well-identified case has limits ======
s = prs.slides.add_slide(blank)
add_title_block(s, "06",
                "06 · MC2 NaCl REVISITED — limits of 'well-identified'",
                "what if we also sweep B (in addition to σ × Lp at fixed B)?")
add_stage_badge(s, 1, "Well-identified")

# Subtitle / framing
add_text(s, 0.45, 1.55, 12.4, 0.40,
         "Argmin (x, y, log₁₀ SSR) per channel × parameter-pair sweep",
         size=12, bold=True, color=ACCENT)
nacl_tbl = s.shapes.add_table(4, 6,
                              Inches(0.45), Inches(2.00),
                              Inches(12.4), Inches(3.55)).table
for i, w in enumerate([2.20, 2.04, 2.04, 2.04, 2.04, 2.04]):
    nacl_tbl.columns[i].width = Inches(w)
nacl_tbl.rows[0].height = Inches(0.55)
for i in range(1, 4):
    nacl_tbl.rows[i].height = Inches(1.00)

for j, h in enumerate(["Pair (x × y)", "Mass", "Perm conc", "Ret conc",
                       "Perm cond", "Ret cond"]):
    set_cell(nacl_tbl.cell(0, j), h, size=11, bold=True, fill=HEAD_BG, align=PP_ALIGN.CENTER)

sigLp_row = [
    "σ × Lp\n(B fixed = 9.73)",
    "σ=1.00\nLp=9.95\nlog₁₀=0.96",
    "σ=1.00\nLp=9.95\nlog₁₀=2.01",
    "σ=1.00\nLp=9.95\nlog₁₀=3.63",
    "σ=1.00\nLp=15.93 ↑\nlog₁₀=3.33",
    "σ=1.00\nLp=9.95\nlog₁₀=1.60",
]
BLp_row = [
    "B × Lp\n(σ fixed = 1.00)",
    "B=4.74\nLp=10.95\nlog₁₀=0.76",
    "B=26.84\nLp=6.97\nlog₁₀=1.95",
    "B=30.00 ↑\nLp=10.95\nlog₁₀=3.62",
    "B=1.58\nLp=19.91 ↑\nlog₁₀=2.40",
    "B=30.00 ↑\nLp=10.95\nlog₁₀=1.59",
]
Bsig_row = [
    "B × σ\n(Lp fixed = 9.95)",
    "B=0.00\nσ=0.32 (int)\nlog₁₀=0.90",
    "B=9.47\nσ=1.00\nlog₁₀=2.01",
    "B=17.37\nσ=1.00\nlog₁₀=3.63",
    "B=1.58\nσ=0.00\nlog₁₀=2.70",
    "B=15.79\nσ=1.00\nlog₁₀=1.60",
]
for i, row in enumerate([sigLp_row, BLp_row, Bsig_row], start=1):
    fill = ALT_BG if i % 2 == 1 else None
    for j, v in enumerate(row):
        bold = j == 0
        color = ACCENT if "↑" in v else INK
        set_cell(nacl_tbl.cell(i, j), v, size=9, bold=bold, color=color, fill=fill,
                 align=PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER)

add_takeaway_box(s, 5.85,
    "When σ is pinned at the σ=1 wall (B × Lp row), both retentate channels want B = 30 — the upper bound.  "
    "When Lp is pinned at the warm-start (B × σ row), mass argmins at interior σ ≈ 0.32 with B near 0.  "
    "Even the 'well-identified' σ × Lp picture has subtle limits — the σ=1 wall on Stage 1 was hiding "
    "a B-vs-Lp trade-off.  Worth flagging:  the NaCl identifiability is good but not perfect.",
    height=1.30, color=STAGE1_TEXT)


# ====================== Supporting: NaCl B × Lp + B × σ contours ======
for label, title, sub, png in [
    ("07A", "07A · MC2 NaCl B × Lp contour",
     "σ fixed at warm-start (σ=1) · B × Lp swept · supports slide 06",
     ART / "matlab_ports_5ch_allpairs" / "MC2.05.07.24_NaCl" / "objcontour-x_B-y_Lp.png"),
    ("07B", "07B · MC2 NaCl B × σ contour",
     "Lp fixed at warm-start (Lp=9.95) · B × σ swept · supports slide 06",
     ART / "matlab_ports_5ch_allpairs" / "MC2.05.07.24_NaCl" / "objcontour-x_B-y_sigma.png"),
]:
    if not png.exists():
        continue
    s = prs.slides.add_slide(blank)
    add_title_block(s, label, title, sub)
    add_stage_badge(s, 1, "Well-identified")
    add_image_centered(s, png, y_top=1.55, max_height=5.40, max_width=12.6)


# ====================== 3D slice composites ===========================
for label, title, sub, png in [
    ("08A", "08A · 3D L_p × B × σ — NaCl",
     "10 σ-slices × 3 channels.  Confirms the σ × Lp basin structure of slide 02.",
     ART / "contour3d" / "MC2.05.07.24_NaCl" / "MC2.05.07.24_NaCl_grid3d_slices.png"),
    ("08B", "08B · 3D L_p × B × σ — CaCl₂",
     "10 σ-slices × 3 channels.  Note how basin shape morphs between σ=0 and σ=1 — the side-split of slide 03.",
     ART / "contour3d" / "MC2.05.07.24_CaCl2" / "MC2.05.07.24_CaCl2_grid3d_slices.png"),
]:
    if not png.exists():
        continue
    s = prs.slides.add_slide(blank)
    add_title_block(s, label, title, sub)
    add_image_centered(s, png, y_top=1.55, max_height=5.85, max_width=12.6)


# ====================== Discussion ====================================
s = prs.slides.add_slide(blank)
add_title_block(s, "09",
                "09 · QUESTIONS FOR DISCUSSION",
                "what I'd like your sanity check on before publishing")

questions = [
    ("1. ", "Three-stage diagnostic story",
     "Does the well-identified → inconsistent → not-estimable framing match what you'd intuit from "
     "the underlying physics (anionic NF270 + chloride co-ion + increasing cation valence)?"),
    ("2. ", "Permeate conductivity divergence on NaCl",
     "Even the well-identified Stage 1 has permeate conductivity wanting Lp ≈ 15.9 vs all others at 9.95.  "
     "Is this (a) m.cH-vs-cV physics that the warm-start under-fits, or (b) Shedlovsky-forward calibration?"),
    ("3. ", "Per-salt B bounds in the optimizer",
     "We've landed NaCl [0,30], CaCl₂ [0,15], LaCl₃ [0,10] in the library.  "
     "Anything you'd revise before re-running the full fit campaign?"),
    ("4. ", "LaCl₃ ion-pairing fix",
     "Stage 3 confirms slide 12's 'physics gap' diagnosis.  Should we (a) parameterize σ(cF) explicitly, "
     "(b) add an ion-pairing equilibrium block, or (c) fit only the diluting-regime portion?"),
    ("5. ", "Which subset for the next campaign meeting?",
     "Lead with the three-stage roadmap (slide 1) + the three Stage slides?  Or save Stage 1 wrinkles "
     "(slide 06) and the 3D appendix for follow-up?"),
]
y = 1.55
for n, title, body in questions:
    add_text(s, 0.45, y, 0.5, 0.40, n,
             size=14, bold=True, color=ACCENT)
    add_text(s, 0.85, y, 12.0, 0.40, title,
             size=13, bold=True, color=INK)
    add_text(s, 0.85, y + 0.35, 12.0, 0.80, body,
             size=10.5, color=GREY)
    y += 1.05


# Save
prs.save(DST)
print(f"WROTE  {DST}")
print(f"size:   {DST.stat().st_size:>11,} bytes")
print(f"slides: {len(prs.slides)}")
