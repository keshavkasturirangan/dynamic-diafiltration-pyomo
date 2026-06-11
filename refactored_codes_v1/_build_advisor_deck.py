"""DATA3_contour_diagnostics_2026-06-04_v2.pptx — advisor deck (Path A).

Built immediately with what we have on disk + the workflow's plain-language
per-salt analysis (which is grounded in warm-start θ patterns across all 11
sheets, so it covers the campaign as a whole without needing every contour PNG).

Slide order follows the 3-stage narrative:

   STAGE 1 — well-identified         : MC3.07.22.24_SNaCl (interior σ ≈ 0.45)
   STAGE 2 — inconsistent but min'd  : NaCl (Lp tight, σ bimodal), CaCl2 (channels split by side)
   STAGE 3 — not estimable            : LaCl3 (no consistent basin + n=2 problem)
   + cross-sheet warm-start summary, "what data should we collect", questions

PNGs available on disk:
   matlab_ports_5ch_full/MC2.05.07.24_NaCl/objcontour-x_sigma-y_Lp.png   (15x15 from overnight run)
   matlab_ports_5ch/MC2.05.07.24_NaCl/objcontour-x_sigma-y_Lp.png         (20x20 earlier)
   matlab_ports_5ch_allpairs/MC2.05.07.24_NaCl/objcontour-x_{B-y_Lp,B-y_sigma}.png   (B sweeps for NaCl)
   matlab_ports_5ch/MC2.05.07.24_CaCl2/objcontour-x_sigma-y_Lp.png
   matlab_ports_5ch/MC2.05.21.24_LaCl3/objcontour-x_sigma-y_Lp.png
   contour3d/MC2.05.07.24_{NaCl,CaCl2}/*_grid3d_slices.png
"""

from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE


REPO = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo")
DST  = REPO / "refactored_codes_v1" / "DATA3_contour_diagnostics_2026-06-04_v5.pptx"
ART  = REPO / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
FULL = ART / "matlab_ports_5ch_full"
ANIM = ART / "animations3d"          # 3D scrubbing animations (GIFs)
LOGO = REPO / "refactored_codes_v1" / "notre_dame_logo.png"


INK     = RGBColor(0x1A, 0x1A, 0x1A)
GREY    = RGBColor(0x55, 0x55, 0x55)
LIGHT   = RGBColor(0x8A, 0x8A, 0x8A)
ACCENT  = RGBColor(0xC0, 0x39, 0x2B)
SOFT_BG = RGBColor(0xF6, 0xF8, 0xFB)
SOFT_B  = RGBColor(0xCB, 0xD5, 0xE1)
HEAD_BG = RGBColor(0xEC, 0xEC, 0xEC)
ALT_BG  = RGBColor(0xF7, 0xF7, 0xF7)

STAGE1_FILL = RGBColor(0xE6, 0xF4, 0xEA)
STAGE1_TEXT = RGBColor(0x1E, 0x7C, 0x3F)
STAGE2_FILL = RGBColor(0xFD, 0xF3, 0xD7)
STAGE2_TEXT = RGBColor(0xA8, 0x70, 0x00)
STAGE3_FILL = RGBColor(0xFC, 0xE3, 0xE3)
STAGE3_TEXT = RGBColor(0xB0, 0x21, 0x21)


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
    # Strip any leading "<label> · " from the title to avoid double-numbering
    # (the section label is already shown as the big badge on the left).
    # Handles "01", "10A", "A1" — anything alphanumeric followed by a bullet.
    import re as _re
    title_clean = _re.sub(r"^[A-Za-z0-9]+\s*[·•]\s*", "", title).strip()
    add_text(slide, 0.45, 0.30, 0.85, 0.55, label,
             size=24, bold=True, color=LIGHT)
    add_text(slide, 1.35, 0.30, 11.5, 0.55, title_clean,
             size=22, bold=True, color=INK)
    if subtitle:
        add_text(slide, 1.35, 0.90, 11.5, 0.40, subtitle,
                 size=13, italic=True, color=GREY)


def add_nd_logo(slide, *, width_in=1.40):
    """Add the Notre Dame logo as a footer on the bottom-left of the slide."""
    from PIL import Image as _Image
    img = _Image.open(LOGO)
    aspect = img.size[0] / img.size[1]   # 750/175 ≈ 4.286
    h = width_in / aspect
    # Bottom-left corner, with small margin
    x = 0.30
    y = 7.50 - h - 0.15
    slide.shapes.add_picture(str(LOGO), Inches(x), Inches(y),
                              width=Inches(width_in), height=Inches(h))


def add_stage_badge(slide, stage_num, label_text, x=11.5, y=0.30):
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


def add_takeaway_box(slide, y, text, *, height=0.95, color=ACCENT, label="Takeaway"):
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
    for txt, bold, c in [(f"{label}:  ", True, color),
                          (text,        False, INK)]:
        r = p.add_run(); r.text = txt
        r.font.size = Pt(11); r.font.bold = bold; r.font.color.rgb = c


# ---------- build ----------
prs = Presentation()
prs.slide_width  = Inches(13.333333)
prs.slide_height = Inches(7.5)
blank = prs.slide_layouts[6]


# Monkey-patch slides.add_slide so every new slide automatically gets the
# Notre Dame logo footer (matching the v9 deck convention).
_orig_add_slide = prs.slides.add_slide
def _add_slide_with_logo(layout):
    s = _orig_add_slide(layout)
    try:
        add_nd_logo(s)
    except Exception as e:
        print(f"  [logo] skip: {e}")
    return s
prs.slides.add_slide = _add_slide_with_logo


# ======================== Slide 1 — title ============================
s = prs.slides.add_slide(blank)
add_text(s, 0.5, 1.4, 12.3, 1.0,
         "DATA3 — Contour Diagnostics",
         size=44, bold=True, color=INK, align=PP_ALIGN.CENTER)
add_text(s, 0.5, 2.65, 12.3, 0.5,
         "Three stages of parameter identifiability across all 11 NF270 single-salt sheets",
         size=18, italic=True, color=GREY, align=PP_ALIGN.CENTER)

roadmap_y = 3.65
roadmap_h = 2.00
for i, (stage_num, label, who, descr) in enumerate([
    (1, "Well-identified",  "2 sheets · NaCl:  MC3 SNaCl + MC5 S2NaCl",
     "Lp identified two ways:\n(a) interior σ via concentrating run,\n(b) σ-harmless via diluting run"),
    (2, "Inconsistent",      "7 sheets:  4 NaCl + 3 CaCl₂",
     "Channels each find a minimum,\nbut at different (σ, Lp) walls.\nLp tight, σ rail-pinned."),
    (3, "Not estimable",     "2 sheets:  MC2 + MC4 LaCl₃",
     "Same wrong-model bias on both:\nthe constant-(σ,B) Spiegler-Kedem\ncannot represent La³⁺ ion-pairing."),
]):
    x = 0.4 + i * 4.32
    if stage_num == 1:
        fill, text_c = STAGE1_FILL, STAGE1_TEXT
    elif stage_num == 2:
        fill, text_c = STAGE2_FILL, STAGE2_TEXT
    else:
        fill, text_c = STAGE3_FILL, STAGE3_TEXT
    box = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                              Inches(x), Inches(roadmap_y),
                              Inches(4.20), Inches(roadmap_h))
    box.fill.solid(); box.fill.fore_color.rgb = fill
    box.line.color.rgb = text_c; box.line.width = Pt(1.0)
    add_text(s, x, roadmap_y + 0.13, 4.20, 0.35,
             f"STAGE {stage_num}",
             size=11, bold=True, color=text_c, align=PP_ALIGN.CENTER)
    add_text(s, x, roadmap_y + 0.50, 4.20, 0.40,
             label,
             size=15, bold=True, color=INK, align=PP_ALIGN.CENTER)
    add_text(s, x, roadmap_y + 0.90, 4.20, 0.30,
             who,
             size=10, italic=True, color=GREY, align=PP_ALIGN.CENTER)
    add_text(s, x, roadmap_y + 1.22, 4.20, 0.75,
             descr,
             size=10, color=INK, align=PP_ALIGN.CENTER)

add_text(s, 0.5, 6.20, 12.3, 0.45,
         "What we trust most → what we trust least.   Agreement across the 5 measurement channels = trust.",
         size=11, italic=True, color=ACCENT, align=PP_ALIGN.CENTER)
add_text(s, 0.5, 6.85, 12.3, 0.35,
         "Discussion · 2026-06-04",
         size=12, italic=True, color=LIGHT, align=PP_ALIGN.CENTER)


# ======================== Slide 2 — context + bounds =================
s = prs.slides.add_slide(blank)
add_title_block(s, "01", "01 · CONTEXT — parameter bounds in play",
                "L_p ∈ [0.5, 50] and σ ∈ [0, 1] inherited from DATA1/DATA2 KCl; per-salt B bound is the only new knob")

kcl_y = 1.55
kcl_box = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                              Inches(0.45), Inches(kcl_y), Inches(12.4), Inches(0.75))
kcl_box.fill.solid(); kcl_box.fill.fore_color.rgb = SOFT_BG
kcl_box.line.color.rgb = SOFT_B
kcl_tf = kcl_box.text_frame
kcl_tf.margin_left = kcl_tf.margin_right = Emu(72000)
kcl_tf.vertical_anchor = MSO_ANCHOR.MIDDLE
p = kcl_tf.paragraphs[0]; p.alignment = PP_ALIGN.LEFT
for txt, bold in [("KCl baseline (inherited from DATA1/DATA2):  ", True),
                  ("L_p ∈ [0.5, 50] L/(m²·hr·bar),    B ∈ [0, 30] µm/s,    σ ∈ [0, 1]",
                   False)]:
    r = p.add_run(); r.text = txt
    r.font.size = Pt(12); r.font.bold = bold; r.font.color.rgb = INK

add_text(s, 0.45, 2.55, 12.4, 0.30,
         "Per-salt B upper bound — halve per +1 of cation charge:",
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
    ("CaCl₂", "[0, 15]", "divalent — halved"),
    ("LaCl₃", "[0, 10]", "trivalent — third"),
]
for i, (a, b, c) in enumerate(rows, start=1):
    fill = ALT_BG if i % 2 == 1 else None
    set_cell(tbl.cell(i, 0), a, size=11.5, bold=True, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(tbl.cell(i, 1), b, size=11.5, bold=True, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(tbl.cell(i, 2), c, size=10.5, fill=fill)

add_text(s, 0.45, 5.30, 12.4, 1.80,
         "How the contour panels are built:  at every (σ, Lp) grid point we forward-simulate the full "
         "Spiegler-Kedem DAE (B fixed at warm-start, S0/S also at warm-start), then evaluate the weighted "
         "sum-of-squared-residuals per channel.  5 channels per panel:  mass, permeate concentration, "
         "retentate concentration, permeate conductivity, retentate conductivity.  No fitting — pure "
         "forward-simulation diagnostics.",
         size=11, color=INK)


# ======================== Slide 3 — cross-sheet θ summary ============
s = prs.slides.add_slide(blank)
add_title_block(s, "02", "02 · ALL 11 SHEETS — warm-start θ summary",
                "where the optimizer landed on every sheet — color-coded by stage")

add_text(s, 0.45, 1.55, 12.4, 0.30,
         "Warm-start fits across the campaign (sorted by salt, then by sheet ID):",
         size=12, bold=True, color=ACCENT)

# Table: 12 rows (header + 11 sheets) × 7 cols
tbl = s.shapes.add_table(12, 7,
                         Inches(0.45), Inches(1.92),
                         Inches(12.4), Inches(5.10)).table
for i, w in enumerate([2.30, 1.10, 1.10, 1.20, 1.10, 1.40, 4.20]):
    tbl.columns[i].width = Inches(w)
tbl.rows[0].height = Inches(0.40)
for i in range(1, 12):
    tbl.rows[i].height = Inches(0.42)

heads = ["Sheet", "Salt", "L_p", "B [µm/s]", "σ", "WSSE", "Stage assignment"]
for j, h in enumerate(heads):
    set_cell(tbl.cell(0, j), h, size=10.5, bold=True, fill=HEAD_BG, align=PP_ALIGN.CENTER)

# stage colors per row (1=green, 2=amber, 3=red)
sheets_data = [
    ("MC2.05.07.24_NaCl",   "NaCl",  9.95, 9.73,  "1.00 ↑", 4357.3, 2, "σ wall, concentrating regime"),
    ("MC3.07.22.24_SNaCl",  "NaCl",  8.22, 14.68, "0.45",   555.7,  1, "interior σ — concentrating run probed σ-sensitive regime (gold)"),
    ("MC4.07.11.24_SNaCl",  "NaCl",  7.79, 2.02,  "0.00 ↓", 9536.2, 2, "σ=0, diluting regime BUT poor WSSE — data-quality issue"),
    ("MC5.07.23.24_NaCl",   "NaCl",  10.02, 3.70, "1.00 ↑", 5793.6, 2, "σ wall, concentrating regime"),
    ("MC5.07.23.24_SNaCl",  "NaCl",  9.45, 3.07,  "1.00 ↑", 4763.6, 2, "σ wall, concentrating regime"),
    ("MC5.07.23.24_S2NaCl", "NaCl",  7.38, 10.14, "0.00 ↓", 829.6,  1, "σ=0 harmless — diluting run, Lp identified by mass (2nd best WSSE)"),
    ("MC2.05.07.24_CaCl2",  "CaCl₂", 4.00, 0.51,  "1.00 ↑", 6717.1, 2, "channels split by side"),
    ("MC3.07.11.24_SCaCl2", "CaCl₂", 5.53, 9.11,  "1.00 ↑", 4079.7, 2, "σ at wall"),
    ("MC3.07.12.24_S2CaCl2","CaCl₂", 4.63, 30.00, "1.00 ↑", 9251.7, 2, "σ AND B at walls"),
    ("MC2.05.21.24_LaCl3",  "LaCl₃", 2.49, 0.31,  "1.00 ↑", 17711.2, 3, "wrong model — every channel disagrees"),
    ("MC4.07.11.24_SLaCl3", "LaCl₃", 2.78, 0.17,  "1.00 ↑", 21345.1, 3, "wrong model — replicate of MC2"),
]
for i, (run_id, salt, lp, b, sig, wsse, stage, note) in enumerate(sheets_data, start=1):
    if stage == 1:
        fill = STAGE1_FILL; tc = STAGE1_TEXT
    elif stage == 2:
        fill = STAGE2_FILL if i % 2 == 1 else None; tc = STAGE2_TEXT
    else:
        fill = STAGE3_FILL; tc = STAGE3_TEXT
    set_cell(tbl.cell(i, 0), run_id, size=9.5, bold=True, fill=fill, color=tc, align=PP_ALIGN.LEFT)
    set_cell(tbl.cell(i, 1), salt,   size=9.5, fill=fill, color=GREY, align=PP_ALIGN.CENTER)
    set_cell(tbl.cell(i, 2), f"{lp:.2f}", size=10, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(tbl.cell(i, 3), f"{b:.2f}",  size=10, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(tbl.cell(i, 4), sig,         size=10, fill=fill,
             color=ACCENT if "↑" in sig or "↓" in sig else INK, align=PP_ALIGN.CENTER)
    set_cell(tbl.cell(i, 5), f"{wsse:.0f}", size=10, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(tbl.cell(i, 6), f"Stage {stage} · {note}", size=9, fill=fill, color=tc)

add_text(s, 0.45, 7.10, 12.4, 0.35,
         "↑ / ↓ marks parameters pinned at a bound.  Only MC3.07.22.24_SNaCl landed at an interior σ — "
         "that's our reference case for 'well-identified'.",
         size=10, italic=True, color=GREY)


# ======================== Slide 4 — Two paths to well-identified ====
s = prs.slides.add_slide(blank)
add_title_block(s, "03",
                "03 · TWO PATHS TO WELL-IDENTIFIED",
                "the same Lp can be pinned by the data two completely different ways")
add_stage_badge(s, 1, "Well-identified")

add_text(s, 0.45, 1.55, 12.4, 1.20,
         "There are two physically distinct ways the data can pin (Lp, σ, B):  "
         "either we run conditions that PROBE σ directly so the contour shows an interior σ minimum, "
         "OR we run conditions where σ doesn't matter so we don't need to probe it.  Both give us a usable Lp.  "
         "The campaign has one sheet in each category.",
         size=12, color=INK)

# Two large side-by-side cards
panel_y = 2.95
panel_h = 3.30

# Card 1 — Probe σ directly (MC3 SNaCl)
left_box = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                               Inches(0.45), Inches(panel_y),
                               Inches(6.10), Inches(panel_h))
left_box.fill.solid(); left_box.fill.fore_color.rgb = STAGE1_FILL
left_box.line.color.rgb = STAGE1_TEXT; left_box.line.width = Pt(1.0)

add_text(s, 0.65, panel_y + 0.15, 5.70, 0.35,
         "PATH A · probe σ directly",
         size=12, bold=True, color=STAGE1_TEXT)
add_text(s, 0.65, panel_y + 0.50, 5.70, 0.40,
         "MC3.07.22.24_SNaCl",
         size=15, bold=True, color=INK)
add_text(s, 0.65, panel_y + 0.95, 5.70, 0.30,
         "Lp = 8.22    B = 14.68    σ = 0.45    WSSE = 556",
         size=10.5, italic=True, color=INK)
add_text(s, 0.65, panel_y + 1.35, 5.70, 1.85,
         "Concentrating regime, but the run conditions excited σ — probably a wider pressure-driving-force sweep, lower feed, "
         "or a pressure step.  σ leaves a fingerprint in the data, so the optimizer finds an interior minimum at σ = 0.45.  "
         "All 5 channels converge.  This is the textbook well-identified case.",
         size=10.5, color=INK)

# Card 2 — Bypass σ (MC5 S2NaCl)
right_box = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                Inches(6.75), Inches(panel_y),
                                Inches(6.10), Inches(panel_h))
right_box.fill.solid(); right_box.fill.fore_color.rgb = STAGE1_FILL
right_box.line.color.rgb = STAGE1_TEXT; right_box.line.width = Pt(1.0)

add_text(s, 6.95, panel_y + 0.15, 5.70, 0.35,
         "PATH B · bypass σ (diluting regime)",
         size=12, bold=True, color=STAGE1_TEXT)
add_text(s, 6.95, panel_y + 0.50, 5.70, 0.40,
         "MC5.07.23.24_S2NaCl",
         size=15, bold=True, color=INK)
add_text(s, 6.95, panel_y + 0.95, 5.70, 0.30,
         "Lp = 7.38    B = 10.14    σ = 0.00 ↓    WSSE = 830",
         size=10.5, italic=True, color=INK)
add_text(s, 6.95, panel_y + 1.35, 5.70, 1.85,
         "Diluting regime — cF shrinks over the run, so σ · Δπ → 0 by the end no matter what σ is.  "
         "σ literally cannot leave a fingerprint here.  Optimizer drifts to σ = 0 because the gradient is flat, "
         "not because data says σ = 0.  Lp is still pinned tightly by the mass channel.  WSSE = 830 — 2nd best in the campaign.",
         size=10.5, color=INK)

add_takeaway_box(s, 6.45,
    "Same destination, two different roads.  Path A needs careful experimental design; Path B comes free if you happen "
    "to run a diluting experiment.  Either way Lp is trustworthy.  Knowing which path each sheet is on tells us "
    "where the σ value can be trusted (Path A) and where σ is unconstrained but irrelevant (Path B).",
    color=STAGE1_TEXT, height=0.75)


# ======================== Slide 5 — STAGE 1 path A contour ===========
s = prs.slides.add_slide(blank)
add_title_block(s, "04",
                "04 · STAGE 1 (path A) — MC3.07.22.24_SNaCl",
                "interior σ via concentrating run — the gold standard fit (WSSE = 556)")
add_stage_badge(s, 1, "Well-identified")

png = FULL / "MC3.07.22.24_SNaCl" / "objcontour-x_sigma-y_Lp.png"
if png.exists():
    add_image_centered(s, png, y_top=1.55, max_height=4.40, max_width=12.6)

add_takeaway_box(s, 6.10,
    "Look at the mass and permeate-concentration panels — both argmin at interior σ (0.86 and 0.64) near Lp ≈ 8.6.  "
    "Retentate channels still pin at σ=0 (the 0.3 % cF over-weighting effect), but the joint WSSE is dominated by mass + "
    "perm-conc, so the warm-start θ landed interior at σ = 0.45.  This is the only NaCl sheet where σ leaves a real fingerprint.",
    color=STAGE1_TEXT, height=1.10)


# ======================== Slide 6 — STAGE 1 path B contour ===========
s = prs.slides.add_slide(blank)
add_title_block(s, "05",
                "05 · STAGE 1 (path B) — MC5.07.23.24_S2NaCl",
                "diluting run — σ harmless, Lp pinned anyway (WSSE = 830)")
add_stage_badge(s, 1, "Well-identified")

png = FULL / "MC5.07.23.24_S2NaCl" / "objcontour-x_sigma-y_Lp.png"
if png.exists():
    add_image_centered(s, png, y_top=1.55, max_height=4.40, max_width=12.6)

add_takeaway_box(s, 6.10,
    "Four channels argmin at Lp ≈ 7.7 — tight cluster.  σ slides to ~0 because the gradient is flat in the diluting regime "
    "(σ·Δπ → 0 at end of run), but Lp is well-identified by the mass channel.  WSSE = 830 — 2nd best in the entire campaign.  "
    "Diluting experiments are 'free' identifiability for Lp.",
    color=STAGE1_TEXT, height=1.10)


# ======================== Slide 7 — Stage 2 NaCl narrative ===========
s = prs.slides.add_slide(blank)
add_title_block(s, "06",
                "06 · STAGE 2 — NaCl (the 4 concentrating-regime sheets at the σ=1 wall)",
                "Lp tight (8.6 ± 1.1) but σ pinned — these sheets ran concentrating without probing σ-sensitive conditions")
add_stage_badge(s, 2, "Inconsistent")

nacl_png = ART / "matlab_ports_5ch" / "MC2.05.07.24_NaCl" / "objcontour-x_sigma-y_Lp.png"
if nacl_png.exists():
    add_image_centered(s, nacl_png, y_top=1.55, max_height=4.40, max_width=12.6)

add_takeaway_box(s, 6.10,
    "MC2 NaCl shown here — same story on MC5 NaCl, MC5 SNaCl, and the concentrating part of MC4 SNaCl.  "
    "4 of 5 channels converge on Lp ≈ 9.95 (Lp well-identified), but σ pins at the σ=1 wall and the wall is hiding a "
    "B-vs-Lp trade-off.  Concentrating regime + no σ-probing conditions = inconsistent across channels.",
    color=STAGE2_TEXT, height=1.10)


# ======================== Slide 8 — Stage 2 CaCl₂ narrative ==========
s = prs.slides.add_slide(blank)
add_title_block(s, "07",
                "07 · STAGE 2 — CaCl₂ (all 3 sheets)",
                "every sheet pins σ at σ=1 wall in the warm-start fit;  contour shows channels split by membrane side")
add_stage_badge(s, 2, "Inconsistent")

cacl_png = ART / "matlab_ports_5ch" / "MC2.05.07.24_CaCl2" / "objcontour-x_sigma-y_Lp.png"
if cacl_png.exists():
    add_image_centered(s, cacl_png, y_top=1.55, max_height=4.40, max_width=12.6)

add_takeaway_box(s, 6.10,
    "All 3 CaCl₂ sheets pin σ at σ=1 (and one pins B at the upper bound too).  Physics expects interior σ "
    "(divalent counter-ion → weaker Donnan exclusion → lower σ on NF270) — pinning at σ=1 is a weighting problem, "
    "not a model-form problem.  Both permeate channels want σ=0; both retentate channels want σ=1.  "
    "Split by membrane side, not by physics.",
    color=STAGE2_TEXT, height=1.10)


# ======================== Slide 9 — Stage 3 LaCl₃ narrative ==========
s = prs.slides.add_slide(blank)
add_title_block(s, "08",
                "08 · STAGE 3 — LaCl₃ (both sheets)",
                "warm-start θ nearly identical between sheets, but WSSE ~20× worse than NaCl gold standard → shared wrong-model bias")
add_stage_badge(s, 3, "Not estimable")

lacl_png = ART / "matlab_ports_5ch" / "MC2.05.21.24_LaCl3" / "objcontour-x_sigma-y_Lp.png"
if lacl_png.exists():
    add_image_centered(s, lacl_png, y_top=1.55, max_height=4.40, max_width=12.6)

add_takeaway_box(s, 6.10,
    "Both LaCl₃ sheets warm-start to nearly identical θ — but fit the data ~20× worse than NaCl MC3 SNaCl.  "
    "Same wrong-model bias both times, not measurement noise.  The constant-(σ, B) Spiegler-Kedem model "
    "is structurally wrong for trivalent ions — probably because of LaCl²⁺ / LaCl⁺ ion-pairing the model can't see.  "
    "AND we only have 2 sheets:  that's an anecdote, not a sample size.",
    color=STAGE3_TEXT, height=1.20)


# ======================== Slide 8 — what data should we collect ======
s = prs.slides.add_slide(blank)
add_title_block(s, "09", "09 · WHAT DATA SHOULD WE COLLECT",
                "what would unblock the campaign — prioritized by leverage")

# 3 priority panels
priorities = [
    ("PRIORITY 1", STAGE1_FILL, STAGE1_TEXT,
     "Pure-water L_p before and after every run",
     "Single highest-leverage experiment.  Pins L_p independently of the salt-transport data so σ and B don't have "
     "to compensate for it.  Also doubles as a fouling check — if pre-run and post-run L_p drift, we know the "
     "membrane changed during the experiment and the constant-parameter fit was never going to work."),
    ("PRIORITY 2", STAGE2_FILL, STAGE2_TEXT,
     "Multi-pressure sweeps within a single run  +  dilute concentration runs",
     "Pinning σ for CaCl₂ and LaCl₃ requires breaking the σ-B correlation.  σ couples to Δπ, B doesn't — so stepping "
     "ΔP across a wide range during one run changes the relative leverage of σ vs B.  Dilute runs (1-5 mM NaCl, "
     "sub-mM LaCl₃) push into the diffusion-limited regime where σ identifies."),
    ("PRIORITY 3", STAGE3_FILL, STAGE3_TEXT,
     "3-5 more LaCl₃ sheets + retentate-side conductivity probe + UV-Vis / ICP speciation",
     "n = 2 LaCl₃ sheets is not statistics — we need at least 5 before any 'not estimable' claim carries weight.  "
     "Adding high-frequency conductivity on the retentate (we only have it on permeate) eliminates the CaCl₂ side-split "
     "weighting problem.  Speciation directly tests the LaCl²⁺ / LaCl⁺ ion-pairing hypothesis rather than inferring "
     "it from model failure."),
]
y = 1.55
for label, fill, tc, headline, body in priorities:
    box = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                              Inches(0.45), Inches(y),
                              Inches(12.4), Inches(1.65))
    box.fill.solid(); box.fill.fore_color.rgb = fill
    box.line.color.rgb = tc; box.line.width = Pt(1.0)
    add_text(s, 0.65, y + 0.10, 1.6, 0.35,
             label, size=11, bold=True, color=tc)
    add_text(s, 2.30, y + 0.10, 10.45, 0.35,
             headline, size=13, bold=True, color=INK)
    add_text(s, 0.65, y + 0.50, 12.10, 1.05,
             body, size=10.5, color=INK)
    y += 1.85


# ======================== Slide 9 — questions for discussion =========
s = prs.slides.add_slide(blank)
add_title_block(s, "10", "10 · QUESTIONS FOR DISCUSSION",
                "what I'd like your read on before the campaign re-run")

questions = [
    ("1. ", "Does the 3-stage diagnosis match what you'd intuit physically?",
     "Anionic NF270 + chloride co-ion + increasing cation valence → identifiability should degrade with valence.  "
     "That's what we see.  But is the σ=1 wall on the NaCl/CaCl₂ sheets really a regime problem, or is it the "
     "0.3 % cF scaling forcing the optimizer there?"),
    ("2. ", "MC3 SNaCl and MC5 S2NaCl — two different paths to well-identified",
     "MC3 SNaCl probed σ-sensitive conditions and got interior σ = 0.45.  MC5 S2NaCl ran diluting so σ was "
     "irrelevant — same Lp identifiability for free.  What do the run-sheet metadata show that's different from "
     "the other 4 NaCl sheets?  Pressure span, feed concentration, run duration?"),
    ("2b.", "MC4 SNaCl — diluting but poor fit (WSSE = 9536)",
     "MC4 SNaCl looks like a diluting-regime run (σ = 0 wall), but its WSSE is ~10× worse than MC5 S2NaCl.  "
     "Same regime, same salt, same membrane — different fit quality.  Likely a data-quality issue (fouling, leak, "
     "calibration drift) rather than a model issue.  Worth a closer look."),
    ("3. ", "LaCl₃ ion-pairing fix",
     "Stage 3 confirms the v10-deck slide 12 'physics gap' diagnosis.  Should we (a) parameterize σ(cF) explicitly, "
     "(b) add an ion-pairing equilibrium block with LaCl²⁺ / LaCl⁺ partition, or (c) fit only the diluting-regime "
     "portion?  Need at least 3-5 more LaCl₃ sheets either way."),
    ("4. ", "Priority 1 data — pure-water L_p calibration",
     "I think this is the single highest-leverage experiment.  Is it feasible to add to every NF270 run going "
     "forward?  Roughly:  10 minutes of DI water at the operating ΔP before and after the salt run."),
    ("5. ", "Per-salt B bounds in the optimizer",
     "We've landed NaCl [0, 30], CaCl₂ [0, 15], LaCl₃ [0, 10] in the library (Architecture.md §17).  "
     "Anything you'd revise before re-running the full campaign?"),
]
y = 1.50
for n, title, body in questions:
    add_text(s, 0.45, y, 0.5, 0.35,  n,
             size=14, bold=True, color=ACCENT)
    add_text(s, 0.85, y, 12.0, 0.35,  title,
             size=12.5, bold=True, color=INK)
    add_text(s, 0.85, y + 0.32, 12.0, 0.75,  body,
             size=10, color=GREY)
    y += 1.10


# ======================== Slide 12 — Cross-sheet small-multiples =====
# All 11 sheets' σ × Lp contours, arranged 4 columns × 3 rows, stage color tabs
s = prs.slides.add_slide(blank)
add_title_block(s, "11",
                "11 · ALL 11 SHEETS — σ × Lp small-multiples",
                "every sheet's contour, color-tabbed by stage — visual evidence for the categorization")

# 11 sheets ordered: Stage 1 first, then Stage 2 NaCl, then Stage 2 CaCl₂, then Stage 3
sm_layout = [
    ("MC3.07.22.24_SNaCl",   1, "MC3 SNaCl"),
    ("MC5.07.23.24_S2NaCl",  1, "MC5 S2NaCl"),
    ("MC2.05.07.24_NaCl",    2, "MC2 NaCl"),
    ("MC4.07.11.24_SNaCl",   2, "MC4 SNaCl ⚠"),
    ("MC5.07.23.24_NaCl",    2, "MC5 NaCl"),
    ("MC5.07.23.24_SNaCl",   2, "MC5 SNaCl"),
    ("MC2.05.07.24_CaCl2",   2, "MC2 CaCl₂"),
    ("MC3.07.11.24_SCaCl2",  2, "MC3 SCaCl₂"),
    ("MC3.07.12.24_S2CaCl2", 2, "MC3 S2CaCl₂"),
    ("MC2.05.21.24_LaCl3",   3, "MC2 LaCl₃"),
    ("MC4.07.11.24_SLaCl3",  3, "MC4 SLaCl₃"),
]

cols, rows = 4, 3
img_w = 3.10
img_h = 1.75
gap = 0.05
start_x = 0.45
start_y = 1.60

for idx, (sid, stage, label) in enumerate(sm_layout):
    r, c = divmod(idx, cols)
    x = start_x + c * (img_w + gap)
    y = start_y + r * (img_h + 0.35)
    png_path = FULL / sid / "objcontour-x_sigma-y_Lp.png"
    if not png_path.exists():
        continue
    # Stage-colored border via background rectangle
    if stage == 1:
        bcol = STAGE1_TEXT
    elif stage == 2:
        bcol = STAGE2_TEXT
    else:
        bcol = STAGE3_TEXT
    border = s.shapes.add_shape(MSO_SHAPE.RECTANGLE,
                                 Inches(x - 0.03), Inches(y - 0.03),
                                 Inches(img_w + 0.06), Inches(img_h + 0.32))
    border.fill.solid(); border.fill.fore_color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
    border.line.color.rgb = bcol; border.line.width = Pt(1.5)
    border.shadow.inherit = False
    s.shapes.add_picture(str(png_path),
                          Inches(x), Inches(y),
                          width=Inches(img_w), height=Inches(img_h))
    # Label
    add_text(s, x, y + img_h + 0.02, img_w, 0.25,
             label, size=10, bold=True, color=bcol, align=PP_ALIGN.CENTER)


# ======================== Slide 13 — appendix: NaCl B×Lp + B×σ =======
for label, title, sub, png in [
    ("A1", "A1 · APPENDIX — MC2 NaCl B × Lp contour",
     "σ fixed at warm-start (σ=1) · B × Lp swept · supports the 'σ=1 wall hides B-Lp trade-off' claim",
     ART / "matlab_ports_5ch_allpairs" / "MC2.05.07.24_NaCl" / "objcontour-x_B-y_Lp.png"),
    ("A2", "A2 · APPENDIX — MC2 NaCl B × σ contour",
     "Lp fixed at warm-start (Lp=9.95) · B × σ swept · mass argmins at interior σ ≈ 0.32 — interesting",
     ART / "matlab_ports_5ch_allpairs" / "MC2.05.07.24_NaCl" / "objcontour-x_B-y_sigma.png"),
    ("A3", "A3 · APPENDIX — 3D L_p × B × σ slice composite, NaCl",
     "10 σ-slices × 3 channels.  Confirms the σ × Lp basin structure.",
     ART / "contour3d" / "MC2.05.07.24_NaCl" / "MC2.05.07.24_NaCl_grid3d_slices.png"),
    ("A4", "A4 · APPENDIX — 3D L_p × B × σ slice composite, CaCl₂",
     "Basin shape morphs between σ=0 and σ=1 — the side-split on slide 05 made visible across σ.",
     ART / "contour3d" / "MC2.05.07.24_CaCl2" / "MC2.05.07.24_CaCl2_grid3d_slices.png"),
]:
    if not png.exists():
        continue
    s = prs.slides.add_slide(blank)
    add_title_block(s, label, title, sub)
    add_image_centered(s, png, y_top=1.55, max_height=5.70, max_width=12.6)


# ======================== Animation slides — one per sheet, 3 scrub directions ===
# These embed animated GIFs that PowerPoint plays natively when the slide is shown.
# Drop-in animated equivalents of A3/A4 — same data, scrubbing one parameter at a time.
ANIM_SHEETS = []
if ANIM.exists():
    ANIM_SHEETS = sorted(d.name for d in ANIM.iterdir() if d.is_dir())

anim_label_counter = 5  # start at A5
for sid in ANIM_SHEETS:
    for scrub_name, x_var, y_var, scrub_title, scrub_intuition in [
        ("scrub-sigma_BvsLp.gif", "L_p", "B",
         "σ-scrub  ·  L_p × B  contour, animated through σ slices",
         "Watch how the L_p × B basin morphs as σ moves from 0 → 1.  "
         "When channels track each other across σ, the basin is consistent and we trust the fit."),
        ("scrub-Lp_Bvssigma.gif", "σ", "B",
         "L_p-scrub  ·  B × σ  contour, animated through L_p slices",
         "Watch how the B × σ basin morphs as L_p moves from low to high.  "
         "Reveals whether L_p is identifiable independent of (σ, B)."),
        ("scrub-B_Lpvssigma.gif", "σ", "L_p",
         "B-scrub  ·  L_p × σ  contour, animated through B slices",
         "Watch how the L_p × σ basin morphs as B moves from low to high.  "
         "Reveals the σ-B correlation directly."),
    ]:
        gif_path = ANIM / sid / scrub_name
        if not gif_path.exists():
            continue
        label = f"A{anim_label_counter}"
        anim_label_counter += 1

        salt_pretty = "NaCl" if "NaCl" in sid else "CaCl₂" if "CaCl2" in sid else "LaCl₃"
        s = prs.slides.add_slide(blank)
        add_title_block(s, label,
                        f"{label} · APPENDIX — {sid} animated  ·  {scrub_title}",
                        f"{salt_pretty} sheet · animated 3-panel contour (mass, perm conc, ret conc)")
        # Animated GIF embed — PowerPoint plays it in presentation mode
        add_image_centered(s, gif_path, y_top=1.55, max_height=4.40, max_width=12.6)
        # Plain-language caption
        add_takeaway_box(s, 6.10, scrub_intuition,
                         color=ACCENT, height=1.10)


# Save
prs.save(DST)
print(f"WROTE  {DST}")
print(f"size:   {DST.stat().st_size:>11,} bytes")
print(f"slides: {len(prs.slides)}")
