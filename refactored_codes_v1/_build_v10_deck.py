"""Build DATA3_single_salt_analysis_v10.pptx from v9.

Two changes:
  1. Slide 30 ("10 · PROCESS MODEL · CORE EQUATIONS") — bump the stale bound
     'B ∈ (10⁻⁶, 30)' to '(10⁻⁶, 50)' to match the current refactored_ucb_library
     code (line 1683; bumped 2026-05-24 after CaCl₂ fits clipped at 30).
  2. New slide inserted AFTER slide 30, BEFORE slide 31 ("11 · PARAMETER
     ESTIMATION") — '10A · PARAMETER BOUNDS · SOURCES'. The slide carries
     TWO side-by-side views:

        LEFT  | LITERATURE ANCHOR        — Vanýsek D⁰ + Nair 2018 NF270 Pₛ
              | (sparse, multi-ion fit, with adversarial-verification caveats)

        RIGHT | PHYSICS-SUGGESTED STRATEGY — z² Born + Stokes-Einstein scaling
              | with a KCl baseline box on top (the bounds inherited from
              | published DATA1/DATA2) and per-salt B bounds that drop by a
              | factor of ~2 per valence step.

     A single physics paragraph below ties them together (Yaroshchuk 2000,
     Bandini-Vezzani 2003); an amber caveat strip lists the three honest
     caveats; a Source: footer lists every citation.

No other slides touched. Saves to DATA3_single_salt_analysis_v10.pptx.
"""

from pathlib import Path

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE


SRC = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1/DATA3_single_salt_analysis_v9.pptx")
DST = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1/DATA3_single_salt_analysis_v10.pptx")


# ---------- Palette ----------
INK      = RGBColor(0x1A, 0x1A, 0x1A)
GREY     = RGBColor(0x55, 0x55, 0x55)
LIGHT    = RGBColor(0x8A, 0x8A, 0x8A)
ACCENT   = RGBColor(0xC0, 0x39, 0x2B)
HEAD_BG  = RGBColor(0xEC, 0xEC, 0xEC)
ROW_ALT  = RGBColor(0xF7, 0xF7, 0xF7)
CAVEAT   = RGBColor(0xFF, 0xF5, 0xE0)
CAVEAT_B = RGBColor(0xC6, 0x80, 0x00)
KCL_BG   = RGBColor(0xE8, 0xF0, 0xFA)   # soft blue for the KCl baseline box
KCL_B    = RGBColor(0x1F, 0x4E, 0x79)
LEFT_TAG = RGBColor(0x6E, 0x6E, 0x6E)   # neutral grey for "literature" header
RIGHT_TAG= RGBColor(0x2E, 0x7D, 0x32)   # green for "physics-suggested" header


# ---------- Helpers ----------
def add_text(slide, x, y, w, h, text, *, size=11, bold=False, color=INK,
             align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, italic=False, font="Calibri"):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.margin_left = tf.margin_right = Emu(36000)
    tf.margin_top = tf.margin_bottom = Emu(18000)
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    lines = text.split("\n")
    for i, line in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run()
        r.text = line
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.italic = italic
        r.font.name = font
        r.font.color.rgb = color
    return tb


def set_cell(cell, text, *, size=10, bold=False, color=INK, fill=None,
             align=PP_ALIGN.LEFT, font="Calibri", italic=False):
    cell.margin_left = cell.margin_right = Emu(50000)
    cell.margin_top = cell.margin_bottom = Emu(30000)
    cell.vertical_anchor = MSO_ANCHOR.MIDDLE
    if fill is not None:
        cell.fill.solid()
        cell.fill.fore_color.rgb = fill
    tf = cell.text_frame
    tf.word_wrap = True
    tf.text = ""
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run()
        r.text = line
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.italic = italic
        r.font.name = font
        r.font.color.rgb = color


def add_rounded_box(slide, x, y, w, h, fill, border):
    shp = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                 Inches(x), Inches(y), Inches(w), Inches(h))
    shp.fill.solid()
    shp.fill.fore_color.rgb = fill
    shp.line.color.rgb = border
    shp.line.width = Pt(0.75)
    shp.shadow.inherit = False
    return shp


def move_slide(prs, old_idx, new_idx):
    sldIdLst = prs.slides._sldIdLst
    slides = list(sldIdLst)
    sldIdLst.remove(slides[old_idx])
    sldIdLst.insert(new_idx, slides[old_idx])


# ---------- Build ----------
prs = Presentation(str(SRC))
assert abs(prs.slide_width / 914400 - 13.333) < 0.05, f"unexpected slide width"

# (1) Update slide 30 — bump the stale B bound from 30 → 50
slide30 = prs.slides[29]
for shp in slide30.shapes:
    if not shp.has_text_frame:
        continue
    for para in shp.text_frame.paragraphs:
        for run in para.runs:
            if "10⁻⁶, 30" in run.text:
                run.text = run.text.replace("10⁻⁶, 30", "10⁻⁶, 50")
            if "1e-6, 30" in run.text:
                run.text = run.text.replace("1e-6, 30", "1e-6, 50")


# (2) Insert new "10A" slide
blank_layout = prs.slide_layouts[4]
new = prs.slides.add_slide(blank_layout)


# === Title row (matches deck convention) ===
add_text(new, 0.45, 0.30, 0.6, 0.55, "10A",
         size=28, bold=True, color=LIGHT)
add_text(new, 1.15, 0.30, 11.5, 0.55,
         "10A · PARAMETER BOUNDS  ·  SOURCES",
         size=22, bold=True, color=INK)
add_text(new, 1.15, 0.85, 11.5, 0.35,
         "two complementary views  ·  literature anchor (sparse, multi-ion fit)   |   physics-suggested strategy (z² Born + Stokes-Einstein)",
         size=12, italic=True, color=GREY)

# Caption strip
add_text(new, 0.45, 1.25, 12.4, 0.35,
         "L_p and σ are membrane-only — inherited verbatim from DATA1/DATA2 KCl runs.  Only B has a per-salt picture, and that picture has two flavors.",
         size=11.5, color=ACCENT)


# === Section headers above the two tables ===
LEFT_X,  LEFT_W  = 0.45, 6.10
RIGHT_X, RIGHT_W = 6.78, 6.10
HEADER_Y = 1.78

add_text(new, LEFT_X, HEADER_Y, LEFT_W, 0.35,
         "LITERATURE ANCHOR   ·   sparse, multi-ion seawater fit",
         size=12, bold=True, color=LEFT_TAG)
add_text(new, RIGHT_X, HEADER_Y, RIGHT_W, 0.35,
         "PHYSICS-SUGGESTED STRATEGY   ·   z² Born scaling, ~2× per valence step",
         size=12, bold=True, color=RIGHT_TAG)

# WHY rationale strips — one italics sentence per column, explains why these numbers
add_text(new, LEFT_X, HEADER_Y + 0.32, LEFT_W, 0.55,
         "Why these numbers:  Nair et al. (2018) fitted Spiegler-Kedem to a multi-ion seawater dataset on NF270 and reported B per ion (their notation Pₛ).  These are NOT single-salt B values — they are per-ion permeabilities the optimizer settled on with all coupling absorbed into each one.  Use only as order-of-magnitude anchors; the per-ion split itself is unconventional.",
         size=8.5, italic=True, color=LEFT_TAG)

add_text(new, RIGHT_X, HEADER_Y + 0.32, RIGHT_W, 0.55,
         "Why these numbers:  Start from the KCl box (top, inherited from DATA1/DATA2).  Tighten by ~2× per valence step:  D⁰ falls slightly (1.61 → 1.33 → 1.29 × 10⁻⁹),  but the dominant effect is the z² Born self-energy of placing an ion in the low-dielectric pore — partition coefficients scale 1 : 4 : 9 for z = +1 : +2 : +3.  Combined effect: NaCl ≈ KCl, CaCl₂ ≈ ½ KCl, LaCl₃ ≈ ⅓ KCl.",
         size=8.5, italic=True, color=RIGHT_TAG)


# === LEFT TABLE: Literature anchor ===
# Cols: Salt | z | D⁰ at 25 °C [m²/s] | B from Nair 2018 [µm/s]
LT_Y, LT_H = 2.78, 2.50
lt_shape = new.shapes.add_table(4, 4, Inches(LEFT_X), Inches(LT_Y),
                                Inches(LEFT_W), Inches(LT_H))
lt = lt_shape.table
for i, w in enumerate([0.95, 0.50, 1.65, 3.00]):
    lt.columns[i].width = Inches(w)
lt.rows[0].height = Inches(0.55)
for i in range(1, 4):
    lt.rows[i].height = Inches(0.65)

lt_headers = ["Salt", "z", "D⁰ at 25 °C [m²/s]", "B literature anchor (NF270, µm/s)"]
for j, h in enumerate(lt_headers):
    set_cell(lt.cell(0, j), h, size=10.5, bold=True, color=INK,
             fill=HEAD_BG, align=PP_ALIGN.CENTER)

# Note: Nair uses notation Pₛ; we transcribe their Pₛ values as B (same quantity).
lit_rows = [
    ("NaCl",  "+1", "1.61 × 10⁻⁹",
     "B(Cl⁻) ≈ 21,   B(Na⁺) ≈ 1.5\nNair et al. Membranes 8:78 (2018), Pₛ ≡ B"),
    ("CaCl₂", "+2", "1.33 × 10⁻⁹",
     "B(Ca²⁺) ≈ 18.8     ⚠ B(Ca²⁺) > B(Na⁺)\ncontradicts SE+DE; multi-ion fit, not single-salt"),
    ("LaCl₃", "+3", "1.29 × 10⁻⁹",
     "no peer-reviewed NF270 single-salt B located\n→ falls back to physics-suggested (right)"),
]
for i, (salt, z, d0, b) in enumerate(lit_rows, start=1):
    fill = ROW_ALT if i % 2 == 1 else None
    set_cell(lt.cell(i, 0), salt, size=11, bold=True, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(lt.cell(i, 1), z,    size=11, color=GREY, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(lt.cell(i, 2), d0,   size=10.5, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(lt.cell(i, 3), b,    size=9, fill=fill, align=PP_ALIGN.LEFT)


# === RIGHT TABLE: Physics-suggested ===
# First: a KCl baseline box at the top of the right column
KCL_Y, KCL_H = 2.78, 0.55
kcl_box = add_rounded_box(new, RIGHT_X, KCL_Y, RIGHT_W, KCL_H, KCL_BG, KCL_B)
kcl_tf = kcl_box.text_frame
kcl_tf.margin_left = kcl_tf.margin_right = Emu(72000)
kcl_tf.margin_top = kcl_tf.margin_bottom = Emu(28000)
kcl_tf.word_wrap = True
kcl_tf.vertical_anchor = MSO_ANCHOR.MIDDLE
p0 = kcl_tf.paragraphs[0]
p0.alignment = PP_ALIGN.LEFT
for txt, bold, color, size in [
    ("KCl baseline (locked, inherited from DATA1/DATA2):   ", True,  KCL_B, 10.5),
    ("L_p ∈ [0.5, 50] L/(m²·hr·bar),    B ∈ [0, 30] µm/s,    σ ∈ [0, 1]  (physics-fixed)",
                                                              False, INK,   10.5),
]:
    r = p0.add_run(); r.text = txt
    r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = color

# Then: 4-row × 4-col per-salt table below the KCl box (sized to match LT_H total)
PS_Y, PS_H = KCL_Y + KCL_H + 0.08, 2.50 - (KCL_H + 0.08)
ps_shape = new.shapes.add_table(4, 4, Inches(RIGHT_X), Inches(PS_Y),
                                Inches(RIGHT_W), Inches(PS_H))
ps = ps_shape.table
for i, w in enumerate([0.95, 0.50, 1.35, 3.30]):
    ps.columns[i].width = Inches(w)
ps.rows[0].height = Inches(0.42)
for i in range(1, 4):
    ps.rows[i].height = Inches(0.48)

ps_headers = ["Salt", "z", "B bound (µm/s)", "Basis"]
for j, h in enumerate(ps_headers):
    set_cell(ps.cell(0, j), h, size=10.5, bold=True, color=INK,
             fill=HEAD_BG, align=PP_ALIGN.CENTER)

phys_rows = [
    ("NaCl",  "+1", "[0, 30]", "same as KCl;  monovalent, comparable D⁰"),
    ("CaCl₂", "+2", "[0, 15]", "divalent counter-ion;  lower D⁰ + added dielectric exclusion"),
    ("LaCl₃", "+3", "[0, 10]", "trivalent;  strongest dielectric / screening effect"),
]
for i, (salt, z, b, basis) in enumerate(phys_rows, start=1):
    fill = ROW_ALT if i % 2 == 1 else None
    set_cell(ps.cell(i, 0), salt,  size=11, bold=True, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(ps.cell(i, 1), z,     size=11, color=GREY, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(ps.cell(i, 2), b,     size=11, bold=True, color=KCL_B,
             fill=fill, align=PP_ALIGN.CENTER)
    set_cell(ps.cell(i, 3), basis, size=9.5, fill=fill, align=PP_ALIGN.LEFT)


# === Physics bridge (one compact line — the WHY strips above each table already carry
# the per-column explanation; this line just ties them together for the audience) ===
PHYS_Y = 5.40
add_text(new, 0.45, PHYS_Y, 12.4, 0.65,
         "Reading the two columns together:   the literature anchor (left) is what published NF270 fits actually report;  the physics-suggested column (right) is what Stokes-Einstein × z² Born scaling predicts from the KCl baseline.   "
         "For NaCl, both columns converge near ~20–30 µm/s.   For CaCl₂ the LITERATURE number is ANOMALOUS (Pₛ-Ca > Pₛ-Na) — defer to physics.   For LaCl₃ no literature exists — physics is the only defensible bound.",
         size=10, color=INK)


# === Caveat callout ===
CAV_Y = PHYS_Y + 0.75
cav = add_rounded_box(new, 0.45, CAV_Y, 12.4, 0.55, CAVEAT, CAVEAT_B)
cav_tf = cav.text_frame
cav_tf.margin_left = cav_tf.margin_right = Emu(72000)
cav_tf.margin_top = cav_tf.margin_bottom = Emu(28000)
cav_tf.word_wrap = True
cav_tf.vertical_anchor = MSO_ANCHOR.MIDDLE
cav_p = cav_tf.paragraphs[0]; cav_p.alignment = PP_ALIGN.LEFT
for txt, bold, color in [
    ("⚠  Honest caveats:   ", True, CAVEAT_B),
    ("Nair et al. 2018 Pₛ comes from a coupled multi-ion seawater fit — anchors only.  "
     "Lachheb et al. 2025 NF270 NaCl Pₛ was refuted by adversarial verification.  "
     "No peer-reviewed NF270 LaCl₃ single-salt B was located, so the LaCl₃ value relies entirely on the right-hand physics.",
                              False, INK),
]:
    r = cav_p.add_run(); r.text = txt
    r.font.size = Pt(9.5); r.font.bold = bold; r.font.color.rgb = color


# === Source footer ===
SRC_Y = CAV_Y + 0.65
add_text(new, 0.45, SRC_Y, 12.4, 0.40,
         "Source:   KCl baseline (L_p, σ, B) — Kasturi et al., DATA1/DATA2 KCl runs (self-cite).   "
         "D⁰ — Vanýsek, CRC Handbook of Chemistry and Physics;  Rard & Miller, J. Solution Chem. 8, 701 (1979);  Ribeiro et al., Electrochim. Acta (2008).   "
         "NF270 per-ion Pₛ — Nair, Protasova, Strand & Bilstad, Membranes 8(3), 78 (2018).   "
         "z² dielectric exclusion — Yaroshchuk, Adv. Colloid Interface Sci. 85, 193 (2000);  Bandini & Vezzani, Chem. Eng. Sci. 58, 3303 (2003).",
         size=7.5, italic=True, color=LIGHT)


# Move the new slide to position 30 (0-indexed) — after the Process Model slide
move_slide(prs, len(prs.slides) - 1, 30)


# Save
prs.save(str(DST))
print(f"WROTE  {DST}")
print(f"size:   {DST.stat().st_size:>11,} bytes")
print(f"slides: {len(prs.slides)}")
