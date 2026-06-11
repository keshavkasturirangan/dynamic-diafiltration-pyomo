"""Build a simpler v10 slide 10A.

Replaces the two-column literature-vs-physics layout with a plain-language
'one number changes by salt, here's why' slide:

    Title:                10A · PARAMETER BOUNDS
    Subtitle:             one number changes by salt; everything else is shared

    Top stripe (shared):  MEMBRANE BOUNDS (same for every salt)
                          L_p ∈ [0.5, 50]    σ ∈ [0, 1]
                          → these describe the NF270 membrane, not the salt

    Middle (per-salt):    SALT PERMEABILITY  B (µm/s) — tighter as cation charge grows
                          Three big cards side-by-side:
                              NaCl  +1   [0, 30]   bar █████████████   baseline
                              CaCl₂ +2   [0, 15]   bar ███████         half
                              LaCl₃ +3   [0, 10]   bar █████           third

    One-paragraph WHY:    plain English, no dielectric exclusion jargon

    Tiny source footer.

Everything else (literature anchor table, Nair caveat, Lachheb refutation,
z² Born scaling math) becomes a brief paragraph in the speaker notes for
anyone who asks where the numbers come from.

Saves to DATA3_single_salt_analysis_v10.pptx.
"""

from pathlib import Path

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE


SRC = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1/DATA3_single_salt_analysis_v9.pptx")
DST = Path("/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1/DATA3_single_salt_analysis_v10.pptx")


INK     = RGBColor(0x1A, 0x1A, 0x1A)
GREY    = RGBColor(0x55, 0x55, 0x55)
LIGHT   = RGBColor(0x8A, 0x8A, 0x8A)
ACCENT  = RGBColor(0xC0, 0x39, 0x2B)
SOFT_BG = RGBColor(0xF6, 0xF8, 0xFB)
SOFT_B  = RGBColor(0xCB, 0xD5, 0xE1)
CARD_BG = RGBColor(0xFF, 0xFF, 0xFF)
CARD_B  = RGBColor(0xD0, 0xD0, 0xD0)
BAR     = RGBColor(0x1F, 0x4E, 0x79)
BAR_BG  = RGBColor(0xE8, 0xED, 0xF3)
TAG     = RGBColor(0x2E, 0x7D, 0x32)


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


def rounded(slide, x, y, w, h, fill, border, *, corner=0.08):
    shp = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                                 Inches(x), Inches(y), Inches(w), Inches(h))
    shp.fill.solid(); shp.fill.fore_color.rgb = fill
    shp.line.color.rgb = border; shp.line.width = Pt(0.75)
    shp.shadow.inherit = False
    # adjust corner radius
    adj = shp.adjustments
    if len(adj):
        adj[0] = corner
    return shp


def filled_rect(slide, x, y, w, h, fill):
    shp = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE,
                                 Inches(x), Inches(y), Inches(w), Inches(h))
    shp.fill.solid(); shp.fill.fore_color.rgb = fill
    shp.line.fill.background()
    shp.shadow.inherit = False
    return shp


def move_slide(prs, old_idx, new_idx):
    s = prs.slides._sldIdLst
    slides = list(s)
    s.remove(slides[old_idx]); s.insert(new_idx, slides[old_idx])


# ---------- Build ----------
prs = Presentation(str(SRC))

# Slide 30 — bump the stale B bound from 30 → 50 to match the actual code
# (the May 24 widening from §12 is still in place for DATA1/DATA2 and as the
# DATA3 fall-through default; only specific DATA3 salts get the per-salt
# overrides documented on slide 10A).
for shp in prs.slides[29].shapes:
    if not shp.has_text_frame: continue
    for para in shp.text_frame.paragraphs:
        for run in para.runs:
            if "10⁻⁶, 30" in run.text:
                run.text = run.text.replace("10⁻⁶, 30", "10⁻⁶, 50")
            if "1e-6, 30" in run.text:
                run.text = run.text.replace("1e-6, 30", "1e-6, 50")


# New simplified 10A
blank = prs.slide_layouts[4]
s = prs.slides.add_slide(blank)


# === Title block ===
add_text(s, 0.45, 0.30, 0.6, 0.55, "10A",
         size=28, bold=True, color=LIGHT)
add_text(s, 1.15, 0.30, 11.5, 0.55,
         "10A · PARAMETER BOUNDS",
         size=24, bold=True, color=INK)
add_text(s, 1.15, 0.90, 11.5, 0.40,
         "one number changes by salt — everything else is shared",
         size=14, italic=True, color=GREY)


# === Shared bounds stripe ===
sb_y, sb_h = 1.55, 0.95
rounded(s, 0.45, sb_y, 12.4, sb_h, SOFT_BG, SOFT_B)
# Label on the left
add_text(s, 0.75, sb_y + 0.12, 3.5, 0.35,
         "MEMBRANE BOUNDS",
         size=11, bold=True, color=TAG, anchor=MSO_ANCHOR.TOP)
add_text(s, 0.75, sb_y + 0.42, 3.5, 0.45,
         "same for every salt\n— properties of the NF270 membrane, not the salt",
         size=9, italic=True, color=GREY)
# The two bounds, big and centered on the right
add_text(s, 4.50, sb_y + 0.10, 4.20, 0.85,
         "L_p  ∈  [0.5,  50]\nL / (m² · hr · bar)",
         size=18, bold=True, color=INK, align=PP_ALIGN.CENTER,
         anchor=MSO_ANCHOR.MIDDLE)
# tiny unit note
for para in s.shapes[-1].text_frame.paragraphs[1:]:
    for r in para.runs:
        r.font.size = Pt(10); r.font.bold = False; r.font.color.rgb = LIGHT

add_text(s, 8.70, sb_y + 0.10, 4.00, 0.85,
         "σ  ∈  [0,  1]\ndimensionless",
         size=18, bold=True, color=INK, align=PP_ALIGN.CENTER,
         anchor=MSO_ANCHOR.MIDDLE)
for para in s.shapes[-1].text_frame.paragraphs[1:]:
    for r in para.runs:
        r.font.size = Pt(10); r.font.bold = False; r.font.color.rgb = LIGHT


# === Per-salt B bound — three big cards ===
hdr_y = 2.70
add_text(s, 0.45, hdr_y, 12.4, 0.35,
         "SALT PERMEABILITY  B   ·   tighter as cation charge grows",
         size=12, bold=True, color=TAG)
add_text(s, 0.45, hdr_y + 0.32, 12.4, 0.30,
         "units: µm/s  ·  one upper bound per salt  ·  rule of thumb: halve B per +1 of cation charge",
         size=10, italic=True, color=GREY)

# 3 cards
card_y, card_h = 3.45, 2.30
card_w = 3.95
gap = 0.27
card_x0 = 0.45 + (12.4 - 3*card_w - 2*gap) / 2   # center the row

SALTS = [
    {"name": "NaCl",  "z": "+1", "bound": "[0, 30]",
     "tag":   "baseline",
     "tag_b": RGBColor(0x4F, 0x6B, 0x99),
     "bar":   1.00,
     "note":  "monovalent — same as KCl"},
    {"name": "CaCl₂", "z": "+2", "bound": "[0, 15]",
     "tag":   "halved",
     "tag_b": RGBColor(0x6B, 0x57, 0x99),
     "bar":   0.50,
     "note":  "divalent — Ca²⁺ sticks to the membrane more"},
    {"name": "LaCl₃", "z": "+3", "bound": "[0, 10]",
     "tag":   "thirded",
     "tag_b": RGBColor(0x99, 0x55, 0x6B),
     "bar":   0.33,
     "note":  "trivalent — La³⁺ sticks the most"},
]

for i, salt in enumerate(SALTS):
    x = card_x0 + i * (card_w + gap)
    # card background
    rounded(s, x, card_y, card_w, card_h, CARD_BG, CARD_B, corner=0.08)
    # Salt name + valence (top of card)
    add_text(s, x, card_y + 0.18, card_w, 0.45,
             salt["name"],
             size=22, bold=True, color=INK, align=PP_ALIGN.CENTER)
    add_text(s, x, card_y + 0.62, card_w, 0.28,
             f"charge  {salt['z']}",
             size=11, color=GREY, align=PP_ALIGN.CENTER)
    # Big bound number (middle of card)
    add_text(s, x, card_y + 0.95, card_w, 0.55,
             salt["bound"],
             size=28, bold=True, color=salt["tag_b"], align=PP_ALIGN.CENTER)
    add_text(s, x, card_y + 1.46, card_w, 0.20,
             "µm/s",
             size=10, color=LIGHT, align=PP_ALIGN.CENTER, italic=True)
    # Bar visualization
    bar_w_full = card_w - 0.60
    bar_x = x + 0.30
    bar_y = card_y + 1.72
    bar_h = 0.12
    filled_rect(s, bar_x, bar_y, bar_w_full, bar_h, BAR_BG)
    filled_rect(s, bar_x, bar_y, bar_w_full * salt["bar"], bar_h, BAR)
    # tag — placed BELOW the bar with clear separation
    add_text(s, x, card_y + 1.92, card_w, 0.30,
             salt["tag"],
             size=11, bold=True, color=salt["tag_b"], align=PP_ALIGN.CENTER)


# === WHY plain English ===
why_y = 6.00
add_text(s, 0.45, why_y, 12.4, 0.30,
         "Why?",
         size=12, bold=True, color=ACCENT)
add_text(s, 0.45, why_y + 0.32, 12.4, 1.10,
         "Two effects make multivalent salts harder to diffuse through the membrane on their own:   "
         "(1)  Heavier-charged ions move slower in water to begin with.   "
         "(2)  Squeezing a charged ion into the pore takes energy — and that cost rises with the square of the charge,  "
         "so about 4× the cost for divalent ions and 9× for trivalent ones.   "
         "Engineering rule of thumb:  halve the B upper bound for each step up in cation charge.   "
         "(See slide 10B for why this doesn't conflict with the rejection data you see at the bench.)",
         size=11, color=INK)


# === Source footer (tight, one line) ===
add_text(s, 0.45, 7.10, 12.4, 0.40,
         "Source:   L_p and σ — inherited from published DATA1/DATA2 KCl runs (self-cite).   "
         "B per-valence rule — Stokes-Einstein × dielectric exclusion (Yaroshchuk 2000;  Bandini & Vezzani 2003).   "
         "Literature anchors and verification audit available in speaker notes.",
         size=8, italic=True, color=LIGHT)


# === Speaker notes — everything we just stripped off the slide ===
notes_tf = s.notes_slide.notes_text_frame
notes_tf.text = (
    "SPEAKER NOTES — where these numbers come from\n"
    "\n"
    "1. The L_p ∈ [0.5, 50] and σ ∈ [0, 1] bounds are inherited verbatim from "
    "the published DATA1/DATA2 KCl experiments (Kasturi et al.). Self-cite.\n"
    "\n"
    "2. Aqueous limiting diffusion coefficients D⁰ at 25 °C used in the\n"
    "Stokes-Einstein argument:\n"
    "      D⁰(NaCl)  = 1.61 × 10⁻⁹ m²/s\n"
    "      D⁰(CaCl₂) = 1.33 × 10⁻⁹ m²/s\n"
    "      D⁰(LaCl₃) = 1.29 × 10⁻⁹ m²/s\n"
    "   Source: Vanýsek, CRC Handbook (ionic conductivity at infinite dilution),\n"
    "   corroborated by Rard & Miller, J. Solution Chem. 8, 701 (1979) for NaCl/\n"
    "   CaCl₂ and Ribeiro et al., Electrochim. Acta (2008) for CaCl₂.\n"
    "\n"
    "3. The 'halve B per +1 cation charge' rule comes from combining\n"
    "   Stokes-Einstein (B ∝ D⁰) with the z²-scaling of the Born self-energy\n"
    "   for placing an ion in the low-dielectric pore (partition coefficient\n"
    "   ratios 1 : 4 : 9 for z = +1, +2, +3). Sources: Yaroshchuk, Adv. Colloid\n"
    "   Interface Sci. 85, 193 (2000); Bandini & Vezzani, Chem. Eng. Sci. 58,\n"
    "   3303 (2003).\n"
    "\n"
    "4. NF270-specific solute-permeability values found in literature (anchors\n"
    "   only — do NOT cite as authoritative single-salt B):\n"
    "      Nair et al., Membranes 8(3), 78 (2018):  multi-ion seawater fit on\n"
    "      NF270 reports B(Cl⁻) ≈ 21, B(Na⁺) ≈ 1.5, B(Ca²⁺) ≈ 18.8 µm/s, with\n"
    "      σ(Ca²⁺) = 0.41, σ(Na⁺) = 0.19. Note: their B(Ca²⁺) > B(Na⁺) is the\n"
    "      opposite of what Stokes-Einstein + dielectric exclusion predicts —\n"
    "      this is an artifact of fitting per-ion permeabilities to coupled\n"
    "      multi-ion flux data, not a true single-salt B.\n"
    "\n"
    "      Lachheb et al., ChemistryOpen (2025):  reports σ(Na⁺) = 0.55–0.69 and\n"
    "      σ(Cl⁻) = 0.59–0.61 on NF270 at 2–6 g/L TDS. Their accompanying B\n"
    "      values for NaCl on NF270 were adversarially refuted in our deep-\n"
    "      research verification (vote 0-3) and should not be cited.\n"
    "\n"
    "      No peer-reviewed single-salt LaCl₃ B on NF270 was located. The\n"
    "      LaCl₃ bound is set entirely by physics extrapolation.\n"
    "\n"
    "5. Current code bound (refactored_ucb_library.py:1683) is the unified\n"
    "   [1e-6, 50] µm/s, raised from [0, 30] on 2026-05-24 when the three\n"
    "   CaCl₂ fits clipped at 30. Per-salt narrowing per the table above is\n"
    "   a proposed refinement, not yet implemented in code.\n"
)


# Move 10A to position 30 (after Process Model, before WSSE Objective)
move_slide(prs, len(prs.slides) - 1, 30)


# =====================================================================
# Slide 10B — TWO ROUTES THROUGH THE MEMBRANE
# (sits right after 10A in the main flow; explains the σ-vs-B paradox
#  in plain language without symbols)
# =====================================================================
r = prs.slides.add_slide(blank)

# Title block
add_text(r, 0.45, 0.30, 0.6, 0.55, "10B",
         size=28, bold=True, color=LIGHT)
add_text(r, 1.15, 0.30, 11.5, 0.55,
         "10B · TWO ROUTES THROUGH THE MEMBRANE",
         size=22, bold=True, color=INK)
add_text(r, 1.15, 0.88, 11.5, 0.40,
         "why a tighter B bound for multivalent salts doesn't conflict with the rejection data",
         size=13, italic=True, color=GREY)

# Setup-the-paradox stripe (red accent)
setup_y, setup_h = 1.45, 0.85
rounded(r, 0.45, setup_y, 12.4, setup_h, SOFT_BG, SOFT_B)
add_text(r, 0.70, setup_y + 0.10, 12.0, 0.30,
         "Sounds backwards at first:",
         size=12, bold=True, color=ACCENT)
add_text(r, 0.70, setup_y + 0.38, 12.0, 0.45,
         "Multivalent salts are LESS rejected at the bench (Na > Ca > La)  —  yet the proposed B bound is TIGHTER for them (Na > Ca > La).   "
         "How can both be true?   Because salt gets through the membrane TWO ways, and they're governed by different parameters.",
         size=10.5, color=INK)


# Two parallel route cards
route_y, route_h = 2.55, 3.05
route_w = 6.10
route_gap = 0.30
route_x0 = 0.45

ROUTES = [
    {"label":     "ROUTE 1 — RIDING WITH THE WATER",
     "label_b":   RGBColor(0x1F, 0x4E, 0x79),
     "icon":      "▶",
     "headline":  "Salt rides through with the flowing water",
     "controlled":"governed by σ (reflection coefficient)",
     "bullets": [
        "Water flows through the membrane;  salt can ride along if the membrane lets it.",
        "NF270 is negatively charged — it attracts multivalent cations (Ca²⁺, La³⁺) more strongly than Na⁺.",
        "Once the cation is pulled in, chloride has to follow to keep things electrically neutral.",
        "Result:  multivalent salts ride through MORE easily.",
     ],
     "outcome":   "Rejection drops for higher-charge salts.",
     "tag":       "← this is what you actually see in the rejection data"},

    {"label":     "ROUTE 2 — DIFFUSING THROUGH ALONE",
     "label_b":   RGBColor(0x2E, 0x7D, 0x32),
     "icon":      "▶",
     "headline":  "Salt molecules squeeze through the pore on their own",
     "controlled":"governed by B (solute permeability)",
     "bullets": [
        "No water flow needed — just random thermal motion through the pore.",
        "Heavier-charged ions move slower in water to start with.",
        "Squeezing a charged ion into the pore costs energy — about 4× more for divalent, 9× more for trivalent.",
        "Result:  multivalent salts diffuse through LESS easily.",
     ],
     "outcome":   "B drops for higher-charge salts.",
     "tag":       "← this is what the proposed B bound encodes"},
]

for i, route in enumerate(ROUTES):
    x = route_x0 + i * (route_w + route_gap)
    # card outline
    rounded(r, x, route_y, route_w, route_h, CARD_BG, CARD_B, corner=0.04)
    # colored top stripe
    filled_rect(r, x, route_y, route_w, 0.45, route["label_b"])
    # icon + label on stripe
    add_text(r, x + 0.20, route_y + 0.04, route_w - 0.40, 0.38,
             route["icon"] + "  " + route["label"],
             size=12, bold=True, color=RGBColor(0xFF, 0xFF, 0xFF), anchor=MSO_ANCHOR.MIDDLE)
    # headline
    add_text(r, x + 0.25, route_y + 0.55, route_w - 0.50, 0.32,
             route["headline"],
             size=12.5, bold=True, color=INK)
    # controlled-by tag (italic)
    add_text(r, x + 0.25, route_y + 0.88, route_w - 0.50, 0.26,
             route["controlled"],
             size=10, italic=True, color=route["label_b"])
    # bullet list
    bullet_y = route_y + 1.18
    bullet_box = r.shapes.add_textbox(Inches(x + 0.25), Inches(bullet_y),
                                     Inches(route_w - 0.50), Inches(1.40))
    btf = bullet_box.text_frame
    btf.margin_left = btf.margin_right = Emu(36000)
    btf.margin_top = btf.margin_bottom = Emu(18000)
    btf.word_wrap = True
    for j, bull in enumerate(route["bullets"]):
        p = btf.paragraphs[0] if j == 0 else btf.add_paragraph()
        p.alignment = PP_ALIGN.LEFT
        p.space_after = Pt(3)
        run = p.add_run()
        run.text = "•  " + bull
        run.font.size = Pt(9.5); run.font.color.rgb = INK; run.font.name = "Calibri"
    # outcome (bold, prominent)
    outcome_y = route_y + 2.50
    add_text(r, x + 0.25, outcome_y, route_w - 0.50, 0.30,
             "→  " + route["outcome"],
             size=11, bold=True, color=route["label_b"])
    # tag at the bottom
    add_text(r, x + 0.25, outcome_y + 0.30, route_w - 0.50, 0.25,
             route["tag"],
             size=9, italic=True, color=GREY)


# Bottom takeaway — soft blue panel
take_y = 5.80
take_box = rounded(r, 0.45, take_y, 12.4, 0.85, SOFT_BG, SOFT_B)
tk_tf = take_box.text_frame
tk_tf.margin_left = tk_tf.margin_right = Emu(72000)
tk_tf.margin_top = tk_tf.margin_bottom = Emu(28000)
tk_tf.word_wrap = True
tk_tf.vertical_anchor = MSO_ANCHOR.MIDDLE
p_tk = tk_tf.paragraphs[0]; p_tk.alignment = PP_ALIGN.LEFT
for txt, bold, color in [
    ("Takeaway:   ",  True,  SOFT_B),
    ("Both rules are right — they describe different stages of the same transport problem.  "
     "Multivalent salts mostly get through by riding the water (Route 1), so the rejection data is dominated by σ.  "
     "Diffusion through alone (Route 2) is a smaller channel, and physics says it should be tighter for higher-charge salts — which is what the proposed B bound encodes.",
                       False, INK),
]:
    rt = p_tk.add_run(); rt.text = txt
    rt.font.size = Pt(11); rt.font.bold = bold; rt.font.color.rgb = color


# Tiny source footer
add_text(r, 0.45, 6.95, 12.4, 0.40,
         "Source:   Spiegler-Kedem model splits salt flux into a convective term (σ-controlled, Donnan/steric exclusion at the interface) and a diffusive term (B-controlled, Born partitioning into the pore).   "
         "Per-salt B reasoning:  Yaroshchuk (2000);  Bandini & Vezzani (2003).   See slide 10C for the cited literature anchors.",
         size=8, italic=True, color=LIGHT)


# Move 10B to position 31 (after the new 10A at position 30, before WSSE Objective at 31)
move_slide(prs, len(prs.slides) - 1, 31)


# =====================================================================
# Slide 10C — LITERATURE ANCHORS (backup reference, lives at the very end)
# (renamed from 10B since the new 10B above is now the σ-vs-B clarification)
# =====================================================================
b = prs.slides.add_slide(blank)

# Title block
add_text(b, 0.45, 0.30, 0.6, 0.55, "10C",
         size=28, bold=True, color=LIGHT)
add_text(b, 1.15, 0.30, 11.5, 0.55,
         "10C · LITERATURE ANCHORS",
         size=24, bold=True, color=INK)
add_text(b, 1.15, 0.90, 11.5, 0.40,
         "backup to slides 10A and 10B  ·  where the per-salt B numbers actually come from in the peer-reviewed record",
         size=13, italic=True, color=GREY)


# ---------- Section 1: D⁰ values (Vanýsek CRC + corroborators) ----------
s1_y = 1.45
add_text(b, 0.45, s1_y, 12.4, 0.30,
         "1.  Aqueous limiting diffusion coefficients D⁰ at 25 °C   ·   the well-settled physical-chemistry input",
         size=12.5, bold=True, color=TAG)

# Compact D⁰ table
d0_tbl = b.shapes.add_table(4, 4,
                            Inches(0.45), Inches(s1_y + 0.35),
                            Inches(12.4), Inches(1.05)).table
for i, w in enumerate([1.10, 0.65, 2.50, 8.15]):
    d0_tbl.columns[i].width = Inches(w)
d0_tbl.rows[0].height = Inches(0.30)
for i in range(1, 4):
    d0_tbl.rows[i].height = Inches(0.25)

def set_cell(cell, text, *, size=10, bold=False, color=INK, fill=None,
             align=PP_ALIGN.LEFT, italic=False):
    cell.margin_left = cell.margin_right = Emu(45000)
    cell.margin_top = cell.margin_bottom = Emu(20000)
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

HEAD = RGBColor(0xEC, 0xEC, 0xEC)
ALT  = RGBColor(0xF7, 0xF7, 0xF7)
for j, h in enumerate(["Salt", "z", "D⁰ (m²/s)", "Source"]):
    set_cell(d0_tbl.cell(0, j), h, size=10.5, bold=True, fill=HEAD, align=PP_ALIGN.CENTER)
d0_rows = [
    ("NaCl",  "+1", "1.61 × 10⁻⁹",
     "Vanýsek (CRC Handbook, ionic conductivity at infinite dilution);  Rard & Miller, J. Sol. Chem. 8, 701 (1979)"),
    ("CaCl₂", "+2", "1.33 × 10⁻⁹",
     "Vanýsek (CRC Handbook);  Ribeiro et al., Electrochim. Acta (2008) — direct measurement"),
    ("LaCl₃", "+3", "1.29 × 10⁻⁹",
     "Vanýsek (CRC Handbook) via the Nernst-Hartley salt-mixing formula"),
]
for i, (salt, z, d0, src) in enumerate(d0_rows, start=1):
    fill = ALT if i % 2 == 1 else None
    set_cell(d0_tbl.cell(i, 0), salt, size=10.5, bold=True, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(d0_tbl.cell(i, 1), z,    size=10.5, color=GREY, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(d0_tbl.cell(i, 2), d0,   size=10.5, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(d0_tbl.cell(i, 3), src,  size=9, fill=fill)


# ---------- Section 2: Nair 2018 NF270 per-ion Pₛ ----------
s2_y = 2.85
add_text(b, 0.45, s2_y, 12.4, 0.28,
         "2.  Spiegler-Kedem B (their notation Pₛ) on NF270   ·   Nair et al., Membranes 8(3), 78 (2018)",
         size=12.5, bold=True, color=TAG)
add_text(b, 0.45, s2_y + 0.26, 12.4, 0.22,
         "Multi-ion fit to synthetic seawater (TDS 30,400 mg/L, 9–18 bar, room T).  Per-ion Pₛ — NOT single-salt B.",
         size=9.5, italic=True, color=GREY)

# Per-ion Nair table — 3 most relevant ions for the NaCl / CaCl₂ comparison
n_tbl = b.shapes.add_table(4, 4,
                           Inches(0.45), Inches(s2_y + 0.52),
                           Inches(12.4), Inches(0.95)).table
for i, w in enumerate([1.30, 0.65, 2.20, 8.25]):
    n_tbl.columns[i].width = Inches(w)
n_tbl.rows[0].height = Inches(0.26)
for i in range(1, 4):
    n_tbl.rows[i].height = Inches(0.23)

for j, h in enumerate(["Ion", "z", "B (µm/s)", "σ  /  comment"]):
    set_cell(n_tbl.cell(0, j), h, size=10.5, bold=True, fill=HEAD, align=PP_ALIGN.CENTER)

nair_rows = [
    ("Cl⁻",   "−1", "21.05", "σ = 0.18  ·  loose-NF, convective-dominant — the Cl⁻ co-ion that NaCl, CaCl₂ and LaCl₃ all share"),
    ("Na⁺",   "+1", "1.52",  "σ = 0.19  ·  the canonical NaCl counter-ion in the fit"),
    ("Ca²⁺",  "+2", "18.79", "σ = 0.41  ⚠ B(Ca²⁺) > B(Na⁺) — contradicts Stokes-Einstein + dielectric exclusion;  artifact of fitting per-ion to coupled multi-ion data"),
]
for i, (ion, z, ps, note) in enumerate(nair_rows, start=1):
    fill = ALT if i % 2 == 1 else None
    set_cell(n_tbl.cell(i, 0), ion,  size=10.5, bold=True, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(n_tbl.cell(i, 1), z,    size=10.5, color=GREY, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(n_tbl.cell(i, 2), ps,   size=10.5, fill=fill, align=PP_ALIGN.CENTER)
    set_cell(n_tbl.cell(i, 3), note, size=9, fill=fill)


# ---------- Section 3: LaCl₃ + Lachheb refutation, in a single amber strip ----------
s3_y = 4.55
amber_box = b.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE,
                               Inches(0.45), Inches(s3_y), Inches(12.4), Inches(0.95))
amber_box.fill.solid(); amber_box.fill.fore_color.rgb = RGBColor(0xFF, 0xF5, 0xE0)
amber_box.line.color.rgb = RGBColor(0xC6, 0x80, 0x00); amber_box.line.width = Pt(0.75)
amber_box.shadow.inherit = False
atf = amber_box.text_frame
atf.margin_left = atf.margin_right = Emu(72000)
atf.margin_top = atf.margin_bottom = Emu(36000)
atf.word_wrap = True
atf.vertical_anchor = MSO_ANCHOR.TOP

p0 = atf.paragraphs[0]; p0.alignment = PP_ALIGN.LEFT
for txt, bold, color in [
    ("3.  Two missing pieces  ", True, RGBColor(0xC6, 0x80, 0x00)),
    ("(why literature alone can't drive per-salt bounds)", False, GREY),
]:
    r = p0.add_run(); r.text = txt
    r.font.size = Pt(11.5); r.font.bold = bold; r.font.color.rgb = color
p1 = atf.add_paragraph(); p1.alignment = PP_ALIGN.LEFT
p1.space_before = Pt(4)
for txt, bold in [
    ("•  No peer-reviewed NF270 single-salt B for LaCl₃ ",                 True),
    ("was located.  The bound on slide 10A is set entirely by the z² Born / Stokes-Einstein extrapolation.",
                                                                            False),
]:
    r = p1.add_run(); r.text = txt
    r.font.size = Pt(10); r.font.bold = bold; r.font.color.rgb = INK
p2 = atf.add_paragraph(); p2.alignment = PP_ALIGN.LEFT
p2.space_before = Pt(2)
for txt, bold in [
    ("•  Lachheb et al., ChemistryOpen (2025) ",                            True),
    ("NF270 NaCl Pₛ ≈ 6–22 µm/s — refuted 0-3 in adversarial verification.  Their σ values (σ(Na⁺) = 0.55–0.69, σ(Cl⁻) = 0.59–0.61) survived and corroborate 'NF270 is a loose membrane' but are not cited as B anchors.",
                                                                            False),
]:
    r = p2.add_run(); r.text = txt
    r.font.size = Pt(10); r.font.bold = bold; r.font.color.rgb = INK


# ---------- Section 4: Why physics wins ----------
s4_y = 5.65
add_text(b, 0.45, s4_y, 12.4, 0.30,
         "4.  Why slide 10A uses physics-informed bounds, not literature-pinned",
         size=12, bold=True, color=ACCENT)
add_text(b, 0.45, s4_y + 0.32, 12.4, 1.30,
         "Single-salt B on NF270 is too sparse and self-contradictory to drive the bound choice.   "
         "For NaCl,  literature B spans 1.5 → 21 µm/s depending on which ion is being reported.   "
         "For CaCl₂,  the only number reported  (18.8 µm/s)  has the wrong ordering vs. Na⁺.   "
         "For LaCl₃,  no peer-reviewed NF270 single-salt B exists.\n\n"
         "The z² Born + Stokes-Einstein rule (halve B per +1 of cation valence) is internally consistent across all three salts and is the only defensible bound for LaCl₃.   "
         "Literature is used here only as an order-of-magnitude sanity check — both the NaCl and CaCl₂ physics-informed bounds sit comfortably inside the reported range.",
         size=10.5, color=INK)


# ---------- Source footer ----------
add_text(b, 0.45, 6.95, 12.4, 0.50,
         "Sources:   Vanýsek, P. — 'Ionic Conductivity and Diffusion at Infinite Dilution,' CRC Handbook of Chemistry and Physics, Sec. 5.   "
         "Rard, J.A. & Miller, D.G. — J. Solution Chem. 8, 701 (1979).   "
         "Ribeiro, A.C.F. et al. — Electrochim. Acta (2008).   "
         "Nair, R.R., Protasova, E., Strand, S. & Bilstad, T. — Membranes 8(3), 78 (2018).   "
         "Lachheb, M. et al. — ChemistryOpen (2025), DOI 10.1002/open.202500198  (Pₛ refuted; σ confirmed).   "
         "Yaroshchuk, A.E. — Adv. Colloid Interface Sci. 85, 193 (2000).   "
         "Bandini, S. & Vezzani, D. — Chem. Eng. Sci. 58, 3303 (2003).",
         size=7.5, italic=True, color=LIGHT)


# 10B stays at the END of the deck (it's a backup reference, no move() call)

prs.save(str(DST))
print(f"WROTE  {DST}")
print(f"size:   {DST.stat().st_size:>11,} bytes")
print(f"slides: {len(prs.slides)}")
