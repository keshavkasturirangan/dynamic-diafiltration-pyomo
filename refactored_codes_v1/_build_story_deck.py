#!/usr/bin/env python3
"""Build the DATA3 'B(c), transport regime, and multi-experiment identifiability' STORY deck.

Coherent narrative (one idea per slide), assembled from the analyses produced this session.
Figure-aware + idempotent: missing figures are skipped with a placeholder note, so it can be
run repeatedly as the background campaigns finish.  Re-run any time.

Output: refactored_codes_v1/DATA3_STORY_deck.pptx (+ QA pdf via soffice if available)
Usage: python3 _build_story_deck.py
"""
import sys, subprocess
from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

HERE = Path(__file__).resolve().parent
ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270" / "bform_study"
MS = ART / "model_selection"; POOL = ART / "pooling"; PEC = ART / "peclet"
STK = ART / "stack_contours"; PRED = ART / "predictions"
OUT = HERE / "DATA3_STORY_deck.pptx"

NAVY = RGBColor(0x1F, 0x37, 0x5B); GREEN = RGBColor(0x2C, 0xA0, 0x2C); GREY = RGBColor(0x44, 0x44, 0x44)
HEAD = {"NaCl diluting": "MC3.07.22.24_SNaCl", "NaCl concentrating": "MC2.05.07.24_NaCl",
        "CaCl2": "MC2.05.07.24_CaCl2", "LaCl3": "MC2.05.21.24_LaCl3"}


def _slide(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])  # blank


def _txt(slide, x, y, w, h, text, size=18, bold=False, color=GREY, align=PP_ALIGN.LEFT):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)); tf = tb.text_frame
    tf.word_wrap = True
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = line; p.alignment = align
        for r in p.runs:
            r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = color
    return tb


def _pic(slide, img, x, y, w, h, missing="(figure pending — campaign still running)"):
    img = Path(img)
    if img.exists():
        try:
            slide.shapes.add_picture(str(img), Inches(x), Inches(y), Inches(w), Inches(h)); return True
        except Exception:
            pass
    _txt(slide, x, y + h / 2, w, 0.6, missing, size=14, color=RGBColor(0xB0, 0x00, 0x00), align=PP_ALIGN.CENTER)
    return False


def header(slide, kicker, title):
    _txt(slide, 0.5, 0.3, 12.3, 0.4, kicker, size=13, bold=True, color=GREEN)
    _txt(slide, 0.5, 0.7, 12.3, 0.8, title, size=26, bold=True, color=NAVY)


def take(slide, text):
    _txt(slide, 0.5, 6.9, 12.3, 0.7, "→ " + text, size=15, bold=True, color=NAVY)


def build():
    prs = Presentation(); prs.slide_width = Inches(13.333); prs.slide_height = Inches(7.5)

    # 1 — title
    s = _slide(prs)
    _txt(s, 0.7, 2.4, 12, 1.2, "DATA3 NF270 — Choosing B(c), the transport regime,\nand whether experiments inform each other",
         size=34, bold=True, color=NAVY)
    _txt(s, 0.7, 4.2, 12, 0.6, "Single-salt diafiltration: NaCl · CaCl₂ · LaCl₃   |   per-response AIC · Peclet regime · "
         "multi-experiment pooling", size=17, color=GREY)
    _txt(s, 0.7, 6.6, 12, 0.4, "Generated from the session analyses — auto-updates as the contour campaign completes", size=12, color=GREY)

    # 2 — the question
    s = _slide(prs); header(s, "01 · THE QUESTION", "What sets the solute permeability B, and does it depend on concentration?")
    _txt(s, 0.6, 1.8, 12, 3.0,
         "• B is the membrane solute permeability: Jₛ = Jw·B·(c_in − c_H).\n"
         "• DATA1/DATA2 transport theory: B = Dₘ·H/L — a partition (H) and a diffusion length.\n"
         "• Whether B is constant or concentration-dependent is set by the Peclet number Pe = Jw·L/Dₘ:\n"
         "        – Pe → 0 (diffusion-dominated): B = DₘH/L = CONSTANT\n"
         "        – Pe ~ 1+ (convection + diffusion): B = Jw·Σ βᵢ cⁱ (concentration-dependent Taylor form)\n"
         "• Under a constant-B fit, σ and B and Lₚ are coupled and σ rails — the identifiability problem we solve here.",
         size=17)
    take(s, "Three threads: which B(c) form (AIC), why it depends on c (Peclet), and can pooling experiments fix identifiability.")

    # 3 — Peclet regime
    s = _slide(prs); header(s, "02 · WHY B DEPENDS ON c — THE PECLET REGIME", "Regression MSE vs Pe (≙ DATA2 Fig 8): which regime each salt is in")
    _pic(s, PEC / "peclet_mse.png", 0.4, 1.6, 12.5, 4.1)
    take(s, "Diluting NaCl → Pe→0 (diffusion, B≈const). CaCl₂ / concentrating-NaCl / LaCl₃ → Pe≫1 (convection+diffusion, B(c)). "
            "The regime direction is the robust signal.")

    # 4 — AIC B(c) selection
    s = _slide(prs); header(s, "03 · WHICH B(c) ORDER — DATA2 PER-RESPONSE AIC", "Order chosen by the per-channel AIC (eq 39), not raw WSSE3")
    _pic(s, MS / "aic_selection.png", 0.4, 1.6, 12.5, 4.1)
    take(s, "CaCl₂ → quadratic/linear, cubic REJECTED (Δ=8–25). NaCl cubic 'wins' are mass-channel/convergence artifacts. "
            "Use the per-response AIC, not the pooled JSON AIC.")

    # 5 — identifiability baseline
    s = _slide(prs); header(s, "04 · THE IDENTIFIABILITY PROBLEM", "Per-experiment σ rails; the σ–B–Lₚ landscape is ill-conditioned")
    _pic(s, POOL / "stage1_identifiability.png", 0.4, 1.6, 12.5, 4.1)
    take(s, "σ rails in 9/11 single experiments; condition numbers 10⁶–10¹². A single experiment cannot pin σ.")

    # 6 — shared-sigma pooling
    s = _slide(prs); header(s, "05 · POOLING RESOLVES IT — SHARED σ", "One shared σ per salt-regime group costs ~nothing — and comes off the rail")
    _pic(s, POOL / "sharesigma_profile.png", 0.4, 1.6, 12.5, 4.1)
    take(s, "Cost ≈1× (sharing σ is nearly free → σ is a salt property). Diluting-NaCl identifies σ*=0.75, OFF the individual rails (1.0/0.38).")

    # 7 — stacking (experiments inform each other)
    s = _slide(prs); header(s, "06 · EXPERIMENTS INFORM EACH OTHER — STACKED LANDSCAPES",
                            "Summing the objective surfaces = summing Fisher information → a sharper pooled minimum")
    _pic(s, STK / "CaCl2" / "single" / "sigLp" / "stacked.png", 0.4, 1.6, 12.5, 3.0,
         missing="(CaCl2 σ×Lp stacked — pending)")
    _txt(s, 0.5, 4.8, 12.3, 1.6,
         "Per channel across a group: Σ mass → Lₚ, Σ permeate-conc → B, Σ retentate-conc → σ. White ○ = each experiment's "
         "scattered own minimum; red ★ = the pooled minimum (joint-fit seed). The transport model is evaluated at every node "
         "(slice = forward-sim, legacy calc_contour_2d; profile = re-solve B(c)).", size=15)
    take(s, "Where individuals are sloppy, the summed surface localizes — the visual proof that pooling improves identifiability.")

    # 8 — per-B-form contours (CaCl2 example)
    s = _slide(prs); header(s, "07 · PER B(c)-FORM σ×Lₚ MAPS", "Each concentration-function of B reshapes the σ/Lₚ identifiability (CaCl₂)")
    forms = [("single", "constant"), ("poly1", "linear"), ("poly2", "quadratic"),
             ("poly3", "cubic"), ("sat", "saturating"), ("donnan", "Donnan")]
    for i, (f, lab) in enumerate(forms):
        x = 0.4 + (i % 3) * 4.3; y = 1.7 + (i // 3) * 2.55
        _txt(s, x, y - 0.25, 4.2, 0.3, lab, size=12, bold=True, color=NAVY, align=PP_ALIGN.CENTER)
        _pic(s, STK / "CaCl2" / f / "sigLp" / "slice" / "stacked.png", x, y, 4.1, 2.2, missing=f"({lab} pending)")
    take(s, "Comparing forms: constant rails σ; the concentration-dependent forms relax the σ/Lₚ trade-off differently.")

    # 9 — predictions
    s = _slide(prs); header(s, "08 · CONSTANT-B FIT — MEASURED vs PREDICTED", "The first B-form (constant B) at the optimized parameters")
    pcs = []
    for lbl, rid in HEAD.items():
        d = PRED / rid / "single"
        if d.exists():
            pcs += sorted(d.glob("concentration-*.png"))[:1]
    placed = 0
    for lbl, rid in HEAD.items():
        d = PRED / rid / "single"
        figs = sorted(d.glob("concentration-*.png")) if d.exists() else []
        x = 0.4 + (placed % 2) * 6.4; y = 1.7 + (placed // 2) * 2.55
        _txt(s, x, y - 0.25, 6.2, 0.3, lbl, size=12, bold=True, color=NAVY, align=PP_ALIGN.CENTER)
        _pic(s, figs[0] if figs else d / "missing.png", x, y, 6.1, 2.3, missing=f"({lbl} prediction pending)")
        placed += 1
    take(s, "Constant-B tracks NaCl well (flat B); for CaCl₂/LaCl₃ the constant-B residual motivates B(c).")

    # 10 — summary
    s = _slide(prs); header(s, "09 · THE COHERENT STORY", "")
    _txt(s, 0.6, 1.7, 12, 4.6,
         "1.  B is constant when diffusion dominates (Pe→0) and concentration-dependent when convection matters (Pe≫1).\n"
         "      The Peclet analysis separates the salts: diluting-NaCl diffusion-limited; CaCl₂/LaCl₃/conc-NaCl convection.\n\n"
         "2.  Among concentration-functions of B, the DATA2 per-response AIC selects LOW order — quadratic/linear for CaCl₂,\n"
         "      cubic rejected. Raw WSSE3 over-selects; the per-channel likelihood is the honest criterion.\n\n"
         "3.  A single experiment cannot identify σ (it rails; cond 10⁶–10¹²). Pooling experiments of the same salt+regime\n"
         "      shares σ at ~zero cost and resolves the rail (diluting-NaCl σ*=0.75). Stacked objective landscapes show the\n"
         "      experiments informing each other: the summed surface has a sharper minimum than any individual.\n\n"
         "4.  Same NF270 membrane across DATA3 → Lₚ is consistent; the within-group scatter was identifiability noise that\n"
         "      pooling removes. The stacked minimum is the warm-start for the joint multi-experiment fit.",
         size=16)

    prs.save(str(OUT))
    n = len(prs.slides._sldIdLst)
    print(f"saved {OUT}  ({n} slides)", flush=True)
    try:
        subprocess.run(["/opt/homebrew/bin/soffice", "--headless", "--convert-to", "pdf", "--outdir", str(HERE), str(OUT)],
                       timeout=180, check=False)
        print("QA PDF rendered" if (HERE / (OUT.stem + ".pdf")).exists() else "QA PDF not found", flush=True)
    except Exception as e:
        print("soffice skipped:", e, flush=True)
    print("STORY DECK DONE", flush=True)


if __name__ == "__main__":
    build()
