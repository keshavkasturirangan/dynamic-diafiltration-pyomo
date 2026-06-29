#!/usr/bin/env python3
"""Focused presentation deck: FITTING of the 4 representative DATA3 NF270 experiments.
One slide per experiment — measured vs predicted (mass + concentration) at the optimized
parameters, with the parameter box and a one-line takeaway.  Figure-aware + idempotent.

Output: refactored_codes_v1/DATA3_FITTING_4experiments.pptx (+ QA pdf)
Usage: python3 _build_fitting_deck.py
"""
import json, subprocess
from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

HERE = Path(__file__).resolve().parent
ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270" / "bform_study"
PRED = ART / "predictions"; MS = ART / "model_selection"; PEC = ART / "peclet"; PCON = ART / "profile_contours"
OUT = HERE / "DATA3_FITTING_4experiments.pptx"
NAVY = RGBColor(0x1F, 0x37, 0x5B); GREEN = RGBColor(0x2C, 0xA0, 0x2C); GREY = RGBColor(0x44, 0x44, 0x44); RED = RGBColor(0xB0, 0x00, 0x00)

# label -> (rid, fit form tag, form-arg-for-predictions, B-form description, takeaway)
EXP = [
    ("NaCl — diluting regime  (MC3 · 07.22.24)", "MC3.07.22.24_SNaCl", "single", "constant B",
     "B is flat over the range — constant B already tracks mass and both concentrations (WSSE₃=42)."),
    ("NaCl — concentrating regime  (MC2 · 05.07.24)", "MC2.05.07.24_NaCl", "single", "constant B",
     "Constant B is adequate for NaCl (WSSE₃=114); B is essentially concentration-independent."),
    ("CaCl₂  (MC2 · 05.07.24)", "MC2.05.07.24_CaCl2", "poly2", "quadratic  B = Jw(β₀+β₁c+β₂c²)",
     "B(c) is required: quadratic tracks the permeate (direct solute signal) and cuts WSSE₃ 220→147 "
     "(~33%); the retentate is still over-predicted at high c → σ(c)/polarization beyond B(c)."),
    ("LaCl₃  (MC2 · 05.21.24)", "MC2.05.21.24_LaCl3", "poly2", "quadratic  B = Jw(β₀+β₁c+β₂c²)",
     "Model-limited: no B(c) form captures the trivalent — mass under- and retentate over-predicted "
     "(WSSE₃~1190–1244); the residual points past B(c) alone (σ(c), ion-pairing)."),
]


def _slide(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


def _txt(slide, x, y, w, h, text, size=18, bold=False, color=GREY, align=PP_ALIGN.LEFT):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)); tf = tb.text_frame; tf.word_wrap = True
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = line; p.alignment = align
        for r in p.runs:
            r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = color
    return tb


def _pic(slide, globdir, pat, x, y, w, h, fallback=None):
    figs = sorted(Path(globdir).glob(pat)) if Path(globdir).exists() else []
    if not figs and fallback:
        figs = sorted(Path(fallback).glob(pat)) if Path(fallback).exists() else []
    if figs:
        try:
            slide.shapes.add_picture(str(figs[0]), Inches(x), Inches(y), Inches(w), Inches(h)); return True
        except Exception:
            pass
    _txt(slide, x, y + h / 2, w, 0.5, "(fit figure pending)", size=13, color=RED, align=PP_ALIGN.CENTER)
    return False


def _picpath(slide, path, x, y, w, h):
    path = Path(path)
    if path.exists():
        try:
            slide.shapes.add_picture(str(path), Inches(x), Inches(y), Inches(w), Inches(h)); return True
        except Exception:
            pass
    _txt(slide, x, y + h / 2, w, 0.5, "(figure pending — campaign still running)", size=13, color=RED, align=PP_ALIGN.CENTER)
    return False


def contour_fig(rid, tag):
    """Resolve the supporting σ×Lp objective contour for a fit (fall back across form/plane)."""
    cands = [(tag, "sigLp"), ("single", "sigLp"), ("poly1", "sigLp"), ("poly2", "sigLp"),
             ("single", "BLp"), ("single", "sigB")]
    plab = {"sigLp": "σ×Lₚ", "BLp": "B×Lₚ", "sigB": "σ×B"}
    flab = {"single": "constant B", "poly1": "linear B(c)", "poly2": "quadratic B(c)", "sat": "saturating B(c)"}
    for f, pl in cands:
        p = PCON / rid / f / pl / "contour.png"
        if p.exists():
            note = "" if (f == tag and pl == "sigLp") else "  (nearest available form/plane)"
            return p, f"{plab[pl]} objective contour · {flab.get(f, f)}{note}"
    return PCON / rid / "missing.png", "(contour pending)"


def param_box(rid, tag):
    p = json.loads((ART / rid / f"result_{tag}.json").read_text())
    par = p["parameters"]; w = p["WSSE3"]
    lines = [f"Lₚ = {par['Lp']:.2f} L m⁻² h⁻¹ bar⁻¹", f"σ = {par['sigma']:.2f}"]
    if tag == "single":
        lines.insert(1, f"B = {par.get('B'):.2f} µm s⁻¹")
    else:
        lines.insert(1, f"β₀={par.get('beta_0',0):.3f}  β₁={par.get('beta_1',0):.3f}  β₂={par.get('beta_2',0):.4f}")
    lines.append(f"WSSE₃ = {w:.0f}")
    return "\n".join(lines)


def build():
    prs = Presentation(); prs.slide_width = Inches(13.333); prs.slide_height = Inches(7.5)
    # title
    s = _slide(prs)
    _txt(s, 0.7, 2.6, 12, 1.1, "Fitting the four representative DATA3 NF270 experiments",
         size=32, bold=True, color=NAVY)
    _txt(s, 0.7, 3.9, 12, 0.6, "Measured vs predicted at the optimized parameters — NaCl (diluting & concentrating), CaCl₂, LaCl₃",
         size=17, color=GREY)
    # per-experiment slides
    for lbl, rid, tag, desc, take in EXP:
        s = _slide(prs)
        _txt(s, 0.5, 0.3, 12.3, 0.5, lbl, size=24, bold=True, color=NAVY)
        _txt(s, 0.5, 0.85, 12.3, 0.35, f"fit form: {desc}", size=14, bold=True, color=GREEN)
        d = PRED / rid / tag
        dfall = PRED / rid / "single"
        _txt(s, 0.6, 1.3, 6, 0.3, "Mass — measured (points) vs predicted (line)", size=12, color=GREY)
        _pic(s, d, "mass-*.png", 0.4, 1.6, 6.3, 3.5, fallback=dfall)
        _txt(s, 6.9, 1.3, 6, 0.3, "Concentration — retentate / permeate", size=12, color=GREY)
        _pic(s, d, "concentration-*.png", 6.7, 1.6, 6.3, 3.5, fallback=dfall)
        # parameter box
        try:
            pb = param_box(rid, tag)
        except Exception:
            pb = param_box(rid, "single")
        _txt(s, 0.6, 5.35, 5.2, 1.7, "Optimized parameters\n" + pb, size=14, bold=False, color=NAVY)
        _txt(s, 6.2, 5.35, 6.8, 1.7, "→ " + take, size=15, bold=True, color=NAVY)
    # ---- WHY THESE B(c) CHOICES: rationale plots ----
    s = _slide(prs)
    _txt(s, 0.5, 0.3, 12.3, 0.5, "Why B depends on concentration — the Peclet regime", size=24, bold=True, color=NAVY)
    _txt(s, 0.5, 0.85, 12.3, 0.35, "Regression MSE vs Pe (≙ DATA2 Fig 8): which transport regime each salt sits in", size=14, bold=True, color=GREEN)
    _picpath(s, PEC / "peclet_mse.png", 0.3, 1.4, 12.7, 4.1)
    _txt(s, 0.5, 5.8, 12.3, 1.1,
         "→ Diluting NaCl → Pe→0 (diffusion-limited, B = DₘH/L = constant). CaCl₂ / concentrating-NaCl / LaCl₃ → Pe≫1 "
         "(convection+diffusion → B = Jw·Σβᵢcⁱ). The Peclet number is WHICH salts need a concentration-dependent B.",
         size=15, bold=True, color=NAVY)

    s = _slide(prs)
    _txt(s, 0.5, 0.3, 12.3, 0.5, "Which B(c) order — the DATA2 per-response AIC (not raw fit)", size=24, bold=True, color=NAVY)
    _txt(s, 0.5, 0.85, 12.3, 0.35, "Order chosen by the per-channel AIC (eq 39); raw WSSE3 always over-selects more parameters", size=14, bold=True, color=GREEN)
    _picpath(s, MS / "aic_selection.png", 0.3, 1.4, 12.7, 4.1)
    _txt(s, 0.5, 5.8, 12.3, 1.1,
         "→ CaCl₂ → quadratic / linear, cubic REJECTED (Δ=8–25); NaCl B is flat (constant). The per-response AIC is the "
         "honest criterion — the pooled single-n AIC stored in the JSONs over-picks cubic.",
         size=15, bold=True, color=NAVY)

    # ---- SUPPORTING PARAMETER CONTOURS (one per experiment, at the fitted form) ----
    cnote = {"MC3.07.22.24_SNaCl": "clear Lₚ minimum; σ rails (flat valley) but B is flat so the rail is harmless.",
             "MC2.05.07.24_NaCl": "constant-B identifiable in B–Lₚ; σ rails — adequate because B is concentration-independent.",
             "MC2.05.07.24_CaCl2": "with quadratic B(c) profiled, the permeate/retentate channels pin the landscape — the B(c) win.",
             "MC2.05.21.24_LaCl3": "shallow/flat objective — the model-limited case; the landscape itself shows poor identifiability."}
    for lbl, rid, tag, desc, take in EXP:
        cf, cdesc = contour_fig(rid, tag)
        s = _slide(prs)
        _txt(s, 0.5, 0.3, 12.3, 0.5, f"Supporting contour — {lbl.split('(')[0].strip()}", size=22, bold=True, color=NAVY)
        _txt(s, 0.5, 0.85, 12.3, 0.35, cdesc, size=14, bold=True, color=GREEN)
        _picpath(s, cf, 0.25, 1.35, 12.8, 4.0)
        _txt(s, 0.5, 5.6, 12.3, 1.5,
             "Objective landscape per channel (mass | permeate | retentate | combined), other params re-solved at each node; "
             "red ▲ = minimum. A sharp minimum = identified; a flat valley along σ = the σ-rail.\n→ "
             + cnote.get(rid, ""), size=14, bold=False, color=NAVY)

    # summary
    s = _slide(prs)
    _txt(s, 0.5, 0.4, 12.3, 0.6, "What the four fits show", size=26, bold=True, color=NAVY)
    _txt(s, 0.6, 1.4, 12.2, 4.5,
         "• NaCl (both regimes): constant B is sufficient — B is concentration-independent (Peclet → 0, diffusion-limited\n"
         "   for the diluting case). Mass and both concentrations are tracked at WSSE₃ ≈ 42 / 114.\n\n"
         "• CaCl₂: constant B fails (σ rails, WSSE₃=220). A quadratic B(c) is required and cuts the error ~33% to 147 —\n"
         "   consistent with convection+diffusion transport (Pe ≫ 1) and the DATA2 per-response AIC selecting quadratic.\n\n"
         "• LaCl₃: model-limited — every B(c) form sits at WSSE₃ ≈ 1190–1244; the residual points past B(c) alone\n"
         "   (e.g. σ(c) / ion-pairing for the trivalent).\n\n"
         "• Across the campaign: same NF270 membrane → Lₚ is consistent; the next step is the multi-experiment pooled fit\n"
         "   (shared σ + B(c)) that resolves the per-experiment σ-identifiability.",
         size=16)
    prs.save(str(OUT))
    print(f"saved {OUT}  ({len(prs.slides._sldIdLst)} slides)", flush=True)
    try:
        subprocess.run(["/opt/homebrew/bin/soffice", "--headless", "--convert-to", "pdf", "--outdir", str(HERE), str(OUT)],
                       timeout=180, check=False)
        print("QA PDF rendered" if (HERE / (OUT.stem + ".pdf")).exists() else "QA PDF not found", flush=True)
    except Exception as e:
        print("soffice skipped:", e, flush=True)
    print("FITTING DECK DONE", flush=True)


if __name__ == "__main__":
    build()
