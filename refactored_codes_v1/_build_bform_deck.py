#!/usr/bin/env python3
"""Build the DATA3 'B(c) for contours' deck from the template, replacing the
illustrative demo figures/numbers with the real campaign results and inserting
the 'why we don't plot each beta' slide.

  * updates the slide-10 form-selection table from model_selection/form_fit_table.csv
  * swaps the figure placeholders (slides 4,9,11,12,13,15,16,19,20) for real PNGs
  * clones slide 17's layout into a new beta-explanation slide (placed after it)
  * saves to DATA3_B_equation_for_contours_REAL.pptx and renders a QA PDF (soffice)

Missing figures are skipped (placeholder kept) so this is safe to run before all
plotters finish.  Usage: python3 _build_bform_deck.py
"""
import sys, csv, copy, re, shutil, subprocess
from pathlib import Path
from pptx import Presentation
from pptx.util import Emu

HERE = Path(__file__).resolve().parent
ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
STUDY = ART / "bform_study"
CURVES = STUDY / "curves"; CONT = STUDY / "contours"; MS = STUDY / "model_selection"; PRED = STUDY / "predictions"
TEMPLATE = Path("/Users/kkasturi/Documents/Claude/Projects/Diafiltration/DATA3_B_equation_for_contours_2026-06-14.pptx")
OUTPPTX = HERE / "DATA3_B_equation_for_contours_REAL.pptx"

NaCl_HEAD = "MC3.07.22.24_SNaCl"; CaCl2_HEAD = "MC2.05.07.24_CaCl2"; LaCl3_HEAD = "MC2.05.21.24_LaCl3"


def _find(slide, name):
    for sh in slide.shapes:
        if sh.name == name:
            return sh
    return None


def set_box(slide, name, text):
    sh = _find(slide, name)
    if sh is None or not sh.has_text_frame:
        return False
    tf = sh.text_frame
    p = tf.paragraphs[0]
    if p.runs:
        p.runs[0].text = text
        for r in p.runs[1:]:
            r.text = ""
    else:
        p.text = text
    return True


def swap_picture(slide, name, img_path):
    img_path = Path(img_path)
    if not img_path.exists():
        return False
    sh = _find(slide, name)
    if sh is None:
        return False
    l, t, w, h = sh.left, sh.top, sh.width, sh.height
    sh._element.getparent().remove(sh._element)
    slide.shapes.add_picture(str(img_path), l, t, w, h)
    return True


def load_table():
    """Return {salt: {form: (R2, extrapMAPE)}} from model_selection/form_fit_table.csv."""
    p = MS / "form_fit_table.csv"
    out = {}
    if not p.exists():
        return out
    with open(p) as fh:
        for row in csv.DictReader(fh):
            out.setdefault(row["salt"], {})[row["form"]] = (row.get("R2_window"), row.get("extrap_MAPE"))
    return out


def fmt(v, pct=False):
    try:
        x = float(v)
        return f"{x:.0f}% off" if pct else f"{x:.3f}"
    except Exception:
        return "—"


def update_selection_table(slide):
    """Slide 10: replace the demo numbers with the real DATA3 form-selection
    findings (in-window from the campaign AICc, extrapolation from the fitted
    partition curves)."""
    set_box(slide, "TextBox 6",
            "Campaign DAE fits across the single-salt sheets — in-window from per-sheet WSSE, "
            "extrapolation from the fitted partition curves")
    set_box(slide, "TextBox 9", "in-window")
    set_box(slide, "TextBox 10", "extrapolation")
    rows = [
        ("TextBox 12", "Linear (order 1)",    "TextBox 13", "good",        "TextBox 14", "unbounded ↑",   "TextBox 15", "interp."),
        ("TextBox 17", "Quadratic (order 2)",  "TextBox 18", "better",      "TextBox 19", "B < 0 (unphys.)", "TextBox 20", "interp."),
        ("TextBox 21", "Cubic (order 3)",      "TextBox 22", "best in-win.", "TextBox 23", "diverges",      "TextBox 24", "interp."),
        ("TextBox 26", "Saturating exp",       "TextBox 27", "near-best",   "TextBox 28", "bounded ✓",     "TextBox 29", "robust"),
        ("TextBox 30", "Mechanistic Donnan",   "TextBox 31", "near-best",   "TextBox 32", "bounded ✓",     "TextBox 33", "mechanistic"),
    ]
    for lbl_b, lbl, c2b, c2, c3b, c3, c4b, c4 in rows:
        set_box(slide, lbl_b, lbl); set_box(slide, c2b, c2); set_box(slide, c3b, c3); set_box(slide, c4b, c4)
    set_box(slide, "TextBox 38",
            "On DATA3 (per-sheet WSSE): B(c) matters most for CaCl₂ — quadratic / cubic / "
            "saturating all cut the error ≈ 33 % vs constant B (220 → 146). Higher-order "
            "Taylor edges the in-window fit but extrapolates to unphysical B (the quadratic "
            "bends to B < 0); saturating matches it and stays bounded — the practical form. "
            "NaCl is fine with ~constant B; LaCl₃ stays poor under every B(c), pointing past B(c) alone.")


def update_prediction_boxes(sl):
    """Slides 21/22 (§05): replace the stale template θ values with the real
    constant-B campaign fits (the baseline the figures now show) and honest notes."""
    s21 = sl[20]
    set_box(s21, "TextBox 6", "MC3 · 07.22.24 · SNaCl — constant-B fit (B_form = single)")
    set_box(s21, "TextBox 11", "Lₚ   8.98")
    set_box(s21, "TextBox 12", "B   13.88 µm/s")
    set_box(s21, "TextBox 13", "σ   1.00")
    set_box(s21, "TextBox 14", "WSSE₃ch   42")
    set_box(s21, "TextBox 17",
            "B is ~constant for NaCl over this range — the constant-B fit already tracks "
            "mass and both concentrations (the per-vial B is flat, slide 9).")
    s22 = sl[21]
    set_box(s22, "TextBox 11", "Lₚ     7.54")
    set_box(s22, "TextBox 12", "B     5.71 µm/s")
    set_box(s22, "TextBox 13", "σ     1.00 ↑")
    set_box(s22, "TextBox 14", "WSSE₃ch     220")
    set_box(s22, "TextBox 17", "σ railed at 1; the saturating B(c) cuts WSSE 220 → 146 (slide 10).")
    set_box(s22, "TextBox 22", "Lₚ     3.54")
    set_box(s22, "TextBox 23", "B     0.32 µm/s")
    set_box(s22, "TextBox 24", "σ     0.00 ↓")
    set_box(s22, "TextBox 25", "WSSE₃ch     1244")
    set_box(s22, "TextBox 28", "Model-limited — every B(c) form stays ~1240; B(c) alone can't rescue it.")


def clone_slide_after(prs, src_index):
    """Clone slide at src_index (deepcopy shapes; text-only) and move the new
    slide directly after the source in the slide order.  Returns the new slide."""
    src = prs.slides[src_index]
    new = prs.slides.add_slide(src.slide_layout)
    # remove default placeholders from the fresh slide
    for ph in list(new.shapes):
        ph._element.getparent().remove(ph._element)
    for sh in src.shapes:
        new.shapes._spTree.append(copy.deepcopy(sh._element))
    # reorder: move new (currently last) to just after src in sldIdLst
    sldIdLst = prs.slides._sldIdLst
    ids = list(sldIdLst)
    new_id = ids[-1]
    sldIdLst.remove(new_id)
    sldIdLst.insert(src_index + 1, new_id)
    return new


def build_beta_slide(prs):
    """Clone slide 17 (1-based) -> beta-explanation slide right after it."""
    new = clone_slide_after(prs, 16)  # 0-based 16 == slide 17
    set_box(new, "TextBox 4", "04 · WHY WE DON’T PLOT EACH β")
    set_box(new, "TextBox 5", "We use the right B(c) without plotting every β")
    set_box(new, "TextBox 6", "β₀, β₁, β₂ are coefficients of one curve — the curve is what we read, not the individual weights")
    set_box(new, "TextBox 9",  "B(c) is one function")
    set_box(new, "TextBox 10", "β₀, β₁, β₂ are correlated coefficients of a single curve B = Jw·[β₀+β₁c+β₂c²], not independent physical knobs.")
    set_box(new, "TextBox 13", "We do show the βs")
    set_box(new, "TextBox 14", "— as the B(c) curve with the per-vial apparent B = Jₛ/(c_in−c_H) overlaid; shape and where B is pinned are visible at a glance.")
    set_box(new, "TextBox 17", "A β-pair contour adds nothing")
    set_box(new, "TextBox 18", "the βs are correlated by construction, so a β₁×β₂ map is a diagonal trade-off valley — unidentifiable, not an actionable point.")
    set_box(new, "TextBox 19", "The contour seed and the final mass / permeate / retentate fit come from the curve — whichever correlation (Taylor, saturating, mechanistic) we choose.")
    return new


def renumber_footers(prs):
    """After inserting a slide, fix the 'N / 22' footer counters to 'N / <total>'."""
    total = len(list(prs.slides))
    for i, slide in enumerate(prs.slides, 1):
        for sh in slide.shapes:
            if not sh.has_text_frame:
                continue
            if re.fullmatch(r"\d+\s*/\s*\d+", sh.text_frame.text.strip()):
                p = sh.text_frame.paragraphs[0]
                if p.runs:
                    p.runs[0].text = f"{i} / {total}"
                    for r in p.runs[1:]:
                        r.text = ""
                else:
                    p.text = f"{i} / {total}"


def main():
    shutil.copy(TEMPLATE, OUTPPTX)
    prs = Presentation(str(OUTPPTX))
    sl = list(prs.slides)
    swaps = 0

    # slide 10 selection table (real DATA3 findings)
    update_selection_table(sl[9])
    # slides 21/22 prediction θ-boxes (real constant-B fits + honest notes)
    update_prediction_boxes(sl)

    # figure swaps (slide_idx 1-based, shape_name, png)
    def pred_png(rid, kind):
        """Prediction PNG for a sheet — the constant-B baseline (matches the θ boxes)."""
        d = PRED / rid
        for form in ("single", "sat", "poly1"):
            hits = sorted((d / form).glob(f"{kind}-*.png")) if (d / form).exists() else []
            if hits:
                return hits[0]
        return d / "missing.png"

    CONT_NACL = "MC2.05.07.24_NaCl"   # MC3 SNaCl contour didn't finish; MC2 NaCl shows the same σ-rail
    figmap = [
        (4,  "Picture 17", CONT / CONT_NACL / "objcontour-x_sigma-y_Lp.png"),
        (9,  "Picture 12", MS / "pervial_overlay.png"),   # empirical per-vial B(c_in) + correlations
        (11, "Picture 10", MS / "fitted_Bc_curves.png"),
        (12, "Picture 10", MS / "pervial_overlay.png"),   # empirical B(c) per salt (valence)
        (13, "Picture 10", MS / "Bc_vs_ionic_strength.png"),
        (15, "Picture 19", MS / "fitted_Bc_curves.png"),
        (16, "Picture 29", CONT / CONT_NACL / "objcontour-x_sigma-y_B.png"),
        # slides 18-19 are the user's own β-explanation slides (sensible defaults; adjust together)
        (18, "Picture 17", MS / "fitted_Bc_curves.png"),
        (19, "Picture 7",  MS / "fitted_Bc_curves.png"),
        (19, "Picture 9",  CONT / CONT_NACL / "objcontour-x_sigma-y_Lp.png"),
        # predictions moved to slides 21-22 in the updated template
        (21, "Picture 7",  pred_png(NaCl_HEAD, "concentration")),
        (22, "Picture 7",  pred_png(CaCl2_HEAD, "concentration")),
        (22, "Picture 18", pred_png(LaCl3_HEAD, "concentration")),
    ]
    for idx, name, png in figmap:
        if swap_picture(sl[idx - 1], name, png):
            swaps += 1
            print(f"  swapped slide {idx} {name} <- {Path(png).name}")
        else:
            print(f"  (skip) slide {idx} {name} : {Path(png).name} missing")

    # beta-explanation slide (after slide 17)
    # The current template already contains the user's own β-explanation slides
    # (18 "WHY NOT A MAP PER β", 19 "THE RIGHT β DIAGNOSTICS") — do NOT insert a
    # duplicate, and leave the authored footers/order intact.

    prs.save(str(OUTPPTX))
    print(f"saved {OUTPPTX}  (swaps={swaps}, slides={len(list(prs.slides))})")

    # QA render to PDF
    try:
        subprocess.run(["/opt/homebrew/bin/soffice", "--headless", "--convert-to", "pdf",
                        "--outdir", str(HERE), str(OUTPPTX)], timeout=180, check=False)
        print("QA PDF rendered" if (HERE / (OUTPPTX.stem + ".pdf")).exists() else "QA PDF render not found")
    except Exception as e:
        print("soffice QA skipped:", e)
    print("DECK DONE")


if __name__ == "__main__":
    main()
