#!/usr/bin/env python3
"""DATA3 figure-pitch deck: each slide = ONE plot + a publication-style caption (no titles).
Built from existing figures on disk. Caption style follows DATA1/DATA2 (J. Membr. Sci.)."""
import json
from pathlib import Path
from pptx import Presentation
from pptx.oxml import parse_xml
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from PIL import Image

HERE = Path(__file__).resolve().parent
B = HERE.parent / "UnifiedFramework/DATA3/results/paper_artifacts/nf270"
BF = B / "bform_study"
SW, SH = Inches(13.333), Inches(7.5)
NDBLUE = RGBColor(0x0C, 0x23, 0x40)
REPS = ["MC3.07.22.24_SNaCl", "MC2.05.07.24_NaCl", "MC2.05.07.24_CaCl2", "MC2.05.21.24_LaCl3"]

prs = Presentation(); prs.slide_width = SW; prs.slide_height = SH
BLANK = prs.slide_layouts[6]
_n = [0]


def _fg(rid, sub, kind):
    f = sorted((BF / "fit_improvement" / rid / sub).glob(f"{kind}-*.png"))
    return f[0] if f else None


def _reveal_on_click(slide, shape):
    """Inject the canonical PowerPoint <p:timing> so `shape` is HIDDEN when the
    slide first appears and FADES IN on the first mouse click — the audience reads
    the figure first, then clicks once to bring up the caption. This is the exact
    single-fade-entrance ('Appear/Fade', presetID=10, presetClass=entr) timing tree
    PowerPoint itself writes; python-pptx has no animation API so we append the XML
    to the <p:sld> element directly (schema order: cSld, clrMapOvr, transition,
    timing — appending last is valid since this deck sets no transition/extLst)."""
    spid = shape.shape_id
    xml = (
        '<p:timing xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">'
        '<p:tnLst><p:par><p:cTn id="1" dur="indefinite" restart="never" nodeType="tmRoot">'
        '<p:childTnLst><p:seq concurrent="1" nextAc="seek">'
        '<p:cTn id="2" dur="indefinite" nodeType="mainSeq"><p:childTnLst>'
        '<p:par><p:cTn id="3" fill="hold"><p:stCondLst><p:cond delay="indefinite"/></p:stCondLst>'
        '<p:childTnLst><p:par><p:cTn id="4" fill="hold"><p:stCondLst><p:cond delay="0"/></p:stCondLst>'
        '<p:childTnLst><p:par>'
        '<p:cTn id="5" presetID="10" presetClass="entr" presetSubtype="0" fill="hold" grpId="0" nodeType="clickEffect">'
        '<p:stCondLst><p:cond delay="0"/></p:stCondLst><p:childTnLst>'
        '<p:set><p:cBhvr><p:cTn id="6" dur="1" fill="hold"><p:stCondLst><p:cond delay="0"/></p:stCondLst></p:cTn>'
        f'<p:tgtEl><p:spTgt spid="{spid}"/></p:tgtEl>'
        '<p:attrNameLst><p:attrName>style.visibility</p:attrName></p:attrNameLst></p:cBhvr>'
        '<p:to><p:strVal val="visible"/></p:to></p:set>'
        '<p:animEffect transition="in" filter="fade"><p:cBhvr><p:cTn id="7" dur="500"/>'
        f'<p:tgtEl><p:spTgt spid="{spid}"/></p:tgtEl></p:cBhvr></p:animEffect>'
        '</p:childTnLst></p:cTn></p:par></p:childTnLst></p:cTn></p:par></p:childTnLst></p:cTn></p:par>'
        '</p:childTnLst></p:cTn>'
        '<p:prevCondLst><p:cond evt="onPrev" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond></p:prevCondLst>'
        '<p:nextCondLst><p:cond evt="onNext" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond></p:nextCondLst>'
        '</p:seq></p:childTnLst></p:cTn></p:par></p:tnLst></p:timing>'
    )
    slide._element.append(parse_xml(xml))


def _place(slide, im, box_l, box_t, box_w, box_h):
    """Fit image inside the box (inches), preserving aspect, centered."""
    w, h = Image.open(im).size
    ar = w / h
    if box_w / box_h > ar:                 # height-bound
        ph = box_h; pw = box_h * ar
    else:                                  # width-bound
        pw = box_w; ph = box_w / ar
    left = box_l + (box_w - pw) / 2.0
    top = box_t + (box_h - ph) / 2.0
    slide.shapes.add_picture(str(im), Inches(left), Inches(top), width=Inches(pw), height=Inches(ph))


def fig(images, caption):
    images = [Path(p) for p in images if p and Path(p).exists()]
    if not images:
        print(f"  SKIP (missing figs)"); return
    _n[0] += 1
    s = prs.slides.add_slide(BLANK)
    # plot region (top); caption region (bottom 1.85")
    if len(images) == 1:
        _place(s, images[0], 0.4, 0.25, 12.53, 5.15)
    elif len(images) == 2:
        for i, im in enumerate(images):
            _place(s, im, 0.3 + i * 6.43, 0.4, 6.2, 4.9)
    else:  # 4 → 2×2
        for i, im in enumerate(images):
            r, c = divmod(i, 2)
            _place(s, im, 0.4 + c * 6.43, 0.25 + r * 2.6, 6.1, 2.55)
    tb = s.shapes.add_textbox(Inches(0.5), Inches(5.55), SW - Inches(1.0), Inches(1.85))
    tf = tb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]
    rb = p.add_run(); rb.text = f"Figure {_n[0]}. "; rb.font.bold = True; rb.font.size = Pt(13); rb.font.color.rgb = NDBLUE
    rc = p.add_run(); rc.text = caption; rc.font.size = Pt(13); rc.font.color.rgb = RGBColor(0x11, 0x11, 0x11)
    _reveal_on_click(s, tb)   # figure shows first; caption fades in on the first click
    print(f"  Figure {_n[0]} added (caption = click-to-reveal)")


# ---------------------------------------------------------------- title slide
_t = prs.slides.add_slide(BLANK)
_tb = _t.shapes.add_textbox(Inches(0.9), Inches(2.5), SW - Inches(1.8), Inches(2.6))
_tf = _tb.text_frame; _tf.word_wrap = True
_p = _tf.paragraphs[0]; _r = _p.add_run()
_r.text = "Identifiability of solute-transport parameters in NF270 single-salt diafiltration"
_r.font.size = Pt(32); _r.font.bold = True; _r.font.color.rgb = NDBLUE
_p2 = _tf.add_paragraph(); _r2 = _p2.add_run()
_r2.text = "Figure pitch"
_r2.font.size = Pt(18); _r2.font.color.rgb = RGBColor(0x55, 0x55, 0x55)

# ---------------------------------------------------------------- figures + captions
# Representative experiments: (label, run_id, extra note). The numeric Lp/B/σ in
# the captions are read LIVE from each sheet's regenerated single-form fit
# (bform_study/<rid>/result_single.json) by _fit_params(), so the caption numbers
# always match the plotted curves after a refit. _FIT_FALLBACK holds the last
# hard-coded values and is only used if a result_single.json is missing.
REP4 = [
    ("diluting NaCl", "MC3.07.22.24_SNaCl", ""),
    ("concentrating NaCl", "MC2.05.07.24_NaCl", ""),
    ("concentrating CaCl₂", "MC2.05.07.24_CaCl2", ""),
    ("concentrating LaCl₃", "MC2.05.21.24_LaCl3",
     " This trivalent-salt experiment is the most demanding case in the campaign and "
     "shows the largest residuals of the four."),
]

_FIT_FALLBACK = {  # (Lp, B, σ-bound word, σ value) — pre vial-close refit; fallback only
    "MC3.07.22.24_SNaCl": (8.98, 13.88, "upper", "1"),
    "MC2.05.07.24_NaCl":  (9.91, 10.00, "upper", "1"),
    "MC2.05.07.24_CaCl2": (7.54, 5.71,  "upper", "1"),
    "MC2.05.21.24_LaCl3": (3.54, 0.32,  "lower", "0"),
}


def _fit_params(rid):
    """(Lp, B, σ-bound word, σ value-string) read live from result_single.json.

    σ rails to a bound in the single-experiment fit, so it is reported as the
    nearest bound (upper→"1", lower→"0"). Falls back to the last hard-coded
    values (with a warning) if the regenerated fit is not present yet."""
    rp = BF / rid / "result_single.json"
    if rp.exists():
        try:
            p = json.loads(rp.read_text()).get("parameters", {})
            lp, b, sig = float(p["Lp"]), float(p["B"]), float(p.get("sigma", 1.0))
            bound = "upper" if sig >= 0.5 else "lower"
            return lp, b, bound, ("1" if sig >= 0.5 else "0")
        except Exception as e:
            print(f"  WARN {rid}: could not read result_single.json ({e}); using fallback numbers")
    else:
        print(f"  WARN {rid}: result_single.json missing; using fallback caption numbers")
    return _FIT_FALLBACK[rid]


def _pred(rid, kind):
    f = sorted((BF / "predictions" / rid / "single").glob(f"{kind}-*.png"))
    return f[0] if f else None


# (Figures 1-4) representative model-vs-measurement fits
for label, rid, note in REP4:
    lp, b, bound, sval = _fit_params(rid)
    fig([_pred(rid, "mass"), _pred(rid, "concentration")],
        f"Measured and model-predicted time profiles for a representative {label} diafiltration "
        f"experiment. (Left) total retentate mass and (right) solute concentration in the retentate and "
        f"permeate streams versus time; symbols are measurements — the retentate inductively-coupled-plasma "
        f"(ICP) reference points and the time-resolved permeate vial samples — and the curves are the "
        f"lumped-parameter model evaluated at the fitted parameters reported in the inset "
        f"(Lₚ = {lp:.2f} L m⁻² h⁻¹ bar⁻¹, B = {b:.2f} µm s⁻¹, σ = {sval}). A single concentration-independent "
        f"solute permeability reproduces the mass decline and both concentration histories, with the retentate "
        f"and permeate plotted in distinct colors; the reflection coefficient is forced to its {bound} bound "
        f"(σ = {sval}).{note}")

# (Figures 5-8) sigma x B identifiability at the fixed, identified Lp
for label, rid, note in REP4:
    lp, b, bound, sval = _fit_params(rid)
    fig([BF / "Bsigma_fixedLp_nofloor" / rid / "single" / "contour_titled.png"],
        f"Identifiability of the reflection coefficient σ and the solute permeability B for a representative "
        f"{label} experiment, with the hydraulic permeability fixed at its independently identified optimum "
        f"Lₚ = {lp:.2f} L m⁻² h⁻¹ bar⁻¹ (well determined and insensitive to σ and B). The three panels are the base-10 logarithm of the weighted "
        f"sum-of-squared-error (WSSE) surface for each measured response — labelled mass, permeate, and "
        f"retentate — computed by forward simulation over the σ×B grid at that fixed Lₚ; the red triangle "
        f"marks each surface's minimum. Each surface constrains B through a narrow valley but is only weakly "
        f"sensitive to σ, so the three responses reach their minima at different σ values along the boundaries; "
        f"the reflection coefficient is therefore not jointly identified, and the single-experiment optimum "
        f"places σ at its {bound} bound (σ = {sval}). The solute permeability is well determined while σ is not "
        f"recoverable from one experiment.")

fig([BF / "peclet" / "peclet_mse.png", BF / "peclet" / "peclet_jsfit.png"],
    "(Left) Mean squared error of the convection–diffusion partition-coefficient regression for the fitted "
    "coefficients (h₀, h₁) as a function of the Péclet number Pe = J_w L / D_m, and (right) the simulated solute "
    "flux J_s versus interfacial concentration cᵢₙ,f with the convection–diffusion fit overlaid; both panels show "
    "the four salt/flow regimes. Following the framework of the DATA2 study, J_w, cᵢₙ,f, c_h, and J_s are "
    "forward-simulated from each experiment's fitted model and the interfacial partition is regressed against Pe. "
    "The three convection-dominated regimes—concentrating NaCl (Pe*≈77), CaCl₂ (Pe*≈36), and LaCl₃ (Pe*≈36)—"
    "minimize their regression error at high Pe, whereas diluting NaCl is the diffusion-limited counterexample, "
    "with its error minimized at Pe*→0 (≈0.001) in its own panel. Because neither panel plots B directly, this "
    "high-Pe behavior is taken as a mechanistic inference for, rather than a direct demonstration of, the "
    "concentration-dependent solute permeability, and only for the convection-dominated salts.")

fig([BF / "peclet" / "peclet_regimes.png"],
    "Operating Péclet number Pe*=J_w·L/D_m for the four representative experiments, placed on a logarithmic "
    "Péclet axis with the diffusion (Pe<1) and convection (Pe>1) regions shaded and the Pe=1 crossover marked. "
    "Each Pe* is the value that minimizes the convection–diffusion regression error of the preceding figure. "
    "Diluting NaCl is diffusion-limited (Pe*≪1, here ≈10⁻³) while the three concentrating salts—NaCl (Pe*≈77), "
    "CaCl₂ (Pe*≈36), and LaCl₃ (Pe*≈36)—are convection-dominated (Pe*≫1), demonstrating that the same membrane "
    "spans distinct transport regimes—the mechanistic basis for the concentration-dependent solute permeability.")

fig([BF / "taylor_vs_donnan" / "donnan_reconciliation.png"],
    "Normalized Donnan–dielectric solute permeability B/B_plateau versus the similarity variable "
    "u = c/c_knee = 2kc/X — the single master curve onto which the fitted partition "
    "B(c) = J_w·P₀·(−X+√(X²+4k²c²))/(2c) collapses for every salt, with plateau B_plateau = J_w·P₀·k and "
    "knee c_knee = X/(2k). The curve rises approximately linearly through the charge-exclusion regime (u < 1, "
    "where the curvature determines the fixed-charge scale X) and saturates onto the screened plateau for u ≫ 1. "
    "Shaded bands are the actual measured concentration windows of the three salts whose Donnan fit converged, "
    "mapped onto u at their own fitted (X, k). Because the fits force X to its lower bound (X = 10⁻³ mM), the knee "
    "sits at c_knee ≈ 10⁻³ mM, so every window falls at u ≈ 10³–10⁵ and samples B within ~0.07 % of the plateau. "
    "The saturating Donnan shape is therefore the correct mechanistic form, but a single salt probes only its flat "
    "top: the curvature that would identify X is never measured. This is why the fitted Donnan reference is "
    "effectively constant within the measured window, why even a constant or low-order Taylor reproduces it "
    "there, and why X (and the reflection coefficient σ) cannot be recovered from one salt — resolving the knee "
    "would require data at far lower concentration or several salts fitted jointly.")

fig([BF / "taylor_vs_donnan" / "MC2.05.21.24_LaCl3_taylor_vs_donnan.png"],
    "Solute permeability B as a function of interfacial concentration for a representative LaCl₃ diafiltration "
    "experiment, comparing candidate concentration-dependent models, each evaluated at its own fitted "
    "parameters. The mechanistic Donnan–dielectric partition (solid black) is the reference; first-, second-, and "
    "third-order Taylor polynomials (linear, quadratic, cubic) and a saturating exponential are fitted to the same "
    "data. The Donnan reference is nearly flat (≈0.3 µm s⁻¹) because its fitted knee lies far below the measured "
    "window: within the measured range the partition is on its screened plateau, so a constant is the leading "
    "behavior. (Left) Within the window all forms therefore agree closely, scattering by less than ±0.2 µm s⁻¹ "
    "about this plateau. (Right) Extrapolated beyond the window the polynomials diverge unphysically—the cubic "
    "turns negative and the quadratic curls upward—whereas the saturating exponential collapses onto the Donnan "
    "plateau. A low-order Taylor expansion thus captures the partition only within the calibration range, and the "
    "in-window flatness is the direct consequence of the partition sitting on its screened plateau.")

fig([BF / "taylor_vs_donnan" / "MC3.07.11.24_SCaCl2_taylor_vs_donnan.png"],
    "Solute permeability B as a function of interfacial concentration for a representative CaCl₂ diafiltration "
    "experiment, comparing candidate concentration-dependent models, each evaluated at its own fitted "
    "parameters. The mechanistic Donnan–dielectric partition (solid black) is the reference; a single "
    "concentration-independent constant, first-, second-, and third-order Taylor polynomials (linear, "
    "quadratic, cubic), and a saturating exponential are fitted to the same data. The Donnan "
    "reference is on its screened plateau across the whole window (fitted knee far below the data), so it is "
    "nearly constant at B ≈ 8–9 µm s⁻¹. (Left) Within the window all forms therefore agree closely about this "
    "plateau. (Right) Extrapolated beyond the window the polynomials diverge unphysically while the Donnan "
    "reference stays bounded, the correct mechanistic behavior, confirming that a low-order Taylor expansion "
    "captures the partition only within its calibration range.")

fig([BF / "taylor_vs_donnan" / "MC5.07.23.24_NaCl_taylor_vs_donnan.png"],
    "Solute permeability B as a function of interfacial concentration for a representative NaCl diafiltration "
    "experiment, comparing candidate concentration-dependent models, each evaluated at its own fitted "
    "parameters. The mechanistic Donnan–dielectric partition (solid black) is the reference; only the "
    "concentration-independent constant, second-order (quadratic), and third-order (cubic) Taylor "
    "polynomials are shown, because the linear and saturating-exponential fits did not converge for this sheet. "
    "As in Figure 11 the Donnan reference is on its screened plateau across the window, so it is nearly constant "
    "at B ≈ 10 µm s⁻¹. (Left) Within the window the shown forms therefore agree closely about this plateau. "
    "(Right) Extrapolated beyond the window the cubic polynomial diverges unphysically (turning sharply negative) "
    "while the Donnan reference stays bounded, the correct mechanistic behavior, confirming that a low-order "
    "Taylor expansion captures the partition only within its calibration range.")

fig([BF / "conc_flux" / "MC2.05.07.24_CaCl2" / "conc_flux.png"],
    "Forward-simulated membrane-interface concentration cᵢₙ,f and bulk (stirred-cell) concentration c_h (left), "
    "water flux J_w (center), and solute flux J_s (right) versus time for a representative concentrating CaCl₂ "
    "diafiltration experiment, computed from the fitted lumped-parameter model. As the stirred cell concentrates, "
    "both concentrations rise and J_w declines through the increasing osmotic pressure, while J_s increases with "
    "concentration—consistent with a concentration-dependent solute permeability. The gap between the membrane-"
    "interface concentration cᵢₙ,f and the lower bulk concentration c_h is the concentration polarization: solute "
    "accumulates at the membrane face relative to the well-mixed bulk.")

out = HERE / "DATA3_figure_pitch.pptx"
prs.save(str(out))
print(f"\nSAVED {len(prs.slides._sldIdLst)} slides -> {out}")
