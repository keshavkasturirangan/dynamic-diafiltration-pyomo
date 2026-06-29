#!/usr/bin/env python3
"""Inject speaker notes into DATA3_figure_pitch.pptx (one note per slide, shown in
PowerPoint Presenter View) and also write DATA3_figure_pitch_NOTES.md.
Notes explain each plot in plain language (a coffee-filter running example) + the math."""
from pathlib import Path
from pptx import Presentation

HERE = Path(__file__).resolve().parent
DECK = HERE / "DATA3_figure_pitch.pptx"

FRAME = ("Running picture — squeezing salty water through a coffee filter:\n"
         "  • Lp = how EASILY water passes (bigger Lp → more water per push)\n"
         "  • σ  = how well the filter BLOCKS salt by osmosis (σ=1 blocks all, σ=0 lets salt ride along)\n"
         "  • B  = how much salt LEAKS through by diffusion\n"
         "  • ΔP = how hard you push;  Δπ = the osmotic 'back-push' from the salt difference")

CORE_MATH = ("Core model (3 fitted parameters Lp, σ, B):\n"
             "  Water flux   Jw = Lp·(ΔP − σ·Δπ)          [net push = applied − osmotic back-push]\n"
             "  Osmotic      Δπ = i·R·T·(cF − cP)         [i = ions/molecule: NaCl 2, CaCl₂ 3, LaCl₃ 4]\n"
             "  Salt flux    Js = B·(cF − cP) + (1−σ)·Jw·c̄  [diffusion term + convection term]\n"
             "  Mass balance dM/dt = −ρ·A·Jw              [retentate lightens as permeate leaves]")


def fit_note(salt, lp, b, sig, bound, extra=""):
    return (
        f"FIT CHECK — {salt}. Does the model reproduce the experiment?\n\n"
        f"What it shows: lines = model, symbols = measurements. LEFT = retentate mass vs time; "
        f"RIGHT = concentration vs time (TEAL = retentate, RED = permeate; teal diamonds = retentate ICP at start & end).\n\n"
        f"Simple steps: (1) choose Lp, σ, B; (2) march the equations forward in time; (3) overlay on the data. "
        f"The inset box reports the best-fit values.\n\n"
        f"Coffee-filter example: as water is pushed out, the cup (retentate) gets saltier and lighter while the "
        f"drips (permeate) carry a little salt.\n\n"
        f"{CORE_MATH}\n\n"
        f"This run: Lp = {lp}, B = {b} µm/s, σ = {sig} (σ forced to its {bound} bound). {extra}\n\n{FRAME}")


def contour_note(salt, lp, sig, bound):
    return (
        f"IDENTIFIABILITY MAP — {salt}. Can the data tell σ and B apart?\n\n"
        f"What it shows: with Lp FIXED at its identified value (Lp = {lp}), we map the fit error over σ (x-axis) and "
        f"B (y-axis). Three panels = the three things we measure (mass, permeate, retentate). Red triangle = best fit.\n\n"
        f"Simple steps: (1) fix Lp; (2) try every (σ, B) on a grid; (3) simulate, score the error, color the map "
        f"(blue = good fit, red = bad). Titles say mass/permeate/retentate; the color is log₁₀ of the weighted error (WSSE).\n\n"
        f"How to read it: the contours run nearly FLAT along the σ direction → moving σ barely changes the error → "
        f"σ is NOT identifiable and is forced to its {bound} bound (σ = {sig}). The error DOES change with B → B is "
        f"well determined.\n\n"
        f"Coffee-filter example: the filter blocks salt so completely that you can't measure exactly how well — σ just "
        f"pegs at 'fully blocking.'\n\n"
        f"Math: WSSE = Σ[(model − data)/noise]², plotted as log₁₀. A long flat valley = many equally-good answers = "
        f"not identifiable; a tight bowl = identifiable.")


NOTES = {}

# 0 — title
NOTES[0] = (
    "OPENING (30 s). We fit three transport parameters to single-salt NF270 diafiltration data and ask which are "
    "really identifiable.\n\n"
    "Punchline of the whole talk: Lp (water permeability) and B (salt permeability) are well determined; σ (salt "
    "rejection) is NOT — it is forced to its bound; and the concentration-dependence of B cannot be located from a "
    "single salt.\n\n" + FRAME)

# 1-4 — fits
NOTES[1] = fit_note("diluting NaCl", "8.98", "13.88", "1", "upper")
NOTES[2] = fit_note("concentrating NaCl", "9.91", "10.00", "1", "upper")
NOTES[3] = fit_note("concentrating CaCl₂", "7.54", "5.71", "1", "upper")
NOTES[4] = fit_note("concentrating LaCl₃", "3.54", "0.32", "0", "lower",
                    extra="This trivalent salt is the hardest case — largest residuals of the four.")

# 5 — Lp identifiability 3x3
NOTES[5] = (
    "WHY WE CAN FIX Lp. This 3×3 map proves the water permeability Lp is pinned down.\n\n"
    "What it shows: 9 error maps. COLUMNS = the 3 measured responses (mass, permeate, retentate). ROWS = parameter "
    "pairs (σ–Lp, B–Lp, σ–B). Color = log of fit error; red triangle = best.\n\n"
    "Simple steps: (1) sweep two parameters on a grid, hold the third at its best value; (2) simulate & score the "
    "error at each grid point; (3) plot the surface.\n\n"
    "How to read it: a tight HORIZONTAL band in the top two rows means Lp lands at the same value no matter what σ "
    "and B do → Lp is identifiable. The bottom row (σ–B) is a long diagonal VALLEY → σ and B trade off and can't be "
    "separated.\n\n"
    "Coffee-filter example: if many different filter settings give the same drip pattern, you can't tell them apart "
    "— that flat valley is exactly that.\n\n"
    "Why it matters: because Lp is well pinned, we are ALLOWED to fix it and just study σ-vs-B in the next four slides.\n\n"
    "Math: error = WSSE = Σ[(model − data)/noise]²; we plot log₁₀(WSSE).")

# 6-9 — sigma x B contours
NOTES[6] = contour_note("diluting NaCl", "8.98", "1", "upper")
NOTES[7] = contour_note("concentrating NaCl", "9.91", "1", "upper")
NOTES[8] = contour_note("concentrating CaCl₂", "7.54", "1", "upper")
NOTES[9] = contour_note("concentrating LaCl₃", "3.54", "0", "lower")

# 10 — Peclet MSE/Js
NOTES[10] = (
    "WHY B DEPENDS ON CONCENTRATION (diagnostic 1). \n\n"
    "What it shows: LEFT — regression error vs Péclet number for each salt; RIGHT — simulated salt flux Js vs "
    "concentration with a convection–diffusion fit overlaid.\n\n"
    "Math: Péclet number Pe = Jw·L/Dm = (transport by water being pushed through) ÷ (transport by diffusion). "
    "L = membrane thickness, Dm = diffusivity. Pe > 1 means convection wins; Pe < 1 means diffusion wins.\n\n"
    "Simple steps: forward-simulate Jw and the concentrations, then regress how salt partitions vs Pe.\n\n"
    "How to read it: the three concentrating salts minimize their error at HIGH Pe (convection-dominated → B varies "
    "with c); diluting NaCl is the counter-example at LOW Pe (diffusion → B roughly constant).\n\n"
    "Coffee-filter example: a slow drip (diffusion) and a firehose (convection) through the same filter behave "
    "differently — that's why the salt permeability isn't a single fixed number.")

# 11 — Peclet regimes
NOTES[11] = (
    "WHY B DEPENDS ON CONCENTRATION (diagnostic 2 — the one-line version).\n\n"
    "What it shows: one number per salt — its operating Péclet — on a log axis. Blue zone = diffusion (Pe < 1), "
    "orange zone = convection (Pe > 1), line at Pe = 1.\n\n"
    "Punchline: the SAME membrane spans BOTH regimes. Diluting NaCl sits deep in diffusion (Pe ≈ 0.001); the three "
    "concentrating salts sit in convection (Pe ≈ 36–77). That is the mechanistic reason the salt permeability is "
    "concentration-dependent for the concentrating runs.\n\n"
    "Math reminder: Pe = Jw·L/Dm. Crossing Pe = 1 is the switch from diffusion- to convection-controlled transport.")

# 12 — Donnan reconciliation
NOTES[12] = (
    "THE SHAPE OF B(c) — and why our fit looks flat.\n\n"
    "Math: the mechanistic Donnan partition is\n"
    "   B(c) = Jw·P0·(−X + √(X² + 4k²c²)) / (2c).\n"
    "Two knobs: the PLATEAU height P0·k (reached at high salt) and the KNEE location c_knee = X/(2k) (where the "
    "curve bends). At low c it rises ~linearly from zero; at high c it flattens to the plateau.\n\n"
    "Simple steps: rescale every salt onto ONE universal master curve using u = c/c_knee, then mark where each "
    "salt's MEASURED window lands on it.\n\n"
    "How to read it: every measured window sits far out on the FLAT plateau (u ≈ 10³–10⁵). That's because the fit "
    "pushes the charge term X to its lower bound, so the knee falls ~1000× below where we actually measure.\n\n"
    "Punchline: the saturating shape is the correct physics, but a single salt only samples the flat TOP of it — so "
    "we cannot locate the knee (the curvature) from one salt at a time. Honest framing: 'one salt can't see it,' "
    "NOT 'the membrane is flat.'\n\n"
    "Coffee-filter example: we're standing on the flat summit far past the slope — we know we're high up, but can't "
    "see where the hill rose.")

# 13-15 — Taylor vs Donnan
def taylor_note(salt, level):
    return (
        f"WHICH B(c) FORM TO FIT — {salt}.\n\n"
        f"What it shows: the mechanistic Donnan curve (solid black) vs simple polynomial 'Taylor' fits and a "
        f"saturating-exp fit, INSIDE the measured window (left) and EXTRAPOLATED beyond it (right).\n\n"
        f"Math: a Taylor series B(c) ≈ β₀ + β₁c + β₂c² is just a local polynomial approximation — accurate near the "
        f"data, unreliable far from it.\n\n"
        f"How to read it: inside the window all forms agree (the Donnan reference is on its flat plateau ≈ {level}, "
        f"so a constant already captures it). OUTSIDE the window the polynomials blow up unphysically (cubic turns "
        f"negative, quadratic curls up) while Donnan stays bounded.\n\n"
        f"Punchline: use a low-order Taylor INSIDE the calibration range — it's identifiable and sufficient — but "
        f"do NOT extrapolate it; for transfer across conditions use the mechanistic form.\n\n"
        f"Coffee-filter example: a straight-line guess between two known points works between them but goes crazy "
        f"outside them.")
NOTES[13] = taylor_note("LaCl₃", "0.3 µm/s")
NOTES[14] = taylor_note("CaCl₂", "8–9 µm/s")
NOTES[15] = taylor_note("NaCl", "10 µm/s")

# 16 — conc & flux
NOTES[16] = (
    "THE FULL PICTURE OVER TIME (CaCl₂).\n\n"
    "What it shows: three panels vs time. A — interface concentration c_in,f vs bulk c_h (their GAP is concentration "
    "polarization, salt piling up at the membrane face). B — water flux Jw: model (blue) + MEASURED (black dotted). "
    "C — salt flux Js.\n\n"
    "Simple steps: forward-simulate at the fitted parameters; get the measured water flux from how fast permeate "
    "mass is collected in each vial.\n\n"
    "Math: measured Jw = (dm/dt)/(ρ·A) — collection rate ÷ (density × membrane area); model Jw = Lp·(ΔP − σ·Δπ).\n\n"
    "How to read it: as the cell concentrates, both concentrations rise → osmotic back-push grows → water flux Jw "
    "declines (and the model tracks the measured dotted line) → salt flux Js climbs. All consistent with "
    "concentration-dependent transport.\n\n"
    "Coffee-filter example: the saltier the cup gets, the harder osmosis pushes back, so the drip slows down over time.")


def main():
    prs = Presentation(str(DECK))
    md = ["# DATA3 figure-pitch — speaker notes\n",
          "_Running example throughout: squeezing salty water through a coffee filter "
          "(Lp = water passes easily, σ = blocks salt by osmosis, B = salt leaks by diffusion)._\n"]
    # Slide order after the Lp-identifiability 3×3 panel (old note index 5) was
    # removed from the deck: skip 5, so contours..conc&flux shift down one slot.
    order = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    assert len(order) == len(prs.slides._sldIdLst), "note map length != slide count"
    for i, slide in enumerate(prs.slides):
        note = NOTES.get(order[i], "")
        slide.notes_slide.notes_text_frame.text = note
        title = "Title slide" if i == 0 else f"Slide {i+1} — Figure {i}"
        md.append(f"\n## {title}\n\n{note}\n")
    prs.save(str(DECK))
    (HERE / "DATA3_figure_pitch_NOTES.md").write_text("\n".join(md))
    print(f"Injected notes into {prs.slides and len(prs.slides._sldIdLst)} slides")
    print(f"Wrote {HERE/'DATA3_figure_pitch_NOTES.md'}")


if __name__ == "__main__":
    main()
