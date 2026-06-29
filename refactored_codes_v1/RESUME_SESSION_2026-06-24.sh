#!/usr/bin/env bash
# =============================================================================
# RESUME_SESSION_2026-06-24.sh
# Session memory + reproduce script for the DATA3 work done on 2026-06-23/24.
# Safe to read top-to-bottom; safe to run (default prints a summary, nothing heavy
# runs unless you ask for a subcommand).
#
#   bash RESUME_SESSION_2026-06-24.sh            # print the session summary
#   bash RESUME_SESSION_2026-06-24.sh verify     # quick verifications (fast)
#   bash RESUME_SESSION_2026-06-24.sh deck        # rebuild + validate the collab deck
#   bash RESUME_SESSION_2026-06-24.sh figures     # regenerate the deck's source figures
#   bash RESUME_SESSION_2026-06-24.sh holdup      # regenerate curated-holdup table + plots
#   bash RESUME_SESSION_2026-06-24.sh env         # print the two python envs used
#
# Run from refactored_codes_v1/.   Full prose record:
#   - refactored_codes_v1/Architecture.md  §20   (what we did, with file refs)
#   - refactored_codes_v1/DATA3_SESSION_HANDOFF_2026-06-23.md  (what's left)
# =============================================================================
cd "$(dirname "$0")" || exit 1

# ---- environments -----------------------------------------------------------
# Two different pythons were used this session:
PY_SCI="$HOME/miniforge3/envs/my-idaes-env/bin/python"   # numpy, pyomo, IPOPT, openpyxl, matplotlib  (loader, fits, ports)
PY_DECK="$HOME/miniforge3/bin/python"                     # python-pptx, PIL, matplotlib              (deck building)
PANDOC_BIN="$HOME/miniforge3/bin"                         # pandoc (LaTeX -> OMML native equations)
SKILL_PPTX=""   # set to the pptx skill dir if you want soffice/validate (see render_deck)

env_info () {
  echo "PY_SCI (fits/ports/loader): $PY_SCI"
  "$PY_SCI" -c "import numpy,pyomo,openpyxl,matplotlib;print('  -> numpy/pyomo/openpyxl/matplotlib OK')" 2>/dev/null || echo "  -> MISSING deps"
  command -v "$PY_SCI" >/dev/null && "$PY_SCI" - <<'P' 2>/dev/null
import shutil;print("  -> ipopt:",shutil.which("ipopt") or "(add env bin to PATH)")
P
  echo "PY_DECK (deck build): $PY_DECK"
  "$PY_DECK" -c "import pptx,PIL,matplotlib;print('  -> python-pptx/PIL/matplotlib OK')" 2>/dev/null || echo "  -> MISSING deps"
  echo "pandoc: $PANDOC_BIN/pandoc  (needed for native equation-mode OMML)"
}

# =============================================================================
summary () {
cat <<'TXT'
================ DATA3 SESSION SUMMARY (2026-06-23 / 2026-06-24) ================
ALL work is DATA3 / NF270 only. DATA1, DATA2, and the legacy MATLAB are UNTOUCHED.

WORKSTREAM A — Vial-close permeate time correction (DATA2 convention)
  * Retired the invented V_tube tube-transit model; the single ICP permeate value
    is anchored at vial close (compared vs model cV at tau.last()). One shared
    time axis for mass/retentate/permeate.
  * Files: refactored_ucb_library.py (_load_legacy_data_stru_from_excel ~:4783;
    metadata; _nf270_corrected_cv_index/NF270_TUBE_VOLUME_G deprecated),
    _preflight_audit.py, _audit_loader_fidelity.py, Architecture.md.
  * Verified: DATA2/smoke pytest 98 passed (separate .mat path unaffected).

WORKSTREAM B — Curated holdup-boundary table (DATA2-style override)
  * nf270_holdup_boundaries.csv (per-sheet holdup_end_s; blank => auto-detect).
  * Loader honors it: _nf270_curated_holdup_end_s() + the vial-1 split block.
  * Generators: _make_holdup_boundary_table.py (pre-fill + diagnostics),
    _curate_holdup_boundaries.py (collection-line extrapolation to mass=0).
  * Diagnostics: holdup_diagnostics/  and  holdup_diagnostics/curated/.

WORKSTREAM C — MATLAB analysis functions ported to Python (no MATLAB changed)
  * New in refactored_ucb_library.py: calc_contour_3d_py, doe_heatmap_py,
    heatmap_sigma_sensitivity_py  (existing: calc_ind_objectives_py,
    calc_contour_2d_py, sigma_sensitivity_py, calc_FIM).
  * Driver: run_nf270_matlab_ports(); runfile branch 'x' in refactored_ucb_runfile.py.
  * All evaluate through solve_model/calc_FIM => DATA3 model + error spec.
  * Per-sheet mode propagation (was hardcoded "Lag"); 3-agent audit; bugs fixed.

WORKSTREAM D — DATA3-spec toggle (DATA1/DATA2 untouched)
  * NF270_FORCE_DATA3_SPEC (default False) + run_with_data3_spec() honored by
    _nf270_conc_scales; wired into the runfile NF270 dispatch only.

WORKSTREAM E — Collaborator-update deck (built 2026-06-24)
  * DATA3_collaborator_update_2026-06-24.pptx          (8 slides, main)
  * DATA3_collaborator_update_2026-06-24_backup.pptx   (13 slides, backup)
  * Builder: _build_collab_update_deck.py   (assets in _collab_update_build/)
  * Clean ND-navy on WHITE bg, dark text; native EDITABLE equations (OMML via
    pandoc + _eq2omml.py); figures kept as crisp images (edit at source).
  * Flow: Title -> Feedback -> Topics -> sigma-B bridge -> Peclet -> Donnan
    master curve -> Taylor-vs-Donnan LaCl3 -> Taylor-vs-Donnan CaCl2.
  * 5-core (drop the 3 explanatory: Feedback, Topics, sigma-B bridge).
  * Both decks pass the pptx OOXML validator.

PENDING (not yet run — heavy, run on a machine with IPOPT)
  * The vial-close change shifts fitted B/sigma => regenerate DATA3 fits, figures,
    and the figure-pitch deck. See DATA3_SESSION_HANDOFF_2026-06-23.md "Stage 1-3".
  * Open the collab deck in PowerPoint once to confirm equations are editable
    (cannot be verified outside PowerPoint).
================================================================================
TXT
}

# ---- reproduce / continue ---------------------------------------------------
verify () {
  echo ">>> DATA2 + smoke pytest (expect green; confirms DATA1/2 paths untouched)"
  ( cd .. && PATH="$HOME/miniforge3/envs/my-idaes-env/bin:$PATH" python -m pytest \
      pytest_refactored_codes_v1/tests/test_data1_smoke.py \
      pytest_refactored_codes_v1/tests/test_data2_paper_comparison.py -q 2>&1 | tail -3 )
  echo ">>> vial-close sanity (permeate at last index) — MC3.07.22.24_SNaCl"
  PATH="$HOME/miniforge3/envs/my-idaes-env/bin:$PATH" "$PY_SCI" - <<'P'
import sys; sys.path.insert(0,".")
import numpy as np, refactored_ucb_library as lib
fam=lib.NF270_RUN_REGISTRY["MC3.07.22.24_SNaCl"]
ds=lib.loadxlsx(lib.NF270_DEFAULT_ROOT/fam["workbook"],sheet=fam["sheet"])["data_stru"]
print("  permeate_time_correction:",ds.get("permeate_time_correction",{}).get("method"))
print("  holdup_boundary_source :",ds["data_config"].get("holdup_boundary_source"))
P
}

deck () {
  echo ">>> rebuild collaborator deck (native equations via pandoc)"
  PATH="$PANDOC_BIN:$PATH" "$PY_DECK" _build_collab_update_deck.py
  echo ">>> validate (run the pptx skill validator if available):"
  echo "    python <pptx-skill>/scripts/office/validate.py DATA3_collaborator_update_2026-06-24.pptx"
}

figures () {
  echo ">>> regenerate the deck's source figures (plot-only; reads cached result_*.json)"
  export PATH="$HOME/miniforge3/envs/my-idaes-env/bin:$PATH"
  python _make_peclet.py
  python _make_donnan_reconciliation.py
  python _make_taylor_vs_donnan_lacl3.py
  python _make_taylor_vs_donnan.py MC3.07.11.24_SCaCl2
  echo "    (then re-run:  bash $0 deck )"
}

holdup () {
  echo ">>> regenerate curated-holdup table + diagnostics"
  export PATH="$HOME/miniforge3/envs/my-idaes-env/bin:$PATH"
  python _make_holdup_boundary_table.py
  python _curate_holdup_boundaries.py
}

case "${1:-summary}" in
  summary) summary ;;
  env)     env_info ;;
  verify)  verify ;;
  deck)    deck ;;
  figures) figures ;;
  holdup)  holdup ;;
  *) echo "unknown: $1"; echo "use: summary | env | verify | deck | figures | holdup"; exit 1 ;;
esac
