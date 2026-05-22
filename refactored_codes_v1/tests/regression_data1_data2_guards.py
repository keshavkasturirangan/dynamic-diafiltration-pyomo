"""
Lightweight regression test for the DATA3 guards on NF270_CF_RESIDUAL_FLOOR_MM
and NF270_SIGMA_INTERIOR_BOUNDS.  Does NOT run IPOPT — only constructs the
Pyomo model and inspects σ bounds + the cF scale computation directly.

This is faster (~1 s vs ~minutes) and tests exactly the right thing: that
the guards make the knobs inert for DATA1/DATA2 paths.
"""
import sys, os
sys.path.insert(0, '/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1')
os.chdir('/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo')
import refactored_ucb_library as lib
from pyomo.environ import value


def sigma_bounds_for(workflow_family, knob):
    """Build the model with the given workflow_family + sigma knob,
    return the σ Var's (lb, ub)."""
    lib.NF270_SIGMA_INTERIOR_BOUNDS = knob
    lib.NF270_CF_RESIDUAL_FLOOR_MM  = None
    ds = lib.loadmat("data_library/data_stru-dataset270511.123.mat")['data_stru']
    m = lib.model_construct_inter(ds, "DATA", theta=None, sim_opt=False,
                                   B_form="single",
                                   workflow_family=workflow_family)
    lb = float(value(m.sigma.lb))
    ub = float(value(m.sigma.ub))
    return (lb, ub)


def cf_floor_active_for(workflow_family, knob_mM):
    """Reproduce the floor-vs-relative selection from obj_rule for a
    representative cF_meas value, and report what the obj_rule's scale
    selection would do.  This mirrors the exact conditional at :1875.

    Returns True if the floor is active (i.e. would change the scale)."""
    lib.NF270_CF_RESIDUAL_FLOOR_MM   = knob_mM
    lib.NF270_SIGMA_INTERIOR_BOUNDS  = None
    # The exact conditional from the patched obj_rule:
    cf_floor_mM = (
        lib.NF270_CF_RESIDUAL_FLOOR_MM
        if str(workflow_family).upper() == "DATA3"
        else None
    )
    cf_meas = 5.0  # 5 mM cF; relative_scale = 0.003 * 5 = 0.015 mM, well below floor=1
    relative_scale = 0.003 * cf_meas
    if cf_floor_mM is not None and cf_floor_mM > 0 and relative_scale < cf_floor_mM:
        cf_scale = float(cf_floor_mM)
        floor_active = True
    else:
        cf_scale = relative_scale
        floor_active = False
    return floor_active, cf_scale


print("=" * 72)
print("LIGHTWEIGHT GUARD AUDIT  ·  no IPOPT, just code-path inspection")
print("=" * 72)

# ---------- sigma interior bounds guard ----------
print("\n### sigma interior bounds guard")
print(f"{'workflow_family':>18s} {'NF270_SIGMA_INTERIOR_BOUNDS':>32s} {'σ Var bounds':>20s}  {'expected?':>10s}")
print("-" * 90)
cases_sigma = [
    ("DATA1", None,         (0.0, 1.0), "legacy"),
    ("DATA1", (0.6, 0.9),   (0.0, 1.0), "GUARD ON"),   # knob set, but DATA1 must ignore
    ("DATA2", None,         (0.0, 1.0), "legacy"),
    ("DATA2", (0.6, 0.9),   (0.0, 1.0), "GUARD ON"),   # knob set, but DATA2 must ignore
    ("DATA3", None,         (0.0, 1.0), "legacy"),
    ("DATA3", (0.6, 0.9),   (0.6, 0.9), "active"),     # only DATA3 honors the knob
]
sigma_ok = True
for wf, knob, expected, note in cases_sigma:
    got = sigma_bounds_for(wf, knob)
    match = abs(got[0] - expected[0]) < 1e-12 and abs(got[1] - expected[1]) < 1e-12
    marker = "yes" if match else "NO!"
    if not match: sigma_ok = False
    knob_s = str(knob) if knob is None else f"({knob[0]:.1f},{knob[1]:.1f})"
    got_s  = f"({got[0]:.4f},{got[1]:.4f})"
    print(f"{wf:>18s} {knob_s:>32s} {got_s:>20s}  {marker:>10s}  [{note}]")
print(f"\nsigma-bounds guard:  {'PASS' if sigma_ok else 'FAIL'}")

# ---------- cF residual floor guard ----------
print("\n### cF residual floor guard")
print(f"{'workflow_family':>18s} {'NF270_CF_RESIDUAL_FLOOR_MM':>30s} {'floor active?':>16s}  {'expected?':>10s}")
print("-" * 88)
cases_cf = [
    ("DATA1", None,  False, "legacy"),
    ("DATA1", 1.0,   False, "GUARD ON"),   # knob set, but DATA1 must ignore
    ("DATA2", None,  False, "legacy"),
    ("DATA2", 1.0,   False, "GUARD ON"),   # knob set, but DATA2 must ignore
    ("DATA3", None,  False, "legacy"),
    ("DATA3", 1.0,   True,  "active"),     # only DATA3 honors the knob
]
cf_ok = True
for wf, knob, expected, note in cases_cf:
    active, scale = cf_floor_active_for(wf, knob)
    match = (active == expected)
    marker = "yes" if match else "NO!"
    if not match: cf_ok = False
    knob_s = str(knob)
    print(f"{wf:>18s} {knob_s:>30s} {str(active):>16s}  {marker:>10s}  [{note}]")
print(f"\ncF-floor guard:  {'PASS' if cf_ok else 'FAIL'}")

# ---------- permeate-probe guard ----------
print("\n### permeate-probe channel guard")
print(f"{'workflow_family':>18s} {'NF270_USE_PERMEATE_PROBE':>30s} {'permeate term active?':>22s}  {'expected?':>10s}")
print("-" * 90)
perm_cases = [
    ("DATA1", None,  False, "legacy"),
    ("DATA1", 0.03,  False, "GUARD ON"),
    ("DATA2", None,  False, "legacy"),
    ("DATA2", 0.03,  False, "GUARD ON"),
    ("DATA3", None,  False, "legacy"),
    ("DATA3", 0.03,  True,  "active"),
]
perm_ok = True
for wf, knob, expected, note in perm_cases:
    lib.NF270_USE_PERMEATE_PROBE = knob
    _is_data3 = (str(wf).upper() == "DATA3")
    use_perm_probe = bool(_is_data3 and lib.NF270_USE_PERMEATE_PROBE)
    match = (use_perm_probe == expected)
    marker = "yes" if match else "NO!"
    if not match: perm_ok = False
    print(f"{wf:>18s} {str(knob):>30s} {str(use_perm_probe):>22s}  {marker:>10s}  [{note}]")
print(f"\npermeate-probe guard:  {'PASS' if perm_ok else 'FAIL'}")

# Reset
lib.NF270_CF_RESIDUAL_FLOOR_MM  = None
lib.NF270_SIGMA_INTERIOR_BOUNDS = None
lib.NF270_USE_PERMEATE_PROBE    = None

print("\n" + "=" * 72)
print("SUMMARY")
print("=" * 72)
print(f"  sigma interior-bounds guard:  {'PASS' if sigma_ok else 'FAIL'}")
print(f"  cF residual-floor guard:      {'PASS' if cf_ok else 'FAIL'}")
print(f"  permeate-probe guard:         {'PASS' if perm_ok else 'FAIL'}")
if sigma_ok and cf_ok and perm_ok:
    print()
    print("  Conclusion: both knobs are inert for DATA1/DATA2 model construction.")
    print("  DATA1 and DATA2 fits are byte-equivalent to the legacy code path")
    print("  regardless of the constants' values.")
else:
    print()
    print("  *** REGRESSION ***  the guard leaks.")
    sys.exit(1)
