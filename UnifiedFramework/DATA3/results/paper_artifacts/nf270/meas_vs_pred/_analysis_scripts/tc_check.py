import sys, json, numpy as np, io, contextlib
sys.path.insert(0,"/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1")
import refactored_ucb_library as lib
lib.NF270_CF_RESIDUAL_SCALE_FRACTION=0.02; lib.NF270_CF_RESIDUAL_FLOOR_MM=None
STUDY=lib.NF270_DEFAULT_ROOT.parent/"DATA3/results/paper_artifacts/nf270/bform_study"
rid="MC2.05.07.24_NaCl"; fam=lib.NF270_RUN_REGISTRY[rid]
p=json.load(open(STUDY/rid/"result_single.json"))["parameters"]
def load(conv):
    ds=lib.loadxlsx(lib.NF270_DEFAULT_ROOT/fam["workbook"],sheet=fam["sheet"])["data_stru"]
    lib.apply_campaign_time_correction(ds, convention=conv)
    return ds
def sim_traj(ds):
    with contextlib.redirect_stdout(io.StringIO()):
        fit,s,_=lib.solve_model(ds,ds["mode"],theta=dict(p),sim_opt=True,B_form="single",workflow_family="DATA3",nfe=120,solver_max_cpu_time=40)
    t=[];mF=[];mV=[];cIn=[]
    for i in sorted(s.keys()):
        t+=list(np.asarray(s[i]["time"],float)); mF+=list(np.asarray(s[i]["mF"],float)); mV+=list(np.asarray(s[i]["mV"],float)); cIn+=list(np.asarray(s[i]["cIn"],float))
    return np.array(t),np.array(mF),np.array(mV),np.array(cIn)
print("=== cV_avg (permeate ICP) placement index per vial under each convention ===")
for conv in ("vial_close","tube_transit","none"):
    ds=load(conv)
    idxs=[]
    for row in ds["data_raw"]:
        cv=np.asarray(row["cV_avg"],float).reshape(-1); fin=np.where(np.isfinite(cv))[0]
        idxs.append(f"{fin[-1]}/{cv.size}" if fin.size else "-")
    print(f"  {conv:12s}: {idxs}")
print("\n=== do the MODEL trajectories (mF,mV,cIn) differ at FIXED theta across conventions? ===")
base=None
for conv in ("vial_close","tube_transit","none"):
    ds=load(conv); t,mF,mV,cIn=sim_traj(ds)
    if base is None: base=(mF,mV,cIn); print(f"  {conv:12s}: reference"); continue
    dmF=np.max(np.abs(mF-base[0])); dmV=np.max(np.abs(mV-base[1])); dcIn=np.max(np.abs(cIn-base[2]))
    print(f"  {conv:12s}: max|Δ mF|={dmF:.2e} g  max|Δ mV|={dmV:.2e} g  max|Δ cIn|={dcIn:.2e} mM")
print("\n=== does the FIT (theta) differ across conventions? (re-fit single B) ===")
for conv in ("vial_close","tube_transit","none"):
    ds=load(conv)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fit,_,_=lib.solve_model(ds,ds["mode"],sim_opt=False,B_form="single",workflow_family="DATA3",nfe=150,multistart=False,solver_max_cpu_time=120)
        q=fit["parameters"]; print(f"  {conv:12s}: Lp={q['Lp']:.2f} B={q['B']:.2f} sigma={q['sigma']:.3f}  obj_cv={fit['obj_cv']:.1f} obj_cr={fit['obj_cr']:.1f}")
    except Exception as e: print(f"  {conv:12s}: fit failed {str(e)[:40]}")
