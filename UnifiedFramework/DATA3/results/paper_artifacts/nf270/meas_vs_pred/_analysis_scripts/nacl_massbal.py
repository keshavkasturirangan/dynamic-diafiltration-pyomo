import sys, numpy as np
sys.path.insert(0,"/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1")
import refactored_ucb_library as lib
NACL=["MC2.05.07.24_NaCl","MC5.07.23.24_NaCl","MC5.07.23.24_SNaCl","MC4.07.11.24_SNaCl","MC3.07.22.24_SNaCl","MC5.07.23.24_S2NaCl"]
def load(rid):
    fam=lib.NF270_RUN_REGISTRY[rid]
    return lib.loadxlsx(lib.NF270_DEFAULT_ROOT/fam["workbook"], sheet=fam["sheet"])["data_stru"]

def per_vial(ds):
    """Return arrays over vials: dm (permeate mass gain g), cV (permeate ICP mM),
    cF_close (retentate conc at vial close mM)."""
    dm=[]; cV=[]; cF=[]
    for row in ds["data_raw"]:
        m=np.asarray(row["mass"],float).reshape(-1); mf=m[np.isfinite(m)]
        dm.append(float(mf[-1]-mf[0]) if mf.size>=2 else 0.0)
        cv=np.asarray(row.get("cV_avg",[]),float).reshape(-1); cvf=cv[np.isfinite(cv)]
        cV.append(float(cvf[-1]) if cvf.size else np.nan)
        cf=np.asarray(row.get("cF_exp",[]),float).reshape(-1); cff=cf[np.isfinite(cf)]
        cF.append(float(cff[-1]) if cff.size else np.nan)
    return np.array(dm),np.array(cV),np.array(cF)

print(f"{'sheet':22s} {'dir':5s} {'M_F0':>5s} {'M_O(weigh)':>10s} {'mF_end(solute)':>15s} {'ratio':>6s} {'closureR%':>9s}")
results={}
for rid in NACL:
    ds=load(rid); cfg=ds["data_config"]
    M_F0=cfg["M_F0"]; C_F0=cfg["C_F0"]; C_D=cfg["C_D"]; M_O=cfg["M_O"]
    dm,cV,cF=per_vial(ds)
    mV_total=float(np.nansum(dm))
    # permeate solute out (mM*g); use vials with finite cV
    perm_solute=float(np.nansum(np.where(np.isfinite(cV),cV,0.0)*dm))
    cF_final=cF[np.isfinite(cF)][-1]
    direction="CONC" if cF_final>C_F0 else "DIL"
    # solve effective final cell mass from solute+water balance (NOT assuming const mass):
    #   mF_end = [M_F0*(C_F0-C_D) + C_D*mV_total - perm_solute] / (cF_final - C_D)
    num=M_F0*(C_F0-C_D)+C_D*mV_total-perm_solute
    den=(cF_final-C_D)
    mF_end_solute=num/den
    # constant-mass closure residual (feed=mV_total): solute_in - solute_out_accounted, with mF=M_F0
    R=(C_F0*M_F0 + C_D*mV_total - perm_solute) - cF_final*M_F0
    Rpct=R/(C_F0*M_F0+C_D*mV_total)*100
    print(f"{rid:22s} {direction:5s} {M_F0:5.1f} {M_F0+M_O:10.2f} {mF_end_solute:15.2f} {mF_end_solute/M_F0:6.2f} {Rpct:9.1f}")
    results[rid]=dict(M_F0=M_F0,M_O=M_O,C_F0=C_F0,C_D=C_D,dm=dm,cV=cV,cF=cF,mV_total=mV_total,mF_end_solute=mF_end_solute,dir=direction)
import json
np.save("/private/tmp/claude-505/-Users-kkasturi-GitHub-keshav-dev-dynamic-diafiltration-pyomo/d20a5790-00e6-4df6-bd33-bf6c3f39a8e5/scratchpad/nacl_bal.npy", results, allow_pickle=True)
print("\nInterpretation: mF_end(solute)=effective cell mass implied by the retentate+permeate+collected-mass balance.")
print("  ~10 g (=M_F0) => mass conserved, assumption holds.  >>10 g => cell holds more water than weighed (feed accumulated).")
