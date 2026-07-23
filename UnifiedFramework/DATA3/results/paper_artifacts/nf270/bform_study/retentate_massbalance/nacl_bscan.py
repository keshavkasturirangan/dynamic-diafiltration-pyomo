import sys, json, numpy as np, io, contextlib, time
sys.path.insert(0,"/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1")
import refactored_ucb_library as lib
lib.NF270_CF_RESIDUAL_SCALE_FRACTION=0.02; lib.NF270_CF_RESIDUAL_FLOOR_MM=None
STUDY=lib.NF270_DEFAULT_ROOT.parent/"DATA3/results/paper_artifacts/nf270/bform_study"
NACL=["MC2.05.07.24_NaCl","MC5.07.23.24_NaCl","MC5.07.23.24_SNaCl","MC4.07.11.24_SNaCl","MC3.07.22.24_SNaCl","MC5.07.23.24_S2NaCl"]
def load(rid):
    fam=lib.NF270_RUN_REGISTRY[rid]; return lib.loadxlsx(lib.NF270_DEFAULT_ROOT/fam["workbook"],sheet=fam["sheet"])["data_stru"]
def scan(ds,mode,Lp,S0,S,sigma):
    Bg=np.linspace(1.0,30.0,15); rows=[]
    for B in Bg:
        th=dict(Lp=Lp,B=float(B),sigma=sigma,S0=S0,S=S)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                fit,s,_=lib.solve_model(ds,mode,theta=th,sim_opt=True,B_form="single",workflow_family="DATA3",nfe=120,solver_max_cpu_time=25)
            if isinstance(fit,dict): rows.append((B,fit["obj_m"],fit["obj_cv"],fit["obj_cr"]))
        except Exception: pass
    return np.array(rows)
print(f"{'sheet':22s} {'dir':5s} {'thru':>5s} {'C_D':>6s} {'B_mass':>7s} {'B_perm':>7s} {'B_reten':>8s} {'reten-perm':>10s}")
out={}
for rid in NACL:
    ds=load(rid); mode=ds["mode"]; cfg=ds["data_config"]
    anch=json.load(open(STUDY/rid/"result_single.json"))["parameters"]
    thru=float(np.nansum([ (lambda m: (m[np.isfinite(m)][-1]-m[np.isfinite(m)][0]) if np.isfinite(m).sum()>=2 else 0.0)(np.asarray(r["mass"],float)) for r in ds["data_raw"]]))/cfg["M_F0"]
    cf=[np.asarray(r.get("cF_exp",[]),float) for r in ds["data_raw"]]; cfall=np.concatenate(cf); cfall=cfall[np.isfinite(cfall)]
    direction="CONC" if cfall[-1]>cfall[0] else "DIL"
    arr=scan(ds,mode,anch["Lp"],anch["S0"],anch["S"],1.0)
    if arr.size==0: print(rid,"scan failed"); continue
    B=arr[:,0]
    Bm=B[np.argmin(arr[:,1])]; Bp=B[np.argmin(arr[:,2])]; Br=B[np.argmin(arr[:,3])]
    print(f"{rid:22s} {direction:5s} {thru:5.1f} {cfg['C_D']:6.1f} {Bm:7.1f} {Bp:7.1f} {Br:8.1f} {Br-Bp:10.1f}")
    out[rid]=dict(dir=direction,thru=thru,C_D=cfg["C_D"],Bm=float(Bm),Bp=float(Bp),Br=float(Br),arr=arr.tolist())
json.dump(out,open("/private/tmp/claude-505/-Users-kkasturi-GitHub-keshav-dev-dynamic-diafiltration-pyomo/d20a5790-00e6-4df6-bd33-bf6c3f39a8e5/scratchpad/bscan.json","w"))
print("DONE")
