import sys, json, io, contextlib
sys.path.insert(0,"/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1")
import refactored_ucb_library as lib
import _run_bform_campaign as camp
lib.NF270_CF_RESIDUAL_FLOOR_MM=1.0
STUDY=lib.NF270_DEFAULT_ROOT.parent/"DATA3/results/paper_artifacts/nf270/bform_study"
for rid in ["MC5.07.23.24_NaCl","MC5.07.23.24_SNaCl"]:
    ds=camp.load(rid); mode=ds["mode"]
    base=json.load(open(STUDY/rid/"result_single.json"))["parameters"]
    B=base.get("B",5.0); Lp=base["Lp"]; sg=base["sigma"]; S0=base["S0"]; S=base["S"]
    best=None
    for Binf in (B*0.7,B,B*1.5,B*3,B*6,50.0):
        for cstar in (1.0,5.0,20.0,60.0,150.0,300.0):
            seed=dict(Lp=Lp,sigma=sg,S0=S0,S=S,B_inf=float(Binf),c_star=float(cstar))
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    fit,_,_=lib.solve_model(ds,mode,theta=seed,sim_opt=False,B_form="sat",workflow_family="DATA3",nfe=150,multistart=False,solver_max_cpu_time=60)
                if isinstance(fit,dict):
                    w=sum(float(fit.get(k,0.0) or 0.0) for k in ("obj_m","obj_cv","obj_cr"))
                    if best is None or w<best[0]: best=(w,fit)
            except Exception: pass
    if best is None:
        # last resort: multistart
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                fit,_,_=lib.solve_model(ds,mode,theta=dict(Lp=Lp,sigma=sg,S0=S0,S=S,B_inf=B*1.3,c_star=20.0),sim_opt=False,B_form="sat",workflow_family="DATA3",nfe=150,multistart=True,multistart_iterations=12,solver_max_cpu_time=300)
            if isinstance(fit,dict):
                w=sum(float(fit.get(k,0.0) or 0.0) for k in ("obj_m","obj_cv","obj_cr")); best=(w,fit)
        except Exception: pass
    if best:
        w,fit=best; p=fit["parameters"]
        res={"parameters":p,"WSSE3":w,"obj_m":fit.get("obj_m"),"obj_cv":fit.get("obj_cv"),"obj_cr":fit.get("obj_cr"),"run_id":rid,"salt":"NaCl","form":"sat","refit":True}
        json.dump(res,open(STUDY/rid/"result_sat.json","w"),indent=2,default=float)
        print(f"[sat OK] {rid}: WSSE3={w:.1f} B_inf={p['B_inf']:.3f} c_star={p['c_star']:.2f} sigma={p['sigma']:.3f}",flush=True)
    else:
        print(f"[sat FAIL] {rid}",flush=True)
print("SAT DONE")
