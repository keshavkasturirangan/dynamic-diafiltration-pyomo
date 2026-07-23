import sys, json, io, contextlib, math
sys.path.insert(0,"/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1")
import refactored_ucb_library as lib
import _run_bform_campaign as camp
lib.NF270_CF_RESIDUAL_FLOOR_MM=1.0  # campaign cfg
STUDY=lib.NF270_DEFAULT_ROOT.parent/"DATA3/results/paper_artifacts/nf270/bform_study"
MISSING=[("MC4.07.11.24_SNaCl",2,"poly2"),("MC5.07.23.24_NaCl",1,"poly1"),
         ("MC5.07.23.24_NaCl","sat","sat"),("MC5.07.23.24_SNaCl","sat","sat")]
NB={"poly1":2,"poly2":3,"sat":2}
for rid,Bform,tag in MISSING:
    ds=camp.load(rid); mode=ds["mode"]
    base=json.load(open(STUDY/rid/"result_single.json"))["parameters"]
    seed={k:base[k] for k in ("Lp","sigma","S0","S") if k in base}
    B=base.get("B",1.0)
    if tag=="poly1": seed.update(beta_0=B, beta_1=0.0)
    elif tag=="poly2":
        p1=json.load(open(STUDY/rid/"result_poly1.json")).get("parameters")
        if p1: seed.update(beta_0=p1["beta_0"],beta_1=p1["beta_1"],beta_2=0.0)
        else: seed.update(beta_0=B,beta_1=0.0,beta_2=0.0)
    elif tag=="sat": seed.update(B_inf=max(1e-3,B*1.3), c_star=20.0)
    ok=False
    for ms in (False,True):
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                fit,_,_=lib.solve_model(ds,mode,theta=dict(seed),sim_opt=False,B_form=Bform,workflow_family="DATA3",
                                        nfe=150,multistart=ms,multistart_iterations=8,solver_max_cpu_time=200,LOUD=False)
            if isinstance(fit,dict):
                wsse3=sum(float(fit.get(k,0.0) or 0.0) for k in ("obj_m","obj_cv","obj_cr"))
                res={"parameters":fit["parameters"],"WSSE3":wsse3,"obj_m":fit.get("obj_m"),"obj_cv":fit.get("obj_cv"),
                     "obj_cr":fit.get("obj_cr"),"run_id":rid,"salt":"NaCl","form":tag,"refit":True}
                json.dump(res,open(STUDY/rid/f"result_{tag}.json","w"),indent=2,default=float)
                print(f"[refit OK] {rid} {tag}: WSSE3={wsse3:.1f} params={ {k:round(v,4) for k,v in fit['parameters'].items() if k not in ('S0','S')} }",flush=True)
                ok=True; break
        except Exception as e:
            print(f"  {rid} {tag} ms={ms} FAIL {str(e)[:50]}",flush=True)
    if not ok: print(f"[refit FAIL] {rid} {tag}",flush=True)
print("REFIT DONE")
