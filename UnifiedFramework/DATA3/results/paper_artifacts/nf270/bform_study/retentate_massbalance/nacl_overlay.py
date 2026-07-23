import sys, json, numpy as np, io, contextlib
sys.path.insert(0,"/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1")
import refactored_ucb_library as lib
lib.NF270_CF_RESIDUAL_SCALE_FRACTION=0.02; lib.NF270_CF_RESIDUAL_FLOOR_MM=None
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
rid="MC2.05.07.24_NaCl"; fam=lib.NF270_RUN_REGISTRY[rid]
ds=lib.loadxlsx(lib.NF270_DEFAULT_ROOT/fam["workbook"],sheet=fam["sheet"])["data_stru"]; mode=ds["mode"]
anch=json.load(open(lib.NF270_DEFAULT_ROOT.parent/"DATA3/results/paper_artifacts/nf270/bform_study"/rid/"result_single.json"))["parameters"]
Lp=anch["Lp"]; S0=anch["S0"]; S=anch["S"]
def sim(B,sigma):
    th=dict(Lp=Lp,B=B,sigma=sigma,S0=S0,S=S)
    with contextlib.redirect_stdout(io.StringIO()):
        fit,s,_=lib.solve_model(ds,mode,theta=th,sim_opt=True,B_form="single",workflow_family="DATA3",nfe=120,solver_max_cpu_time=40)
    return s
def traj(s,key):
    t=[]; y=[]
    for i in sorted(s.keys()):
        t+=list(np.asarray(s[i]["time"],float)); y+=list(np.asarray(s[i][key],float))
    return np.array(t),np.array(y)
# measured data
mt=[]; mcf=[]; pt=[]; pcv=[]; gt=[]; gm=[]
n_v0=ds["data_config"].get("n_v0",1)
mprev=0.0
for i,row in enumerate(ds["data_raw"]):
    tt=np.asarray(row["time"],float); 
    cf=np.asarray(row.get("cF_exp",[]),float)
    m=np.asarray(row["mass"],float)
    mt+=list(tt); mcf+=list(cf)
    # cumulative permeate mass
    mf=m[np.isfinite(m)]
    if mf.size: gt+=list(tt); gm+=list(mprev+ (m-mf[0] if np.isfinite(m[0]) else m)); mprev+=float(mf[-1]-mf[0])
    cv=np.asarray(row.get("cV_avg",[]),float); fin=np.where(np.isfinite(cv))[0]
    if fin.size: pt.append(tt[fin[-1]]); pcv.append(cv[fin[-1]])
mt=np.array(mt); mcf=np.array(mcf)
t0=mt[np.isfinite(mt)][0]
cases=[("mass-opt  B=1.7 σ=0.4",1.70,0.40,"tab:green"),
       ("perm-opt  B=6.4 σ=1.0",6.36,1.0,"tab:blue"),
       ("reten-opt B=17.2 σ=1.0",17.23,1.0,"tab:red"),
       ("joint fit B=10.0 σ=1.0",9.997,1.0,"0.4")]
fig,ax=plt.subplots(1,4,figsize=(22,5))
for lab,B,sg,col in cases:
    s=sim(B,sg)
    for j,key in enumerate(["cF","cV","mV","mF"]):
        t,y=traj(s,key); ax[j].plot((t-t0)/60,y,color=col,lw=1.6,label=lab)
# overlay measured
ax[0].plot((mt-t0)/60,mcf,"k.",ms=2,alpha=0.5,label="measured cF (probe)")
ax[1].plot((np.array(pt)-t0)/60,pcv,"ks",ms=6,label="permeate ICP")
ax[2].plot((np.array(gt)-t0)/60,gm,"k.",ms=2,alpha=0.4,label="measured mass")
ax[3].axhline(10.0,color="k",ls=":",label="M_F0=10 g")
titles=["Retentate cF [mM]","Permeate cV [mM]","Permeate mass mV [g]","Cell mass mF(t) [g]"]
for j,tt in enumerate(titles):
    ax[j].set_title(tt); ax[j].set_xlabel("time [min]"); ax[j].legend(fontsize=7); ax[j].grid(alpha=0.3)
fig.suptitle(f"{rid} — forward sim at each channel's optimal B (Lp={Lp:.1f} fixed).  Retentate wants HIGH B; mass wants LOW B.",fontsize=12)
fig.tight_layout(rect=(0,0,1,0.96))
out="/private/tmp/claude-505/-Users-kkasturi-GitHub-keshav-dev-dynamic-diafiltration-pyomo/d20a5790-00e6-4df6-bd33-bf6c3f39a8e5/scratchpad/nacl_overlay.png"
fig.savefig(out,dpi=140); print("WROTE",out)
# also print the model mF range per case
for lab,B,sg,col in cases:
    s=sim(B,sg); t,mF=traj(s,"mF"); print(f"{lab}: model mF range {mF.min():.2f}-{mF.max():.2f} g")
