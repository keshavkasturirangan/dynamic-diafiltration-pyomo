import sys, json, numpy as np, io, contextlib
sys.path.insert(0,"/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1")
import refactored_ucb_library as lib
lib.NF270_CF_RESIDUAL_SCALE_FRACTION=0.02; lib.NF270_CF_RESIDUAL_FLOOR_MM=None
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
STUDY=lib.NF270_DEFAULT_ROOT.parent/"DATA3/results/paper_artifacts/nf270/bform_study"
OUT=STUDY/"Bsigma_coeffs"
CONC=["MC2.05.07.24_NaCl","MC4.07.11.24_SNaCl","MC5.07.23.24_NaCl","MC5.07.23.24_SNaCl"]
def load(rid,conv):
    fam=lib.NF270_RUN_REGISTRY[rid]
    ds=lib.loadxlsx(lib.NF270_DEFAULT_ROOT/fam["workbook"],sheet=fam["sheet"])["data_stru"]
    lib.apply_campaign_time_correction(ds,convention=conv); return ds
def refit(ds,seed):
    with contextlib.redirect_stdout(io.StringIO()):
        fit,_,_=lib.solve_model(ds,ds["mode"],theta=dict(seed),sim_opt=False,B_form="single",workflow_family="DATA3",nfe=150,multistart=False,solver_max_cpu_time=150)
    return fit["parameters"] if isinstance(fit,dict) else None
def sim(ds,theta):
    with contextlib.redirect_stdout(io.StringIO()):
        fit,s,_=lib.solve_model(ds,ds["mode"],theta=dict(theta),sim_opt=True,B_form="single",workflow_family="DATA3",nfe=120,solver_max_cpu_time=40)
    def tr(k):
        t=[];y=[]
        for i in sorted(s.keys()): t+=list(np.asarray(s[i]["time"],float)); y+=list(np.asarray(s[i][k],float))
        return np.array(t),np.array(y)
    return tr
# fit both conventions
res={}
print("=== fitted theta per convention ===")
for rid in CONC:
    seed=json.load(open(STUDY/rid/"result_single.json"))["parameters"]
    thetas={}
    for conv in ("vial_close","tube_transit"):
        ds=load(rid,conv)
        if conv=="vial_close":
            th=seed   # canonical DATA2 vial-close fit already
        else:
            th=refit(ds,seed) or seed
        thetas[conv]=th
    res[rid]=thetas
    v=thetas["vial_close"]; t=thetas["tube_transit"]
    print(f"  {rid:20s} vial_close: Lp={v['Lp']:.2f} B={v['B']:.2f} σ={v['sigma']:.3f} | tube_transit: Lp={t['Lp']:.2f} B={t['B']:.2f} σ={t['sigma']:.3f}")
json.dump({r:{c:{k:float(x) for k,x in th.items()} for c,th in d.items()} for r,d in res.items()},
          open("/private/tmp/claude-505/-Users-kkasturi-GitHub-keshav-dev-dynamic-diafiltration-pyomo/d20a5790-00e6-4df6-bd33-bf6c3f39a8e5/scratchpad/tc_theta.json","w"),indent=2)

def make_fig(conv, title, fname):
    fig,axes=plt.subplots(2,2,figsize=(15,9.5)); axes=axes.flatten()
    for ax,rid in zip(axes,CONC):
        ds=load(rid,conv); th=res[rid][conv]; tr=sim(ds,th)
        t,mF=tr("mF"); _,mV=tr("mV"); _,cIn=tr("cIn"); _,cF=tr("cF")
        t0=t[0]; tm=(t-t0)/60
        mVc=mV.copy(); off=0.0; last=0.0
        for k in range(len(mVc)):
            if mV[k]<last-1e-6: off+=last
            last=mV[k]; mVc[k]=mV[k]+off
        ax.plot(tm,mF,color="tab:red",lw=1.8,label="retentate cell mass m$_F$(t) [g]")
        ax.plot(tm,mVc,color="tab:blue",lw=1.8,label="permeate collected m$_V$(t) [g]")
        ax.set_xlabel("time [min]"); ax.set_ylabel("fluid mass [g]"); ax.grid(alpha=0.3)
        ax2=ax.twinx()
        ax2.plot(tm,cIn,color="tab:green",lw=2.0,label="interfacial conc c$_{in}$(t) [mM]")
        ax2.plot(tm,cF,color="0.4",lw=1.0,ls="--",label="bulk retentate c$_F$(t) [mM]")
        ax2.set_ylabel("concentration [mM]",color="tab:green"); ax2.tick_params(axis="y",labelcolor="tab:green")
        h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
        ax.legend(h1+h2,l1+l2,fontsize=8,loc="lower right",framealpha=0.9)
        ax.set_title(f"{rid}   fit: Lp={th['Lp']:.1f}, B={th['B']:.1f}, σ={th['sigma']:.3f}",fontsize=10)
    fig.suptitle(title,fontsize=13); fig.tight_layout(rect=(0,0,1,0.96))
    fig.savefig(str(OUT)+"/"+fname,dpi=145); plt.close(fig); print("WROTE",fname)
make_fig("vial_close","Concentrating NaCl — cell mass / permeate mass / interfacial conc — WITH DATA2 vial-close time correction (default)","concNaCl_cellconc_WITH_vialclose.png")
make_fig("tube_transit","Concentrating NaCl — cell mass / permeate mass / interfacial conc — tube-transit (V_tube) lag time correction","concNaCl_cellconc_tubetransit.png")
print("DONE")
