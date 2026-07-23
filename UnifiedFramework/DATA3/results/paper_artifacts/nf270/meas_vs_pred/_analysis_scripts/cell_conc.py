import sys, json, numpy as np, io, contextlib
sys.path.insert(0,"/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1")
import refactored_ucb_library as lib
lib.NF270_CF_RESIDUAL_SCALE_FRACTION=0.02; lib.NF270_CF_RESIDUAL_FLOOR_MM=None
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
STUDY=lib.NF270_DEFAULT_ROOT.parent/"DATA3/results/paper_artifacts/nf270/bform_study"
CONC=["MC2.05.07.24_NaCl","MC4.07.11.24_SNaCl","MC5.07.23.24_NaCl","MC5.07.23.24_SNaCl"]
def load(rid):
    fam=lib.NF270_RUN_REGISTRY[rid]; return lib.loadxlsx(lib.NF270_DEFAULT_ROOT/fam["workbook"],sheet=fam["sheet"])["data_stru"]
def traj(s,key):
    t=[];y=[]
    for i in sorted(s.keys()): t+=list(np.asarray(s[i]["time"],float)); y+=list(np.asarray(s[i][key],float))
    return np.array(t),np.array(y)
fig,axes=plt.subplots(2,2,figsize=(15,9.5)); axes=axes.flatten()
for ax,rid in zip(axes,CONC):
    ds=load(rid); mode=ds["mode"]
    p=json.load(open(STUDY/rid/"result_single.json"))["parameters"]
    with contextlib.redirect_stdout(io.StringIO()):
        fit,s,_=lib.solve_model(ds,mode,theta=dict(p),sim_opt=True,B_form="single",workflow_family="DATA3",nfe=120,solver_max_cpu_time=40)
    t,mF=traj(s,"mF"); _,mV=traj(s,"mV"); _,cIn=traj(s,"cIn"); _,cF=traj(s,"cF"); _,cH=traj(s,"cH")
    t0=t[0]; tm=(t-t0)/60
    # cumulative permeate: mV resets per vial -> accumulate
    mVc=mV.copy(); off=0.0; last=0.0
    for k in range(len(mVc)):
        if mVc[k]<last-1e-6: off+=last  # reset detected
        last=mV[k]; mVc[k]=mV[k]+off
    ax.plot(tm,mF,color="tab:red",lw=1.8,label="retentate cell mass  m$_F$(t) [g]")
    ax.plot(tm,mVc,color="tab:blue",lw=1.8,label="permeate collected  m$_V$(t) [g]")
    ax.set_xlabel("time [min]"); ax.set_ylabel("fluid mass [g]"); ax.grid(alpha=0.3)
    ax2=ax.twinx()
    ax2.plot(tm,cIn,color="tab:green",lw=2.0,label="interfacial conc  c$_{in}$(t) [mM]")
    ax2.plot(tm,cF,color="0.4",lw=1.0,ls="--",label="bulk retentate c$_F$(t) [mM]")
    ax2.set_ylabel("concentration [mM]",color="tab:green"); ax2.tick_params(axis="y",labelcolor="tab:green")
    h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
    ax.legend(h1+h2,l1+l2,fontsize=8,loc="lower right",framealpha=0.9)
    pol=cIn/np.clip(cF,1e-6,None)
    ax.set_title(f"{rid}  (Lp={p['Lp']:.1f}, B={p['B']:.1f}, σ={p['sigma']:.2f})\n"
                 f"conc.-polarization c$_{{in}}$/c$_F$ ≈ {np.nanmedian(pol[len(pol)//3:]):.2f}", fontsize=10)
fig.suptitle("Concentrating NaCl — stirred-cell (retentate) mass, permeate mass, and interfacial concentration vs time",fontsize=13)
fig.tight_layout(rect=(0,0,1,0.96))
OUT="/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/paper_artifacts/nf270/bform_study/Bsigma_coeffs"
fig.savefig(OUT+"/concNaCl_cell_mass_interfacial_conc.png",dpi=145); print("WROTE",OUT+"/concNaCl_cell_mass_interfacial_conc.png")
