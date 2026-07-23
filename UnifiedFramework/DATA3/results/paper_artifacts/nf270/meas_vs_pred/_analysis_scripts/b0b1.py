import sys, json, numpy as np, io, contextlib, time
sys.path.insert(0,"/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/refactored_codes_v1")
import refactored_ucb_library as lib
lib.NF270_CF_RESIDUAL_SCALE_FRACTION=0.02; lib.NF270_CF_RESIDUAL_FLOOR_MM=None
lib.NF270_BFORM_SLICE_UNBOUNDED_B=True
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
STUDY=lib.NF270_DEFAULT_ROOT.parent/"DATA3/results/paper_artifacts/nf270/bform_study"
OUT=STUDY/"Bsigma_coeffs"
CONC=["MC2.05.07.24_NaCl","MC4.07.11.24_SNaCl","MC5.07.23.24_NaCl","MC5.07.23.24_SNaCl"]
def load(rid):
    fam=lib.NF270_RUN_REGISTRY[rid]; return lib.loadxlsx(lib.NF270_DEFAULT_ROOT/fam["workbook"],sheet=fam["sheet"])["data_stru"]
D=6
fig,axes=plt.subplots(2,2,figsize=(14,11)); axes=axes.flatten()
for ax,rid in zip(axes,CONC):
    ds=load(rid); mode=ds["mode"]
    p=json.load(open(STUDY/rid/"result_poly1.json"))["parameters"]
    Lp=p["Lp"]; sg=p["sigma"]; b0f=p["beta_0"]; b1f=p["beta_1"]; S0=p["S0"]; S=p["S"]
    # warm cIn for pre-filter
    with contextlib.redirect_stdout(io.StringIO()):
        _,sw,_=lib.solve_model(ds,mode,theta=dict(Lp=Lp,beta_0=b0f,beta_1=b1f,sigma=sg,S0=S0,S=S),sim_opt=True,B_form=1,workflow_family="DATA3",nfe=120,solver_max_cpu_time=18)
    cIn=np.concatenate([np.asarray(sw[i]["cIn"],float).reshape(-1) for i in sw]); cIn=cIn[np.isfinite(cIn)]
    b0ax=np.linspace(max(1e-3,b0f-max(abs(b0f),0.5)),b0f+max(abs(b0f),0.5),D)
    b1ax=np.linspace(b1f-max(abs(b1f),0.03),b1f+max(abs(b1f),0.03),D)
    Z=np.full((D,D),np.nan)
    for j,b1 in enumerate(b1ax):
        for i,b0 in enumerate(b0ax):
            Bc=b0+b1*cIn
            if not (Bc.min()>1e-4 and Bc.max()<250): continue
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    fit,s,_=lib.solve_model(ds,mode,theta=dict(Lp=Lp,beta_0=float(b0),beta_1=float(b1),sigma=sg,S0=S0,S=S),sim_opt=True,B_form=1,workflow_family="DATA3",nfe=120,solver_max_cpu_time=6,skip_sim_init=True)
                if isinstance(fit,dict): Z[j,i]=fit["obj_m"]+fit["obj_cv"]+fit["obj_cr"]
            except Exception: pass
    X,Y=np.meshgrid(b0ax,b1ax)
    if np.isfinite(Z).sum()>=6:
        Zl=np.log10(np.clip(Z,1e-9,None)); cs=ax.contour(X,Y,Zl,levels=12,cmap="turbo",linewidths=1.1); ax.clabel(cs,inline=True,fontsize=6,fmt="%.1f")
        jm,im=np.unravel_index(np.nanargmin(Z),Z.shape); ax.plot(X[jm,im],Y[jm,im],"^",color="red",ms=11,mec="k",zorder=6)
    ax.plot(b0f,b1f,"o",mfc="none",mec="k",ms=10,mew=1.3,label="fitted (β₀,β₁)")
    ax.axhline(0,color="0.6",lw=0.6,ls=":")
    ax.set_xlabel("β₀  [µm/s]  (B at c→0)"); ax.set_ylabel("β₁  [µm/s/mM]  (slope of B vs c$_{in}$)")
    ax.set_title(f"{rid}\nlinear B(c)=β₀+β₁·c$_{{in}}$  ·  Lp={Lp:.1f}, σ={sg:.2f} fixed",fontsize=10); ax.legend(fontsize=8,loc="lower left",framealpha=0.9); ax.grid(alpha=0.3)
fig.suptitle("Concentrating NaCl — β₀ vs β₁ square-slice WSSE contour (linear B–interfacial-conc correlation)",fontsize=13)
fig.tight_layout(rect=(0,0,1,0.96))
fig.savefig(str(OUT)+"/concNaCl_beta0_vs_beta1.png",dpi=145); print("WROTE",str(OUT)+"/concNaCl_beta0_vs_beta1.png")
