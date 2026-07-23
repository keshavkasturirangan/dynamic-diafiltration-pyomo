import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SCR="/private/tmp/claude-505/-Users-kkasturi-GitHub-keshav-dev-dynamic-diafiltration-pyomo/d20a5790-00e6-4df6-bd33-bf6c3f39a8e5/scratchpad"
bscan=json.load(open(SCR+"/bscan.json"))
bal=np.load(SCR+"/nacl_bal.npy",allow_pickle=True).item()
order=["MC2.05.07.24_NaCl","MC5.07.23.24_NaCl","MC5.07.23.24_SNaCl","MC4.07.11.24_SNaCl","MC3.07.22.24_SNaCl","MC5.07.23.24_S2NaCl"]
short=[s.replace("MC","").replace("_NaCl","·N").replace("_SNaCl","·SN").replace("_S2NaCl","·S2N") for s in order]
fig,ax=plt.subplots(1,3,figsize=(19,5.2))
# panel A: per-channel B optimum
x=np.arange(len(order)); w=0.25
Bm=[bscan[s]["Bm"] for s in order]; Bp=[bscan[s]["Bp"] for s in order]; Br=[bscan[s]["Br"] for s in order]
ax[0].bar(x-w,Bm,w,label="mass",color="tab:green"); ax[0].bar(x,Bp,w,label="permeate",color="tab:blue"); ax[0].bar(x+w,Br,w,label="retentate",color="tab:red")
ax[0].set_xticks(x); ax[0].set_xticklabels(short,rotation=40,ha="right",fontsize=8)
ax[0].set_ylabel("B at channel optimum [µm/s] (σ=1, Lp*)"); ax[0].set_title("(a) Per-channel best-fit B — retentate scatters / rails\n(bars at 1 or 30 = railing to grid edge)")
ax[0].legend(); ax[0].grid(alpha=0.3,axis="y")
for xi,s in zip(x,order):
    d="CONC" if bscan[s]["dir"]=="CONC" else "DIL"
    ax[0].text(xi,-3.0,d,ha="center",fontsize=7,color=("0.3" if d=="CONC" else "tab:purple"))
# panel B: retentate-permeate deviation vs throughput
thru=[bscan[s]["thru"] for s in order]; dev=[bscan[s]["Br"]-bscan[s]["Bp"] for s in order]
cols=["tab:red" if bscan[s]["dir"]=="CONC" else "tab:purple" for s in order]
ax[1].axhline(0,color="k",lw=0.8)
for xi,yi,c,s in zip(thru,dev,cols,short):
    ax[1].scatter(xi,yi,c=c,s=90); ax[1].annotate(s,(xi,yi),fontsize=7,xytext=(4,4),textcoords="offset points")
ax[1].set_xlabel("throughput = permeate collected / M_F0  (feed volumes)")
ax[1].set_ylabel("B_retentate − B_permeate  [µm/s]")
ax[1].set_title("(b) Retentate deviation vs feed throughput\nNO clean correlation → not explained by feed accumulation")
ax[1].grid(alpha=0.3)
# panel C: solute-balance closure + effective cell mass
Rpct=[]; ratio=[]
for s in order:
    b=bal[s]; Rpct.append(b["mF_end_solute"]) # actually store ratio
    ratio.append(b["mF_end_solute"]/b["M_F0"])
ax[2].axhline(1.0,color="k",ls=":",label="M_F0 (no accumulation)")
ax[2].bar(x,ratio,0.5,color=["tab:red" if bscan[order[i]]["dir"]=="CONC" else "tab:purple" for i in range(len(order))])
ax[2].set_xticks(x); ax[2].set_xticklabels(short,rotation=40,ha="right",fontsize=8)
ax[2].set_ylabel("effective cell mass / M_F0  (from solute balance)")
ax[2].set_title("(c) Cell mass implied by solute balance\n≈1 for most; scatter reflects sensitivity, not systematic accumulation")
ax[2].axhline(1.0,color="k",ls=":"); ax[2].grid(alpha=0.3,axis="y"); ax[2].set_ylim(0,1.8)
fig.suptitle("NaCl single-salt: is the retentate-contour deviation caused by non-constant stirred-cell mass (feed accumulation)?",fontsize=13)
fig.tight_layout(rect=(0,0,1,0.95))
out=SCR+"/nacl_summary.png"; fig.savefig(out,dpi=140); print("WROTE",out)
