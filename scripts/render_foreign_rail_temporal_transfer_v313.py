"""Render aggregate-only v313 figures; no raw species/site output."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,default=ROOT/"data/derived/v313_foreign_rail_temporal_transfer.json")
    parser.add_argument("--output-prefix",type=Path,default=ROOT/"docs/manuscript_assets/v313_foreign_rail_temporal_transfer")
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    panels = {(p["endpoint"],p["penalty"]):p for p in data["panels"]}
    habitats = ("Track","Interior Plot 1","Interior Plot 2","balanced")
    display = ("Track","Interior 1","Interior 2","Balanced")
    colors = ("#296c91","#4d8a70","#85a886","#273443")
    fig,axes = plt.subplots(2,2,figsize=(11.4,8.8))
    rows = []
    for col,(endpoint,title) in enumerate((("persistence","Recorded persistence"),("new_record","New records"))):
        p = panels[endpoint,.01]
        ax = axes[0,col]
        for j,h in enumerate(habitats):
            s = p["comparisons"]["foreign_common"][h]["log_loss"]
            ax.bar(j,s["mean"],color=colors[j],width=.6,zorder=2)
            ax.scatter(j+np.linspace(-.16,.16,5),s["block_values"],s=20,
                       facecolor="white",edgecolor=colors[j],linewidth=.9,zorder=3)
            for b,v in enumerate(s["block_values"]):
                rows.append({"endpoint":endpoint,"penalty":.01,"quantity":"M0_minus_M1_log_loss","scope":h,"block":b,"value":v})
        ax.axhline(0,color="#58616a",linewidth=.9)
        ax.set_xticks(range(4),display)
        ax.set_title(f"{'AB'[col]}. {title}: common foreign score",loc="left",fontsize=11,fontweight="bold")
        ax.set_ylabel("Held-out log-loss gain (nat)")
        ax.grid(axis="y",alpha=.18,zorder=0)
        ax = axes[1,col]
        penalties = (.001,.01,.1)
        for key,label,color in (
            ("Track_minus_Interior_Plot_1","Track minus Interior 1","#296c91"),
            ("Track_minus_Interior_Plot_2","Track minus Interior 2","#b76931")):
            values = [panels[endpoint,q]["foreign_slopes"][key]["mean"] for q in penalties]
            ax.plot(penalties,values,"o-",label=label,color=color,linewidth=1.5,markersize=5)
            for q in penalties:
                blockvalues=panels[endpoint,q]["foreign_slopes"][key]["block_values"]
                ax.scatter([q]*5,blockvalues,color=color,alpha=.35,s=16)
                for b,v in enumerate(blockvalues):
                    rows.append({"endpoint":endpoint,"penalty":q,"quantity":"M2_foreign_raw_slope_difference","scope":key,"block":b,"value":v})
        ax.axhline(0,color="#58616a",linewidth=.9)
        ax.set_xscale("log")
        ax.set_xticks(penalties,["0.001","0.01 (primary)","0.1"])
        ax.set_xlabel("Ridge penalty")
        ax.set_ylabel("Foreign-score slope difference")
        ax.set_title(f"{'CD'[col]}. {title}: Track contrast sensitivity",loc="left",fontsize=11,fontweight="bold")
        ax.grid(axis="y",alpha=.18)
        ax.legend(fontsize=8,frameon=False,loc="best")
    for ax in axes.flat:
        ax.spines[["top","right"]].set_visible(False)
        ax.tick_params(labelsize=9)
    fig.suptitle("Foreign railway frequency and Kashmir recorded transitions, 2014–2017",
                 fontsize=15,fontweight="bold",y=.98)
    fig.text(.5,.933,"51 taxa known at baseline · 27 paired transects · 5 spatial holdout blocks",ha="center",fontsize=10)
    fig.text(.08,.038,"A–B: bars = equal-block means; dots = block gains; positive values mean better predictions. Primary λ = 0.01.\n"
             "C–D: lines = mean fitted slope contrasts; dots = overlapping-training-fold estimates, not independent replicates.\n"
             "Retrospective associations. Balanced loss weights Track 1/2 and each interior 1/4; no survival or transport inference.",
             fontsize=9,linespacing=1.5,color="#3c454d")
    fig.subplots_adjust(left=.08,right=.97,top=.87,bottom=.17,hspace=.42,wspace=.28)
    args.output_prefix.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output_prefix.with_suffix(".png"),dpi=200,facecolor="white")
    fig.savefig(args.output_prefix.with_suffix(".svg"),facecolor="white")
    with args.output_prefix.with_suffix(".csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


if __name__ == "__main__":
    main()

