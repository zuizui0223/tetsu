"""Render aggregate-only v314 comparison and machine-readable panel tables."""
from pathlib import Path
import argparse
import csv
import json


LABELS={"SI":"Slovenia", "HRBA":"Croatia + Bosnia", "SK":"Slovakia",
        "IN":"Kashmir, India", "AU_WA":"Western Australia*"}
ORDER=tuple(LABELS)


def write_csv(path,rows):
    keys=list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=keys,lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v
                             for k,v in row.items()})


def render(source,out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(out_dir/"v314_matched_foreign_rail_road_panels.csv",source["panels"])
    contrasts=[]
    for row in source["background_contrasts"]:
        flat={k:v for k,v in row.items() if k not in {"rail","background"}}
        for habitat in ("rail","background"):
            flat.update({habitat+"_"+key:value for key,value in row[habitat].items()})
        contrasts.append(flat)
    write_csv(out_dir/"v314_matched_foreign_rail_road_backgrounds.csv",contrasts)
    panels={(row["study"],row["recorded_only"]):row for row in source["panels"]
            if row["pool"]=="shared" and row["donor_mode"]=="both"}
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,
                         "axes.spines.right":False,"axes.titleweight":"bold"})
    fig,(ax,coverage)=plt.subplots(1,2,figsize=(14,7.2),gridspec_kw={"width_ratios":[1.45,1]})
    fig.subplots_adjust(left=.17,right=.97,top=.77,bottom=.22,wspace=.38)
    fig.text(.04,.95,"Do foreign railway frequencies transfer better than roadside frequencies?",
             fontsize=18,weight="bold",ha="left")
    fig.text(.04,.902,"Fixed Czech + Italian donors  |  Equal country weights  |  No recipient fitting",
             color="#475569",fontsize=12)
    colors={False:"#176B87",True:"#CE742D"}
    offsets={False:-.12,True:.12}
    finite=[row["delta_rho"] for row in panels.values() if row["estimable"]]
    lo=min([0]+finite)-.08
    hi=max([0]+finite)+.08
    ax.axvline(0,color="#64748B",lw=1)
    for i,study in enumerate(ORDER):
        for recorded in (False,True):
            row=panels[study,recorded]
            if row["estimable"]:
                ax.scatter(row["delta_rho"],i+offsets[recorded],s=64,
                           marker="o" if not recorded else "D",c=colors[recorded],zorder=3)
            else:
                ax.text(hi-.01,i+offsets[recorded],"not estimable",ha="right",va="center",
                        color=colors[recorded],fontsize=9)
        a,b=panels[study,False],panels[study,True]
        if a["estimable"] and b["estimable"]:
            ax.plot([a["delta_rho"],b["delta_rho"]],[i+offsets[False],i+offsets[True]],
                    color="#CBD5E1",lw=1.2,zorder=1)
        cov=a["named_taxon_coverage"] or 0
        coverage.barh(i,100*cov,height=.48,color="#91A9B8" if study!="AU_WA" else "#BEC7CF")
        coverage.text(100*cov+2,i,f'{a["n_recorded"]}/{a["all_target_named_taxa"]}  ({cov:.1%})',
                      va="center",fontsize=10)
    for axis in (ax,coverage):
        axis.axhline(3.5,color="#94A3B8",lw=1,ls=(0,(3,3)))
        axis.set_ylim(4.55,-.55)
        axis.set_yticks(range(len(ORDER)))
        axis.grid(axis="x",color="#E2E8F0",lw=.8)
        axis.set_axisbelow(True)
    ax.set_yticklabels([LABELS[study] for study in ORDER])
    coverage.set_yticklabels([])
    coverage.tick_params(axis="y",length=0)
    ax.set_xlim(lo,hi)
    coverage.set_xlim(0,110)
    coverage.set_xticks([0,25,50,75,100],["0","25","50","75","100%"])
    ax.set_title("A. Railway advantage in rank correlation",loc="left",fontsize=12,pad=17)
    ax.set_xlabel(r"$\Delta\rho = \rho_{rail\ donor} - \rho_{road\ donor}$",labelpad=10)
    coverage.set_title("B. Named target flora covered",loc="left",fontsize=12,pad=17)
    coverage.set_xlabel("Recorded target taxa in the shared donor pool",labelpad=10)
    handles=[plt.Line2D([],[],ls="",marker="o",color=colors[False],markersize=7,
                        label="All shared donor taxa, including target nondetections"),
             plt.Line2D([],[],ls="",marker="D",color=colors[True],markersize=7,
                        label="Only taxa recorded in the target survey")]
    fig.legend(handles=handles,loc="lower left",bbox_to_anchor=(.035,.075),frameon=False,
               fontsize=10,ncol=1)
    summary=source["study_summaries"][0]
    mean=summary["mean_delta_rho"]
    text="Not estimable" if mean is None else f'{mean:+.3f}'
    fig.text(.55,.128,f'Primary study mean: {text}  |  Positive: {summary["n_positive"]}/{summary["n_studies"]}',
             fontsize=11,weight="bold",color="#334155")
    fig.text(.04,.035,"*2500 m² railway-reserve quadrats; separate extrapolation, excluded from the primary mean.\n"
             "Unknown taxon keys are excluded. These are descriptive rank contrasts, with no causal interpretation or confidence intervals.",
             fontsize=9,color="#475569",va="bottom")
    path=out_dir/"v314_matched_foreign_rail_road_comparison.png"
    fig.savefig(path,dpi=190,facecolor="white")
    fig.savefig(path.with_suffix(".svg"),facecolor="white")
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    render(json.loads(args.input.read_text()),args.output_dir)


if __name__=="__main__":
    main()
