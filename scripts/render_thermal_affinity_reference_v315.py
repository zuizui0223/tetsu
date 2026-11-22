"""Render the distinct thermal estimands without implying confidence intervals."""
from pathlib import Path
import argparse
import csv
import json

REFERENCES=[("CZ","road:CZ","Czechia / mountain roads"),
            ("IT","road:IT","Italy / mountain roads"),
            ("SI","background:SI:urban","Slovenia / tree pits"),
            ("HRBA","background:HRBA:urban","Croatia–Bosnia / tree pits"),
            ("SK","background:SK:urban","Slovakia / urban plots"),
            ("IN","background:IN:interior","Kashmir / paired interiors"),
            ("IN","road:IN","Kashmir / MIREN roads")]
RAIL_LABELS={"CZ":"Czechia","IT":"Italy","SI":"Slovenia","HRBA":"Croatia–Bosnia",
             "SK":"Slovakia","IN":"Kashmir*","AU_WA":"Western Australia*"}


def write_csv(path,rows):
    columns=list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,lineterminator="\n")
        writer.writeheader()
        writer.writerows({key:json.dumps(value,ensure_ascii=False) if isinstance(value,(dict,list)) else value
                         for key,value in row.items()} for row in rows)


def render(result,directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    directory.mkdir(parents=True,exist_ok=True)
    rows=[]
    for row in result["reference_contrasts"]:
        flat={key:value for key,value in row.items() if key not in {"rail","reference"}}
        for side in ("rail","reference"):
            flat.update({side+"_"+key:value for key,value in row[side].items()})
        rows.append(flat)
    write_csv(directory/"v315_thermal_reference_contrasts.csv",rows)
    write_csv(directory/"v315_thermal_assemblages.csv",result["assemblages"])
    write_csv(directory/"v315_eive_unknown_score_envelopes.csv",result["eive_unknown_envelopes"])
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
                         "axes.spines.right":False,"axes.titleweight":"bold"})
    fig,axes=plt.subplots(1,3,figsize=(16,8),gridspec_kw={"width_ratios":[1,1,1.13]})
    fig.subplots_adjust(left=.235,right=.98,top=.76,bottom=.20,wspace=.32)
    fig.text(.035,.95,"Thermal affinity of railway floras and regional references",fontsize=19,weight="bold")
    fig.text(.035,.905,"Descriptive comparisons  |  Five European studies and one external region  |  Indicator units, not °C",
             fontsize=12,color="#475569")
    index={(row["study"],row["reference_dataset"],row["scenario"]):row for row in result["reference_contrasts"]}
    titles=["A. Definite-affinity subset\nTichý temperature scale 1–12",
            "B. Broader indicator definition\nEIVE temperature scale 0–10",
            "C. Hypothetical all-named EIVE\nUnknown species assigned 0–10"]
    for axis,title in zip(axes,titles):
        axis.set_title(title,loc="left",fontsize=11,pad=14)
        axis.axvline(0,color="#64748B",lw=1)
        axis.axhline(4.5,color="#94A3B8",ls=(0,(3,3)),lw=1)
        axis.set_ylim(6.65,-.65)
        axis.set_yticks(range(7))
        axis.grid(axis="x",color="#E2E8F0",lw=.8)
        axis.set_axisbelow(True)
    for axis,main,common,color in zip(axes[:2],("tichy_numeric","eive_all"),
                ("tichy_common_numeric","eive_common_numeric"),("#176B87","#A65337")):
        values=[]
        for i,(study,reference,_) in enumerate(REFERENCES):
            for scenario,offset,filled in ((main,-.10,True),(common,.10,False)):
                row=index[study,reference,scenario]
                if row["estimable"]:
                    value=row["difference"]
                    values.append(value)
                    axis.scatter(value,i+offset,s=57 if filled else 43,facecolor=color if filled else "white",
                                 edgecolor=color,lw=1.2,zorder=3)
                else:
                    axis.text(.02,i+offset,"NA",va="center",fontsize=9,color=color)
        low,high=min([0]+values),max([0]+values)
        pad=max(.12,(high-low)*.18)
        axis.set_xlim(low-pad,high+pad)
        axis.set_xlabel("Railway − reference mean",labelpad=10)
    envelopes={(row["study"],row["reference_dataset"]):row for row in result["eive_unknown_envelopes"]}
    endpoints=[]
    for i,(study,reference,_) in enumerate(REFERENCES):
        row=envelopes[study,reference]
        if not row["estimable"]: continue
        low,high=row["lower"],row["upper"]
        color="#176B87" if row["sign"]=="positive" else "#A65337" if row["sign"]=="negative" else "#64748B"
        axes[2].plot([low,high],[i,i],color=color,lw=3,solid_capstyle="butt")
        axes[2].plot([low,high],[i,i],ls="",marker="|",color=color,ms=11)
        endpoints.extend((low,high))
    low,high=min([0]+endpoints),max([0]+endpoints)
    pad=max(.12,(high-low)*.12)
    axes[2].set_xlim(low-pad,high+pad)
    axes[2].set_xlabel("Admissible full-score difference",labelpad=10)
    axes[0].set_yticklabels([entry[2] for entry in REFERENCES])
    for axis in axes[1:]:
        axis.set_yticklabels([])
        axis.tick_params(axis="y",length=0)
    handles=[plt.Line2D([],[],ls="",marker="o",color="#475569",markersize=7,label="Each system's numeric taxa"),
             plt.Line2D([],[],ls="",marker="o",color="#475569",markerfacecolor="white",markersize=7,
                        label="Same trait-eligible species universe in both systems")]
    fig.legend(handles=handles,loc="lower left",bbox_to_anchor=(.03,.09),frameon=False,ncol=2,fontsize=10)
    fig.text(.035,.035,"Panel C is a separate hypothetical estimand, not an error bar or a confidence interval for panel B.\n"
             "Reference habitats, climate distributions and sampling designs differ. No between-region climate interaction is estimated.",
             fontsize=10,color="#475569")
    for extension in ("png","svg"):
        fig.savefig(directory/("v315_thermal_reference_comparison."+extension),dpi=180,facecolor="white")
    plt.close(fig)

    # Indifference is an observed source category, not an imputed temperature.
    primary={row["dataset"]:row for row in result["assemblages"] if row["scenario"]=="tichy_numeric"}
    fig,axis=plt.subplots(figsize=(12,6.7))
    fig.subplots_adjust(left=.22,right=.95,top=.78,bottom=.22)
    fig.text(.035,.94,"What the primary temperature index represents",fontsize=19,weight="bold")
    fig.text(.035,.89,"Fractions of identifiable railway species-incidence mass  |  Tichý et al. dataset v2",fontsize=12,color="#475569")
    categories=[("numeric","Numeric temperature","#176B87"),
                ("indifferent","Temperature-indifferent (x)","#D5A142"),
                ("other","Unmatched or unavailable","#CBD5E1")]
    left=np.zeros(len(RAIL_LABELS))
    for category,label,color in categories:
        values=[]
        for key in RAIL_LABELS:
            row=primary["rail:"+key]
            fractions=row["status_frequency_fraction"]
            value=fractions.get(category,0.) if category!="other" else 1-fractions.get("numeric",0.)-fractions.get("indifferent",0.)
            values.append(100*value)
        axis.barh(range(len(values)),values,left=left,color=color,height=.65,label=label)
        left+=values
    for i,key in enumerate(RAIL_LABELS):
        row=primary["rail:"+key]
        value=100*(row["numeric_frequency_coverage"] or 0.)
        if value>8:
            axis.text(value/2,i,f"{value:.1f}%",ha="center",va="center",color="white",fontsize=10,weight="bold")
        else:
            axis.text(101,i,f"{value:.1f}% numeric",ha="left",va="center",fontsize=9,color="#475569")
    axis.set_yticks(range(len(RAIL_LABELS)),list(RAIL_LABELS.values()))
    axis.invert_yaxis()
    axis.axhline(4.5,color="#94A3B8",ls=(0,(3,3)),lw=1)
    axis.set_xlim(0,116)
    axis.set_xticks([0,25,50,75,100],["0","25","50","75","100%"])
    axis.set_xlabel("Recorded-frequency mass within each railway inventory",labelpad=10)
    axis.legend(loc="upper center",bbox_to_anchor=(.48,-.18),frameon=False,ncol=3,fontsize=10)
    fig.text(.035,.04,"*Application outside the European indicator region. The mean excludes x and unavailable values; x is never assigned zero or a midpoint.\n"
             "Frequency weighting is not equal-site weighting; all frequencies refer to recorded named taxa, not detection-adjusted occupancy.",
             fontsize=9,color="#475569")
    for extension in ("png","svg"):
        fig.savefig(directory/("v315_thermal_trait_coverage."+extension),dpi=180,facecolor="white")
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    render(json.loads(args.input.read_text()),args.output_dir)


if __name__=="__main__":
    main()
