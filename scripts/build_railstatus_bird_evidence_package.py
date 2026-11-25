#!/usr/bin/env python3
"""Build publication evidence tables from verified frozen artifacts.
No additional model fitting; retain unsupported panels and explicit limits.
"""
import argparse,csv,hashlib,json,math
from pathlib import Path
RUNS={"temporal":45579119530,"guild":45583858584,"conditional":45584474048,"historical":45590411848}
ARTIFACTS={"temporal":15004679743,"guild":15007417808,"conditional":15008459726,"historical":15010704848}
VERSIONS={"temporal":"railstatus-bird-temporal-joint-influence-v1","guild":"railstatus-avonet-guild-temporal-audit-v1","conditional":"railstatus-avonet-conditional-hwi-connectivity-v1","historical":"railstatus-avonet-heldback-2018-22-hwi-connectivity-v1"}
SHA="a8f0e0c889a05485f306029433cd45068836324e2bc6cf5b1ed8e34417f3d5fe"
def read(p,key):
    raw=p.read_bytes();z=json.loads(raw)
    if z["version"]!=VERSIONS[key]:raise ValueError("Wrong artifact version: "+key)
    if sorted(x["land"] for x in z["panels"])!=[.5,.8]:raise ValueError("Wrong land support")
    return {x["land"]:x for x in z["panels"]},{"run":RUNS[key],"artifact":ARTIFACTS[key],"json_SHA256":hashlib.sha256(raw).hexdigest(),"version":z["version"]}
def row(key,label,land,model,n,term,estimate_group="descriptive_all_years",sign=1,period="2023-25"):
    q={"source":key,"label":label,"land":land,"n_pairs":n,"period":period,"term":term,
       "status":model["status"],"estimate":None,"lower":None,"upper":None,"n_blocks":None}
    if model["status"]!="estimated":return q
    if model["n_pairs"]!=n:raise ValueError("Unexpected model support: "+label)
    e=model[estimate_group][term];b=e["estimate"];l,h=e["interval95_approx"]
    if not all(math.isfinite(x) for x in (b,l,h)) or not l<=b<=h:raise ValueError("Invalid interval")
    q.update(estimate=sign*b,lower=l if sign==1 else -h,upper=h if sign==1 else -l,n_blocks=model["n_blocks"])
    return q
def check(q,estimate):
    if q["status"]!="estimated" or abs(q["estimate"]-estimate)>5e-6:raise ValueError("Previously verified estimate not replayed: "+q["label"])
def build(t,g,c):
    out=[]
    for land in (.5,.8):
        n={.5:819,.8:649}[land]
        if t[land]["n_bird_pairs"]!=n:raise ValueError("Frozen bird geography changed")
        for m,desc in (("bird_endpoint","Bird: endpoint FE"),("bird_block","Bird: block FE")):
            for term in ("active","active_minus_inactive"):
                out.append(row("temporal",desc,land,t[land][m],n,term,"validation_2023_2024",period="2023-24"))
        z=g[land]["axes"]["HWI"];ng={.5:537,.8:424}[land]
        if z["n_all_year_common"]!=ng:raise ValueError("Frozen HWI support changed")
        for mode in ("block_FE","endpoint_FE"):
            audit=z["audit"][mode]
            for target,desc,sign in (("g0_jaccard","Low HWI: Jaccard",1),("g1_jaccard","High HWI: Jaccard",1),("difference_jaccard","Low minus high HWI: Jaccard",-1),("difference_simpson","Low minus high HWI: Simpson",-1)):
                for term in ("active","active_minus_inactive"):
                    out.append(row("guild",desc+" "+mode,land,audit[target],ng,term,sign=sign))
        for ax,label in (("habitat","Open minus dense"),("migration","Migrant minus resident")):
            v=g[land]["axes"][ax];nn=v["n_all_year_common"]
            out.append(row("guild",label+" block FE",land,v["audit"]["block_FE"]["difference_jaccard"],nn,"active"))
        for lab,nn in (("habitat_x_migration",{.5:217,.8:166}[land]),("resident_only",{.5:159,.8:123}[land]),("dense_resident",{.5:33,.8:23}[land])):
            a=c[land]["conditioned"][lab]
            if a["n_all_year_common"]!=nn:raise ValueError("Conditional support changed")
            for mode in ("block_FE","endpoint_FE"):
                q=row("conditional","Conditional HWI: "+lab+" "+mode,land,a[mode]["difference_jaccard"],nn,"active",sign=-1)
                if q["status"]=="estimated":raise ValueError("Previously unidentified panel changed")
                out.append(q)
    def get(label,land,term):
        v=[x for x in out if (x["label"],x["land"],x["term"])==(label,land,term)]
        if len(v)!=1:raise ValueError("Evidence selection ambiguous")
        return v[0]
    check(get("Bird: endpoint FE",.5,"active"),.016763820323536363)
    check(get("Low minus high HWI: Jaccard block_FE",.5,"active"),.02998709167574483)
    check(get("Low minus high HWI: Jaccard block_FE",.8,"active"),.0128861824804614)
    return out
def historical_row(label,land,model,n,term,sign=-1):
    q={"source":"historical","label":label,"land":land,"n_pairs":n,"period":"2018-22 pooled","term":term,
       "status":model["status"],"estimate":None,"lower":None,"upper":None,"n_blocks":model.get("n_blocks")}
    if model.get("n_pairs")!=n:raise ValueError("Historical support changed: "+label)
    if model["status"]!="estimated":return q
    e=model["effects"]["cross_active_normal"] if term=="active" else model["active_minus_inactive"]
    if e is None:
        q["status"]="contrast not identifiable";return q
    b=e["estimate"];l,h=e["interval95_approx"]
    if not all(math.isfinite(x) for x in (b,l,h)) or not l<=b<=h:raise ValueError("Invalid historical interval")
    q.update(estimate=sign*b,lower=l if sign==1 else -h,upper=h if sign==1 else -l)
    return q

def historical_rows(h):
    out=[]
    for land in (.5,.8):
        if h[land]["n_full_historical"]!={.5:811,.8:641}[land]:raise ValueError("Frozen historical support changed")
        for key,lab,n in (("historical_native","Early native",{.5:693,.8:551}[land]),
                          ("historical_on_original_hwi_geography","Early discovery-geography",{.5:535,.8:422}[land])):
            v=h[land][key]
            if v["n_raw"]!=n:raise ValueError("Historical HWI support changed")
            for mode,models in v["unadjusted"].items():
                for target,title,sign in (("g0_jaccard","Low HWI Jaccard",1),("g1_jaccard","High HWI Jaccard",1),("difference_jaccard","Low minus high HWI Jaccard",-1),("difference_simpson","Low minus high HWI Simpson",-1)):
                    for term in ("active","active_minus_inactive"):
                        out.append(historical_row(lab+" / "+title+" / "+mode,land,models[target],n,term,sign))
            for cond,c in v["conditioned"].items():
                for mode,models in c["models"].items():
                    for target,title in (("difference_jaccard","Conditional Jaccard"),("difference_simpson","Conditional Simpson"),("raw_difference_jaccard","Same-design raw Jaccard"),("standardized_minus_raw_jaccard","Conditional minus raw Jaccard")):
                        for term in ("active","active_minus_inactive"):
                            out.append(historical_row(lab+" / "+cond+" / "+title+" / "+mode,land,models[target],c["n"],term))
    def get(lab,land,term):
        v=[r for r in out if (r["label"],r["land"],r["term"])==(lab,land,term)]
        if len(v)!=1:raise ValueError("Ambiguous historical row")
        return v[0]
    check(get("Early native / Low minus high HWI Jaccard / block_FE",.5,"active"),.03170)
    check(get("Early native / Low minus high HWI Jaccard / endpoint_FE",.5,"active"),-.01841)
    check(get("Early native / habitat_x_migration / Conditional Jaccard / block_FE",.5,"active"),.00431)
    check(get("Early native / habitat_x_migration / Same-design raw Jaccard / block_FE",.5,"active"),.01899)
    check(get("Early discovery-geography / habitat_x_migration / Conditional Jaccard / block_FE",.5,"active"),-.00997)
    return out

def plot_historical(rows,root):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    wanted=[
        ("Early native / Low minus high HWI Jaccard / block_FE", "All eligible: block FE"),
        ("Early native / habitat_x_migration / Same-design raw Jaccard / block_FE", "Common support: unadjusted"),
        ("Early native / habitat_x_migration / Conditional Jaccard / block_FE", "Common support: standardized"),
        ("Early native / Low minus high HWI Jaccard / endpoint_FE", "All eligible: endpoint FE"),
        ("Early discovery-geography / Low minus high HWI Jaccard / block_FE","Discovery geography: block FE"),
        ("Early discovery-geography / habitat_x_migration / Conditional Jaccard / block_FE","Discovery geography: standardized")]
    fig,ax=plt.subplots(figsize=(7.4,4.2));labels=[]
    for i,(key,label) in enumerate(wanted):
        z=[r for r in rows if r["source"]=="historical" and r["land"]==.5 and r["label"]==key and r["term"]=="active"]
        if len(z)!=1 or z[0]["estimate"] is None:raise ValueError("Historical figure row unavailable")
        r=z[0];ax.errorbar(r["estimate"],i,xerr=[[r["estimate"]-r["lower"]],[r["upper"]-r["estimate"]]],fmt="o",capsize=3,color="#245D76" if i==0 else "#748896")
        labels.append(label+f' (n={r["n_pairs"]})')
    ax.set_yticks(range(len(labels)),labels);ax.invert_yaxis();ax.axvline(0,color="#777777",linestyle="--")
    ax.set_xlabel("Low-HWI minus high-HWI crossing coefficient")
    ax.set_title("Historical retest: association depends on adjustment",loc="left",fontsize=12)
    ax.spines[["right","top"]].set_visible(False);ax.grid(axis="x",alpha=.15)
    fig.text(.02,.015,"2018–22 pooled; land ≥50%; approximate 95% intervals. Conditional effect not established.",fontsize=8)
    fig.tight_layout(rect=(0,.05,1,1));fig.savefig(root/"figure_historical_hwi_falsification.svg",bbox_inches="tight");plt.close(fig)

def write_report(rows,meta,root):
    root.mkdir(parents=True,exist_ok=True)
    with (root/"bird_evidence.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (root/"source_receipts.json").write_text(json.dumps(meta,indent=2,allow_nan=False)+"\n")
    lines=["# Frozen bird connectivity evidence","",
           "Retrospective observational association, not measured movements. HWI is a morphology proxy. All-years trait contrasts are exploratory. 2023–24 bird overall is a retrospective heldback mean excluding the 2025 discovery year. Same platform and fixed 2026 mapped rail geometry. Approximate block-clustered intervals. No conditional effect is available when support gates fail.","",
           "| Source / contrast | Land fraction | Pairs | Term | Estimate [95% interval] |","|---|---:|---:|---|---|"]
    for r in rows:
        effect=(f'{r["estimate"]:+.5f} [{r["lower"]:+.5f}, {r["upper"]:+.5f}]' if r["estimate"] is not None else "Not identifiable")
        lines.append(f'| {r["label"]} ({r["period"]}) | {r["land"]:.0%} | {r["n_pairs"]} | {r["term"]} | {effect} |')
    lines+=["","Historical 2018–22 pooled retest: unadjusted block-FE HWI contrast recurs on native support, but the preregistered habitat/migration-standardized contrast is uncertain; endpoint-FE contrast reverses. Common-support attenuation itself is uncertain. No independent mobility mechanism established.","","AVONET: verified public derivative, not authenticated publisher archive. Source SHA and workflow receipts are recorded."]
    (root/"EVIDENCE.md").write_text("\n".join(lines)+"\n",encoding="utf8")
def plot(rows,root):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    select=[r for r in rows if r["source"]=="guild" and r["label"].startswith("Low minus high HWI:") and r["label"].endswith("block_FE") and r["estimate"] is not None]
    if len(select)!=8:raise ValueError("Frozen figure panel count changed")
    fig,ax=plt.subplots(figsize=(7.2,4.2))
    for i,r in enumerate(select):
        color="#245D76" if r["land"]==.5 else "#748896"
        ax.errorbar(r["estimate"],i,xerr=[[r["estimate"]-r["lower"]],[r["upper"]-r["estimate"]]],fmt="o",color=color,capsize=3,markersize=4,linewidth=1.4)
    labels=[f'{r["label"].split(": ")[1].split(" ")[0]} · land ≥ {r["land"]:.0%} · '+("active" if r["term"]=="active" else "active − inactive")+f' (n={r["n_pairs"]})' for r in select]
    ax.set_yticks(range(len(select)),labels);ax.invert_yaxis();ax.axvline(0,color="#777777",linestyle="--",linewidth=.9)
    ax.set_xlabel("Low-HWI minus high-HWI coefficient (2023–25 mean)")
    ax.set_title("Wing morphology and mapped railway connectivity",loc="left",fontsize=12)
    ax.spines[["right","top"]].set_visible(False);ax.grid(axis="x",alpha=.15)
    fig.text(.02,.015,"Exploratory block fixed effects; approximate 95% intervals; no movement or causal inference.",fontsize=8)
    fig.tight_layout(rect=(0,.045,1,1));fig.savefig(root/"figure_hwi_contrasts.svg",bbox_inches="tight");plt.close(fig)
def main():
    p=argparse.ArgumentParser()
    for k in (*RUNS,"out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();d={};meta={}
    for k in RUNS:d[k],meta[k]=read(getattr(a,k),k)
    rows=build(d["temporal"],d["guild"],d["conditional"])+historical_rows(d["historical"])
    write_report(rows,meta,a.out);plot(rows,a.out);plot_historical(rows,a.out)
    print(json.dumps({"evidence_rows":len(rows),"estimated":sum(r["status"]=="estimated" for r in rows),"not_identifiable":sum(r["status"]!="estimated" for r in rows),"source_receipts":meta},allow_nan=False))
if __name__=="__main__":main()
