#!/usr/bin/env python3
"""Independent vegetation-plot species unions on original 200km grid.
No plot IDs, species lists, coordinates or contributor names exported."""
import argparse,csv,datetime,hashlib,io,json,math,re,statistics,subprocess,zipfile
from collections import Counter,defaultdict
from pathlib import Path,PurePosixPath
from pyproj import Transformer
from scan_world_railstatus_verified_plant_size import parse_scaled
from scan_world_railstatus_disturbance_traits import norm,fetch_fig,csv_table,number,LIFE_SHA,SIZE_SHA
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import source_coarse_cell
URL="https://idata.idiv.de/ddm/Data/DownloadZip/3474?version=5047"
DOI="10.25829/idiv.3474-40-3292"
W=200000
PANELS=("1990_2015","2000_2015","1990_2015_explicit_complete")
MISSING={"","na","n/a","nan","null"}
def file_sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""):h.update(b)
    return h.hexdigest()
def year_of(v):
    v=v.strip()
    for fmt in ("%d-%m-%Y","%Y-%m-%d"):
        try:return datetime.datetime.strptime(v,fmt).year
        except ValueError:pass
    return None
def member(z,stem):
    matches=[p for p in z.namelist() if re.fullmatch(r"sPlotOpen_"+stem+r"(?:\s*\(\d+\))?\.txt",PurePosixPath(p).name,re.I)]
    if len(matches)!=1:raise ValueError("Missing/ambiguous sPlotOpen "+stem+" member: "+str(matches))
    return matches[0]
def reader(z,p):
    return io.TextIOWrapper(z.open(p),encoding="utf-8-sig",newline="")
def table_header(r,required):
    if not r.fieldnames or not set(required)<=set(r.fieldnames):raise ValueError("Missing columns: "+str(sorted(set(required)-set(r.fieldnames or []))))
def unique_rows(rows,col):
    d=defaultdict(list)
    for r in rows:d[norm(r[col])].append(r)
    return {k:r[0] for k,r in d.items() if k and len(r)==1}
def make_traits(life_rows,size_rows):
    life=unique_rows(life_rows,"Accepted Name");sz=unique_rows(size_rows,"species")
    known={};height={};qa=Counter()
    for k,r in life.items():
        w=norm(r["Woodiness"]);a=norm(r["LifeCycle"])
        if a=="annual" and w=="woody":qa["contradictory"]+=1;continue
        if w not in ("non-woody","woody"):continue
        known[k]=w;qa[w]+=1
        if w!="non-woody" or k not in sz:continue
        v=number(sz[k]["ph"])
        if v is not None:height[k]=v
    qa["height_nonwoody_exact_names"]=len(height)
    if len(height)<3000 or len(known)<100000:raise ValueError("Insufficient exact source trait pool")
    return known,height,dict(qa)
def read_traits():
    cols,rows=csv_table(fetch_fig(20764777,37005529,LIFE_SHA))
    if cols!=["Family","Genus","Accepted Name","LifeCycle","Woodiness"] or len(rows)!=257401:raise ValueError("LifeForm source schema changed")
    if {norm(r["Woodiness"]) for r in rows}-MISSING-{"non-woody","woody"}:raise ValueError("Unexpected woodiness codes")
    sz,audit=parse_scaled(fetch_fig(25600650,45643398,SIZE_SHA))
    known,heights,qa=make_traits(rows,sz)
    return known,heights,{"LifeForm_SHA256":LIFE_SHA,"height_SHA256":SIZE_SHA,"size_source":audit,"name_mapping":"unique exact full source names; no iNat taxonomy intersection, no synonyms or imputation","name_counts":qa}
def header(z,p):
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    qa=Counter();seen=set();selected={};categories=Counter();datasets=set()
    with reader(z,p) as f:
        r=csv.DictReader(f,delimiter="\t")
        required=("PlotObservationID","Dataset","Date_of_recording","Latitude","Longitude","Location_uncertainty","Releve_area","Plant_recorded")
        table_header(r,required);cols=r.fieldnames
        for row in r:
            qa["source_plots"]+=1;pid=row["PlotObservationID"].strip()
            if not pid or pid in seen:raise ValueError("Duplicate/missing plot ID")
            seen.add(pid)
            typ=norm(row["Plant_recorded"]);categories[typ or "(missing)"]+=1
            year=year_of(row["Date_of_recording"])
            if year is None or not 1990<=year<=2015:continue
            qa["dated_1990_2015"]+=1
            if typ not in ("all vascular plants","not specified"):continue
            qa["complete_or_unspecified"]+=1
            area=number(row["Releve_area"]);u=number(row["Location_uncertainty"])
            lon=number(row["Longitude"]);lat=number(row["Latitude"])
            d=row["Dataset"].strip()
            if area is None or not .01<=area<=40000 or u is None or not 0<u<=1000 or lon is None or lat is None or not d or norm(d) in MISSING:continue
            k=source_coarse_cell(lon,lat,tf,W)
            if k is None:continue
            datasets.add(d);selected[pid]=(k,year,math.log(area),d,typ=="not specified")
            qa["eligible_plots"]+=1
    if qa["source_plots"]!=95104:raise ValueError("Published 95,104 plot census differs")
    return selected,seen,{"qa":dict(qa),"plant_recorded_categories":dict(categories),"columns":cols,"n_selected_datasets":len(datasets)}
def plot_species(z,p,selected,seen):
    by=defaultdict(set);qa=Counter();all_ids=set()
    with reader(z,p) as f:
        r=csv.DictReader(f,delimiter="\t");table_header(r,("PlotObservationID","Species"));cols=r.fieldnames
        for row in r:
            qa["species_records"]+=1;pid=row["PlotObservationID"].strip()
            if pid not in seen:raise ValueError("Species table plot ID missing from header")
            all_ids.add(pid)
            if pid not in selected:continue
            n=norm(row["Species"])
            if not n or n in MISSING:continue
            # Retain resolved taxa for denominator of known life-form coverage;
            # exact joins never guess an unresolved taxon or synonym.
            by[pid].add(n);qa["eligible_species_records"]+=1
    if not 1900000<=qa["species_records"]<=2000000 or len(all_ids)<95000:raise ValueError("Published species-table census differs")
    qa["plots_with_species_rows"]=len(all_ids);qa["selected_plots_with_names"]=len(by)
    return by,{"qa":dict(qa),"columns":cols}
def summarize(selected,by,life,height):
    cells={p:defaultdict(lambda:{"species":set(),"years":[],"areas":[],"datasets":Counter(),"unspecified":0}) for p in PANELS}
    for pid,ss in by.items():
        k,y,a,d,un=selected[pid]
        ps=[PANELS[0]]+([PANELS[1]] if y>=2000 else [])+([PANELS[2]] if not un else [])
        for p in ps:
            c=cells[p][k];c["species"].update(ss);c["years"].append(y);c["areas"].append(a);c["datasets"][d]+=1;c["unspecified"]+=int(un)
    out={}
    audit={}
    for p,grid in cells.items():
        rr=[];qa=Counter()
        for k,c in sorted(grid.items()):
            ss=c["species"];n=len(c["years"])
            den=sum(life.get(x)=="non-woody" for x in ss)
            h=[height[x] for x in ss if x in height]
            lk=sum(x in life for x in ss)
            qa["cells_with_plots"]+=1
            if n<5:continue
            qa["cells_at_least_5_plots"]+=1
            if len(h)<20 or den==0 or len(h)/den<.5:continue
            if len(h)>den or den>lk or lk>len(ss):raise ValueError("Trait denominator inconsistency")
            qa["cells_pass_trait_gate"]+=1
            rr.append({"key":k,"n_plots":n,"species":len(ss),"n_datasets":len(c["datasets"]),"dominant_dataset_share":max(c["datasets"].values())/n,
                       "median_log_area":statistics.median(c["areas"]),"mean_year":statistics.fmean(c["years"]),
                       "unspecified_fraction":c["unspecified"]/n,"known_nonwoody":den,"known_height":len(h),
                       "life_coverage":lk/len(ss),"coverage":len(h)/den,"trait_value":statistics.fmean(h)})
        out[p]=rr;audit[p]=dict(qa)
    return out,audit
def main():
    p=argparse.ArgumentParser()
    for n in ("archive","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();a.archive.parent.mkdir(parents=True,exist_ok=True)
    subprocess.run(["curl","-LfsS","--retry","2","--retry-all-errors","--connect-timeout","30","--max-time","1800","--max-filesize","2000000000",URL,"-o",str(a.archive)],check=True)
    receipt={"url":URL,"doi":DOI,"bytes":a.archive.stat().st_size,"sha256":file_sha(a.archive)}
    print(json.dumps({"splot_source_receipt":receipt}),flush=True)
    if not zipfile.is_zipfile(a.archive):raise ValueError("Official source did not return a ZIP")
    life,height,trait_audit=read_traits()
    with zipfile.ZipFile(a.archive) as z:
        hp,dt=member(z,"header"),member(z,"DT")
        selected,seen,ha=header(z,hp);print(json.dumps({"splot_header_audit":ha}),flush=True)
        by,da=plot_species(z,dt,selected,seen)
    grids,qa=summarize(selected,by,life,height)
    o={"version":"splotopen-independent-stature-grid-v1","source":receipt,"traits":trait_audit,"header_audit":ha,"species_audit":da,"panel_support":qa,"width":W,"grids":grids}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(o,allow_nan=False))
    print(json.dumps({k:v for k,v in o.items() if k!="grids"},allow_nan=False),flush=True)
if __name__=="__main__":main()
