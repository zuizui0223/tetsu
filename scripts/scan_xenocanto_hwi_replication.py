#!/usr/bin/env python3
"""Independent recording-platform support and pair aggregates from source DwC-A.
No audio, occurrence IDs, names of recordists or fine locations exported."""
import argparse,csv,hashlib,io,json,math,re,subprocess,unicodedata,zipfile
from pathlib import Path,PurePosixPath
from collections import Counter,defaultdict
from xml.etree import ElementTree as ET
from pyproj import Transformer
from scan_global_grip4_road_crossings import candidates,key_sha
from scan_world_railstatus_bird_conditional_hwi import (validate,csv_rows,DERIVATIVE,request,names,conditioned_match,W,pair_key,taxonomy,source_client,source_coarse_cell,bucket,add,pair_rows,SOURCE_SHA,AXES,HWI,GUILDS,group_pair,MIN_GUILD_HALF,norm)
URL="https://api.biodiversitydata.nl/v2/specimen/dwca/getDataSet/wildlife-sounds-birds"
DOI="10.15468/qv0ksn"
PANELS=("2018_2025_broad","2018_2025_known_accuracy","2018_2022_broad")
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""):h.update(b)
    return h.hexdigest()
def local(tag):return tag.rsplit("}",1)[-1]
def child(parent,name):
    return next((x for x in parent if local(x.tag)==name),None)
def esc(v):
    return {"\\t":"\t","\\n":"\n","\\r":"\r","\\r\\n":"\r\n"}.get(v,v)
def core_schema(z):
    metas=[n for n in z.namelist() if PurePosixPath(n).name=="meta.xml"]
    if len(metas)!=1:raise ValueError("Expected unique DwC-A meta.xml")
    mp=PurePosixPath(metas[0]);root=ET.fromstring(z.read(str(mp)))
    core=child(root,"core")
    if core is None or core.get("rowType","").rsplit("/",1)[-1]!="Occurrence":raise ValueError("Not occurrence core")
    loc=child(child(core,"files"),"location")
    identifier=child(core,"id")
    if loc is None or identifier is None:raise ValueError("Missing source core file or ID")
    file=str(mp.parent/loc.text.strip())
    if file not in z.namelist():raise ValueError("Missing occurrence member")
    fields={}
    for f in core:
        if local(f.tag)!="field":continue
        term=f.get("term","").rsplit("/",1)[-1].rsplit("#",1)[-1]
        if not term or term in fields:raise ValueError("Duplicate or absent field term")
        fields[term]=(int(f.attrib["index"]) if "index" in f.attrib else None,f.get("default",""))
    required={"decimalLatitude","decimalLongitude","recordedBy","scientificName"}
    if not required<=set(fields) or not {"eventDate","year"}&set(fields):raise ValueError("Missing required occurrence fields: "+str(sorted(required-set(fields))))
    sep=esc(core.get("fieldsTerminatedBy","\\t"));quote=esc(core.get("fieldsEnclosedBy",'"'))
    if len(sep)!=1 or len(quote)>1:raise ValueError("Unsupported DwC delimiter")
    return file,fields,int(identifier.attrib["index"]),core.get("encoding","UTF-8"),sep,quote,int(core.get("ignoreHeaderLines","0"))
def field(row,fields,key):
    if key not in fields:return ""
    i,d=fields[key]
    if i is None:return d
    if i>=len(row):raise ValueError("Truncated source row")
    return row[i].strip() or d
def year_of(date,year):
    if date:
        if "/" in date or not re.match(r"^\d{4}(?:$|-\d{2}(?:$|-\d{2}(?:$|T)))",date):return None
        y=int(date[:4])
        if year and (not year.isdigit() or int(year)!=y):return None
        return y
    return int(year) if re.fullmatch(r"\d{4}",year) else None
def binomial(get):
    g,e=get("genus"),get("specificEpithet")
    sn=get("scientificName");rank=get("taxonRank").casefold()
    if rank and rank not in ("species","subspecies","variety","form","infraspecies"):return None
    if g and e:s=g+" "+e
    else:
        parts=sn.split()
        if len(parts)==2:s=sn
        elif rank in ("species","subspecies","infraspecies") and len(parts)>=2:s=" ".join(parts[:2])
        else:return None
    if not re.fullmatch(r"[A-Z][a-z]+ [a-z][a-z-]+",s):return None
    if re.search(r"\b(?:cf|aff|sp|hybrid)\b|×",sn,re.I):return None
    return norm(s)
def observer(v):
    v=unicodedata.normalize("NFKC",v).casefold()
    if any(x in v for x in (";","|")):return None
    v=" ".join(v.split())
    return v if v and v not in ("unknown","anonymous","not provided","n/a","null") else None
def absent(v):return not v or v.casefold() in ("none","no","not applicable","n/a")
def scan_archive(path,keys,traits):
    cells={k for pair in keys for k in pair}
    full={p:defaultdict(bucket) for p in PANELS}
    guild={p:{g:defaultdict(bucket) for g in GUILDS} for p in PANELS}
    tf=Transformer.from_crs(4326,6933,always_xy=True)
    qa=Counter();status=Counter();basis=Counter();matched=set();seen=set()
    with zipfile.ZipFile(path) as z:
        emls=[n for n in z.namelist() if PurePosixPath(n).name.lower()=="eml.xml"]
        if len(emls)!=1:raise ValueError("Missing unique source EML")
        eml=ET.fromstring(z.read(emls[0]))
        titles=[(e.text or "") for e in eml.iter() if local(e.tag)=="title"]
        if not any("xeno-canto" in v.casefold() and "bird" in v.casefold() for v in titles):raise ValueError("Wrong source EML title")
        published=[(e.text or "") for e in eml.iter() if local(e.tag)=="pubDate"]
        member,fields,idcol,enc,sep,quote,headers=core_schema(z)
        with z.open(member) as f,io.TextIOWrapper(f,encoding=enc,newline="") as stream:
            rd=csv.reader(stream,delimiter=sep,quotechar=quote or None,quoting=csv.QUOTE_MINIMAL if quote else csv.QUOTE_NONE)
            for _ in range(headers):next(rd)
            for row in rd:
                qa["rows"]+=1
                if idcol>=len(row) or not row[idcol] or row[idcol] in seen:raise ValueError("Missing/duplicate occurrence core ID")
                seen.add(row[idcol])
                get=lambda k:field(row,fields,k)
                y=year_of(get("eventDate"),get("year"))
                if y is None or not 2018<=y<=2025:continue
                qa["in_period"]+=1
                b=get("basisOfRecord");basis[b or "(missing)"]+=1
                st=get("identificationVerificationStatus").casefold();status[st or "(missing)"]+=1
                if any(t in st for t in ("questioned","uncertain","unverified","rejected")) and not any(t in st for t in ("not questioned","unquestioned","resolved")):
                    qa["excluded_id_status"]+=1;continue
                if get("occurrenceStatus").casefold() not in ("","present"):continue
                if not absent(get("dataGeneralizations")) or not absent(get("informationWithheld")):
                    qa["excluded_generalization_or_withholding"]+=1;continue
                datum=re.sub(r"[^a-z0-9]","",get("geodeticDatum").casefold())
                if datum not in ("","wgs84","epsg4326","wgs1984"):qa["excluded_datum"]+=1;continue
                t=binomial(get);o=observer(get("recordedBy"))
                if t is None or o is None:qa["excluded_name_or_recordist"]+=1;continue
                try:lon=float(get("decimalLongitude"));lat=float(get("decimalLatitude"))
                except ValueError:qa["excluded_coordinates"]+=1;continue
                if not (math.isfinite(lon) and math.isfinite(lat) and -180<=lon<=180 and -90<=lat<=90) or (lon==0 and lat==0):continue
                u=get("coordinateUncertaintyInMeters")
                known=False
                if u:
                    try:v=float(u)
                    except ValueError:qa["excluded_uncertainty"]+=1;continue
                    if not math.isfinite(v) or v<0 or v>1000:qa["excluded_uncertainty"]+=1;continue
                    known=v>0
                qa["eligible_global_records"]+=1
                k=source_coarse_cell(lon,lat,tf,W)
                if k not in cells:continue
                qa["eligible_frozen_records"]+=1
                if not known:qa["frozen_unknown_uncertainty"]+=1
                groups=traits.get(t,())
                if groups:matched.add(t);qa["frozen_trait_matched_records"]+=1
                ps=[PANELS[0]]+([PANELS[1]] if known else [])+([PANELS[2]] if y<=2022 else [])
                for p in ps:
                    add(full[p][k],t,o)
                    for g in groups:add(guild[p][g][k],t,o)
        if qa["rows"]<100000:raise ValueError("Implausibly small full source archive")
    out={};support={}
    for p in PANELS:
        rr=[]
        for r in pair_rows(full[p]):
            k=pair_key(r)
            if k not in keys or "split_similarity" not in r:continue
            a,b=k;gg={}
            for g in GUILDS:
                A=guild[p][g].get(a);B=guild[p][g].get(b)
                if A is not None and B is not None:
                    m=group_pair(A,B)
                    if m is not None:gg[g]=m
            rr.append({"a":a,"b":b,"full":r,"guilds":gg})
        out[p]=rr
        support[p]={"n_candidate_cells_observed":len(full[p]),"n_full_split_pairs":len(rr),"n_HWI_pairs":sum(all(g in r["guilds"] for g in HWI) for r in rr),"n_joint_stratum_pairs":sum(any(all(g in r["guilds"] for g in gs) for gs in AXES.values()) for r in rr)}
    return out,{"source_core_member":member,"source_fields":sorted(fields),"eml_titles":titles,"eml_publication_dates":published,"counts":dict(qa),"basis_counts":dict(basis),"identification_status_counts":dict(status),"n_trait_matched_species_on_frozen_geography":len(matched),"support":support}
def main():
    p=argparse.ArgumentParser()
    for n in ("historical","archive","out"):p.add_argument("--"+n,type=Path,required=True)
    a=p.parse_args();keys=candidates(a.historical);hist=json.loads(a.historical.read_text())
    a.archive.parent.mkdir(parents=True,exist_ok=True)
    subprocess.run(["curl","-LfsS","--retry","2","--retry-all-errors","--connect-timeout","30","--max-time","1800","--max-filesize","1500000000",URL,"-o",str(a.archive)],check=True)
    raw=request(DERIVATIVE,cap=10000000)
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:raise ValueError("Trait SHA mismatch")
    rows=validate(csv_rows(raw));client=source_client();canonical,_,tsha=taxonomy(client)
    if tsha!=hist["taxa_sha256"]:raise ValueError("Original taxonomy SHA mismatch")
    n=names(client,canonical,tsha);trait,qa=conditioned_match(rows,n)
    if qa!=hist["source"]["match"]:raise ValueError("Frozen trait assignment mismatch")
    traits={k:trait[ids[0]] for k,ids in n.items() if len(ids)==1 and ids[0] in trait}
    pairs,audit=scan_archive(a.archive,keys,traits)
    obj={"version":"xc-independent-hwi-pairs-v1","url":URL,"doi":DOI,"source_sha256":sha(a.archive),"source_bytes":a.archive.stat().st_size,"trait_sha256":SOURCE_SHA,"trait_match_replayed":True,"candidate_sha256":key_sha(keys),"n_candidates":len(keys),"guild_half_gate":list(MIN_GUILD_HALF),"audit":audit,"pairs":pairs}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(obj,allow_nan=False));print(json.dumps({k:v for k,v in obj.items() if k!="pairs"},allow_nan=False),flush=True)
if __name__=="__main__":main()
