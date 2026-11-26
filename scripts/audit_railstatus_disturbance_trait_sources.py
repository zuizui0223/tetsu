#!/usr/bin/env python3
"""Audit authentic public trait tables before the rail-status trait analysis.

No occurrence scan, no rail outcomes, no guessed category semantics.
Exact active iNaturalist species-name support and trait-column coverage only.
"""
import argparse,csv,io,json,hashlib,math,re,zipfile
from collections import Counter,defaultdict
from pathlib import Path
from audit_world_NHM_rail_trait_match import norm,request
from scan_world_railstatus_cross_taxa import taxonomy,source_client,open_original_source
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import TAXA_ETAG,TAXA_BYTES
from scan_world_railstatus_bird_avonet_guilds import DERIVATIVE
from scan_world_railstatus_bird_conditional_hwi import SOURCE_SHA

LIFE_ARTICLE=20764777
LIFE_FILE=37005529
SIZE_ARTICLE=25600650
LEP_COMMIT="f73317d28b38f80923ad8e02cc06d35f5c16ed1f"
LEP_PATH="consensus/consensus.csv"
LEP_BLOB="b2bc7bd87eca199b62f08ac539660aebc9fbaf79"
LEP_ROOT="https://raw.githubusercontent.com/RiesLabGU/LepTraits/"+LEP_COMMIT+"/"
MISSING={"","na","n/a","nan","null","none","-","unknown"}
NAME_RX=re.compile(r"^[A-Za-z][A-Za-z0-9. -]* [A-Za-z][A-Za-z0-9. -]*$")

def git_blob_sha(raw):
    return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()

def figshare_metadata(article):
    q=json.loads(request(f"https://api.figshare.com/v2/articles/{article}",cap=2_000_000))
    if q.get("id")!=article or not q.get("files"):raise ValueError("Figshare article identity/files missing")
    return q

def get_figshare(f,cap):
    size=f["size"]
    if not isinstance(size,int) or not 0<size<=cap:raise ValueError("File outside frozen size cap")
    url=f["downloadUrl"]
    if not url.startswith("https://"):raise ValueError("Unexpected file URL")
    raw=request(url,cap=cap)
    md5=f.get("computed_md5") or f.get("supplied_md5")
    if len(raw)!=size or not md5 or hashlib.md5(raw).hexdigest()!=md5.lower():raise ValueError("Figshare size or MD5 mismatch")
    return raw

def receipt(raw,**kw):
    return dict(kw,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())

def csv_table(raw):
    text=raw.decode("utf-8-sig")
    sample=text[:100000]
    try:dialect=csv.Sniffer().sniff(sample,delimiters=",\t;")
    except csv.Error:dialect=csv.excel
    r=csv.DictReader(io.StringIO(text,newline=""),dialect=dialect)
    columns=r.fieldnames or []
    if len(columns)<2 or len(columns)!=len(set(columns)):raise ValueError("Malformed/duplicate table columns")
    rows=list(r)
    if not rows or any(None in r or any(v is None for v in r.values()) for r in rows):raise ValueError("Incomplete rows")
    return columns,rows

def species_names(client,canonical,sha):
    names={g:defaultdict(list) for g in ("plants","insects","birds")}
    rr,raw,stream=open_original_source(client,"taxa.csv.gz",TAXA_ETAG,TAXA_BYTES)
    try:
        for r in rr:
            t=int(r["taxon_id"]);hit=canonical.get(t)
            if hit is not None and hit[0] in names and hit[1]==t and r["rank"].strip().lower()=="species":
                names[hit[0]][norm(r["name"])].append(t)
    finally:stream.close()
    if raw.bytes!=TAXA_BYTES or raw.hash.hexdigest()!=sha:raise ValueError("Exact original taxonomy failed")
    return names

def label_stats(rows,c):
    counts=Counter(norm(r[c]) for r in rows)
    valid={k:v for k,v in counts.items() if k not in MISSING}
    num=[]
    for k,n in valid.items():
        try:x=float(k)
        except (TypeError,ValueError):continue
        if math.isfinite(x):num.append((x,n))
    # Do not export species names or long free text; only low-cardinality
    # category codes and numerical summaries for codebook review.
    cat={k:v for k,v in sorted(valid.items())} if len(valid)<=30 and all(len(k)<=80 for k in valid) else None
    return {"nonmissing":sum(valid.values()),"distinct_nonmissing":len(valid),
            "numeric_rows":sum(n for x,n in num),"numeric_range":[min(x for x,n in num),max(x for x,n in num)] if num else None,
            "categories_if_low_cardinality":cat}

def audit_table(raw,names,source):
    columns,rows=csv_table(raw)
    # Identify possible name columns by exact matches, not by guessing
    # what a trait column means. No synonyms, no genus-level assignment.
    options=[]
    for c in columns:
        if not re.search(r"species|taxon|name|binomial|accepted|scientific|^sp$|^taxa$",c,re.I):continue
        by=Counter(norm(r[c]) for r in rows if norm(r[c]) not in MISSING)
        matched={k:names[k][0] for k,n in by.items() if n==1 and len(names.get(k,[]))==1}
        options.append({"column":c,"unique_exact_species":len(set(matched.values())),"duplicate_source_names":sum(n>1 for n in by.values()),
                        "matched_rows_with_nonmissing_by_column":{k:sum(norm(r[c]) in matched and norm(r[k]) not in MISSING for r in rows) for k in columns}})
    return {"source":source,"n_rows":len(rows),"columns":columns,"name_candidates":options,"column_audit":{c:label_stats(rows,c) for c in columns}}

def inspect_size_source():
    # Official paper data availability identifies this article. Check
    # metadata only initially; avoid arbitrary archive extraction.
    q=figshare_metadata(SIZE_ARTICLE)
    files=[{k:f.get(k) for k in ("id","name","size","computed_md5","supplied_md5","downloadUrl")} for f in q["files"]]
    return {"article_id":SIZE_ARTICLE,"title":q.get("title"),"doi":q.get("doi"),"version":q.get("version"),"files":files,
            "interpretation":"Metadata only. Inspect documented observed trait table and units; do not use imputed values."}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--taxa",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    frozen=json.loads(a.taxa.read_text());client=source_client()
    canonical,_,sha=taxonomy(client)
    if sha!=frozen["taxa_sha256"]:raise ValueError("Frozen iNat taxonomy changed")
    names=species_names(client,canonical,sha)
    result={"version":"railstatus-disturbance-trait-source-audit-v1","taxonomy_sha256":sha,"sources":{}}
    def run(label,fn):
        try:v=fn();z={"status":"retrieved_and_audited",**v}
        except Exception as e:z={"status":"blocked","type":type(e).__name__,"message":str(e)[:300]}
        result["sources"][label]=z
        print(json.dumps({"trait_source":label,**z},allow_nan=False),flush=True)
    def life():
        q=figshare_metadata(LIFE_ARTICLE)
        fs=[f for f in q["files"] if f["id"]==LIFE_FILE and f["name"].lower().endswith(".csv")]
        if len(fs)!=1:raise ValueError("Published life-form file missing")
        raw=get_figshare(fs[0],30_000_000)
        z=audit_table(raw,names["plants"],receipt(raw,article_id=LIFE_ARTICLE,file_id=LIFE_FILE,doi=q.get("doi"),name=fs[0]["name"],version=q.get("version")))
        if z["n_rows"]<200000:raise ValueError("Incomplete original life-form table")
        return z
    def lep():
        raw=request(LEP_ROOT+LEP_PATH,cap=8_000_000)
        if len(raw)!=2687158 or git_blob_sha(raw)!=LEP_BLOB:raise ValueError("Pinned LepTraits blob mismatch")
        z=audit_table(raw,names["insects"],receipt(raw,commit=LEP_COMMIT,path=LEP_PATH,blob=LEP_BLOB,doi="10.1038/s41597-022-01473-5"))
        if z["n_rows"]<10000:raise ValueError("Incomplete LepTraits table")
        key=request(LEP_ROOT+"records/habitat_recordKey.csv",cap=100000)
        if git_blob_sha(key)!="6dad5ae6f047bc36eab47671d9fc49e557cf9b19":raise ValueError("Pinned habitat codebook mismatch")
        cols,rr=csv_table(key)
        z["habitat_codebook"]={"sha256":hashlib.sha256(key).hexdigest(),"columns":cols,"rows":rr}
        return z
    def bird():
        raw=request(DERIVATIVE,cap=10_000_000)
        if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:raise ValueError("Pinned AVONET derivative changed")
        return audit_table(raw,names["birds"],receipt(raw,url=DERIVATIVE,doi="10.1111/ele.13898"))
    run("plant_life_cycle",life);run("butterfly_traits",lep);run("bird_traits",bird);run("plant_size_publication_files",inspect_size_source)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False),encoding="utf-8")
    if all(v["status"]=="blocked" for v in result["sources"].values()):raise ValueError("All source routes blocked")
if __name__=="__main__":main()
