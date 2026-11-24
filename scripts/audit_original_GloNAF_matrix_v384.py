#!/usr/bin/env python3
"""v384: audit GloNAF original nonoverlap matrix vs exact direct regional lists."""
import argparse,csv,hashlib,json,statistics
from collections import Counter,defaultdict
from pathlib import Path
from audit_global_GloNAF_rail2048_v378a import verified,load,unique_map,norm,valid_name,audit_documented_rows,read_original_nol
from analyze_global_rail_documented_naturalized_v379 import source_lists

TX_MD5="af2b5ea09f3fe6c677c7a2ceadd2b543"
DECISION_MD5="1e2547cdb3150aa571cdb006b5f25f68"
DICT_MD5="2e6aedb340ec66e14e4d6c85c5ca9bb8"
def md5(path,expect):
    with path.open("rb") as f:actual=hashlib.file_digest(f,"md5").hexdigest()
    if actual!=expect:raise ValueError("Original publisher MD5 mismatch: "+str(path))
def num_summary(values):
    v=sorted(values)
    if not v:return None
    return {"n":len(v),"min":v[0],"median":statistics.median(v),"max":v[-1],"mean":statistics.mean(v)}
def compare(a,b):
    common=len(a&b);union=len(a|b)
    return {"direct":len(a),"matrix":len(b),"shared":common,
        "direct_only":len(a-b),"matrix_only":len(b-a),"Jaccard":common/union if union else None}
def run(root):
    verified(root)
    for filename,expect in (("glonaf_TxR.csv",TX_MD5),("glonaf_decision_nol_region.csv",DECISION_MD5),("glonaf_flora2_datadictionary.csv",DICT_MD5)):
        md5(root/filename,expect)
    print((root/"glonaf_flora2_datadictionary.csv").read_text(encoding="utf-8-sig"),flush=True)
    g,direct,qa=source_lists(root)
    regions=unique_map(load(root/"glonaf_region.csv",{"id","code","name","OBJIDsic"}),"id")
    bycode={v["code"]:k for k,v in regions.items()}
    taxa=load(root/"glonaf_taxon_wcvp.csv",{"accepted_plant_name_id","taxa_accepted","wcvp_taxon_rank"})
    byname=defaultdict(set);byid=defaultdict(set)
    for t in taxa:
        n=norm(t["taxa_accepted"]);aid=norm(t["accepted_plant_name_id"])
        if valid_name(n) and norm(t["wcvp_taxon_rank"]).lower()=="species" and aid.lower() not in ("","na","nan"):
            byname[n].add(aid);byid[aid].add(n)
    strict={aid:next(iter(n)) for aid,n in byid.items() if len(n)==1 and len(byname[next(iter(n))])==1}
    flora=load(root/"glonaf_flora2.csv",{"status","region_id","list_id","taxon_wcvp_id"})
    status=Counter(r["status"] for r in flora)
    tx=root/"glonaf_TxR.csv"
    matrix=defaultdict(set);states=Counter();rows=0;exact_rows=0;freq_errors=0
    ids=set();distinct=set()
    with tx.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f);fields=reader.fieldnames
        keys=[k for k in fields if k in bycode]
        if len(keys)!=984 or set(fields)-set(keys)!={"wcvp12_accepted_plant_name_id","taxa_accepted","frequency"}:
            raise ValueError("Unexpected official TxR schema")
        for r in reader:
            rows+=1
            if None in r or any(v is None for v in r.values()):raise ValueError("Malformed TxR row")
            aid=norm(r["wcvp12_accepted_plant_name_id"]);name=norm(r["taxa_accepted"])
            if not aid or aid in ids:raise ValueError("Duplicate or absent original accepted ID")
            ids.add(aid)
            vals=[norm(r[k]) for k in keys]
            states.update(vals)
            # Official incidence is encoded 0/1. Do not treat unspecified
            # values, NA or arbitrary nonzero counts as presences.
            if not set(vals)<={"0","1","0.0","1.0"}:
                raise ValueError("Unrecognized original matrix incidence values: "+str(Counter(vals)))
            present=[k for k,v in zip(keys,vals) if v in ("1","1.0")]
            if float(r["frequency"])!=len(present):freq_errors+=1
            if strict.get(aid)!=name:continue
            exact_rows+=1;distinct.add(name)
            for k in present:matrix[bycode[k]].add(name)
    if rows!=15892 or freq_errors:
        raise ValueError("Original matrix row count or occurrence-frequency audit failed")
    geom_ids=set(g["region_id"])-{None}
    results=[];missing=[]
    for code in keys:
        rid=bycode[code]
        a=direct[rid];b=matrix[rid]
        d=compare(a,b)
        d.update(region_code=code,has_unambiguous_geometry=rid in geom_ids)
        results.append(d)
        if rid not in geom_ids:missing.append({"code":code,"name":regions[rid]["name"]})
    selected=[r for r in results if r["has_unambiguous_geometry"]]
    # Original two suspect polygon rows: expose only public region metadata.
    raw=read_original_nol(root)
    byobj={norm(v["OBJIDsic"]):k for k,v in regions.items() if norm(v["OBJIDsic"]) and norm(v["OBJIDsic"]).lower()!="na"}
    rawids=[byobj.get(norm(v)) for v in raw["OBJIDsic"]];counts=Counter(rawids)
    suspect=[{k:str(row[k]) for k in ("OBJIDsic","IDregion","ISOcountry","name")}
             for i,(_,row) in enumerate(raw.iterrows()) if counts[rawids[i]]>1 or rawids[i] is None]
    # The publisher decision table is logged at schema level; no guessed
    # hierarchy aggregation or manual reassignment of a disputed polygon.
    with (root/"glonaf_decision_nol_region.csv").open(encoding="utf-8-sig",newline="") as f:
        d=csv.DictReader(f);decisions=list(d);decision_fields=d.fieldnames
    out={"version":"v384-original-GloNAF-direct-vs-official-nonoverlap-incidence",
         "source_original_TXR_MD5":TX_MD5,"source_decision_table_MD5":DECISION_MD5,
         "source_flora_dictionary_MD5":DICT_MD5,"source_flora_status":dict(status),
         "source_matrix_total_rows":rows,"source_matrix_exact_binomial_rows":exact_rows,
         "source_matrix_exact_species_n":len(distinct),"source_matrix_values":dict(states),
         "source_matrix_frequency_mismatches":freq_errors,
         "unambiguous_geometry_regions":len(selected),"matrix_regions_without_safe_geometry":missing,
         "original_suspect_polygon_metadata":suspect,
         "decision_table_original_fields":decision_fields,"decision_table_rows":len(decisions),
         "n_equal_exact_lists_all":sum(r["direct_only"]==0 and r["matrix_only"]==0 for r in results),
         "n_equal_exact_lists_unambiguous":sum(r["direct_only"]==0 and r["matrix_only"]==0 for r in selected),
         "n_unambiguous_with_matrix_extra":sum(r["matrix_only"]>0 for r in selected),
         "n_unambiguous_with_direct_extra":sum(r["direct_only"]>0 for r in selected),
         "sum_direct_memberships_unambiguous":sum(r["direct"] for r in selected),
         "sum_matrix_memberships_unambiguous":sum(r["matrix"] for r in selected),
         "sum_matrix_only_memberships_unambiguous":sum(r["matrix_only"] for r in selected),
         "sum_direct_only_memberships_unambiguous":sum(r["direct_only"] for r in selected),
         "unambiguous_region_Jaccard":num_summary([r["Jaccard"] for r in selected if r["Jaccard"] is not None]),
         "region_comparisons":results,
         "interpretation":"source classification comparison ONLY; no new iNat fraction or railway effect estimated; unresolved polygons remain excluded"}
    print(json.dumps({k:v for k,v in out.items() if k!="region_comparisons"},allow_nan=False),flush=True)
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    r=run(a.source);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,allow_nan=False))
if __name__=="__main__":main()
