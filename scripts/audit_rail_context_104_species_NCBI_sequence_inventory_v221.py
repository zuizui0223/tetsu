#!/usr/bin/env python3
"""v22.1 inventory ONLY: 104 historic broad-text q>0 Czech railway-context taxa.

Uses original-source authenticated 761-taxon v17.7 sign-unfiltered CSV.
Queries NCBI Entrez Assembly for genome assembly, Nucleotide TSA-master for
assembled transcriptome, SRA WGS strategy as *lead* only.
NO inference about rail adaptation, selection, gene copy number, or
GPS-verified population genomes. Missing hits are NOT evidence of absence.
"""
from __future__ import annotations
import argparse,csv,io,json,random,time,urllib.parse,urllib.request
from collections import Counter
from pathlib import Path

BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
APP="tetsu_public_plant_genomics_inventory_v221"
def request(db,term,limit=2,retries=5):
    params={"db":db,"term":term,"retmode":"json","retmax":str(limit),
            "tool":APP,"email":"none"}
    # No API key: obey NCBI 3 req/s rate guideline, with generous allowance.
    url=BASE+"?"+urllib.parse.urlencode(params)
    why=None
    for i in range(retries):
        time.sleep(0.45 if i==0 else 1.8*(i+1))
        try:
            req=urllib.request.Request(url,headers={"User-Agent":APP+"/0.1 (academic metadata audit)","Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=25) as f:
                j=json.load(f)
            result=j["esearchresult"]
            if result.get("errorlist") or result.get("error"):
                raise ValueError("NCBI reported query error: "+str(result.get("errorlist") or result.get("error")))
            return {"count":int(result["count"]),"ids":result.get("idlist",[]),
                    "translatedquery":result.get("querytranslation",""),"status":"verified"}
        except Exception as ex:
            why=type(ex).__name__+": "+str(ex)[:170]
    return {"count":None,"ids":[],"translatedquery":"","status":"query_error","error":why}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-v177",required=True)
    p.add_argument("--output-csv",required=True)
    p.add_argument("--output-json",required=True)
    p.add_argument("--max",type=int,default=104)
    a=p.parse_args()
    with open(a.source_v177,encoding="utf-8",newline="") as f:
        allrows=list(csv.DictReader(f))
    if not allrows or not {"species","q","q_sign"}.issubset(allrows[0]):
        raise ValueError("Wrong source original 17.7 signed 761-flora taxon schema")
    pos=sorted([x for x in allrows if x["q_sign"]=="positive"],key=lambda x:(-float(x["q"]),x["species"]))
    if len(pos)!=104:
        raise ValueError("Original 104 q-positive broad-text species panel changed")
    assert all(float(x["q"])>0 and int(x["actual_Preslia_station_occupancy_n39"])>=5 for x in pos)
    results=[]
    for idx,r in enumerate(pos[:a.max],1):
        n=r["species"]
        if len(n.split())<2:raise ValueError("Not a valid original species label")
        term='"'+n+'"[Organism]'
        genome=request("assembly",term+" AND latest[filter]",1)
        tsa=request("nuccore",'tsa-master[prop] AND '+term,1)
        wgs=request("sra",term+" AND WGS[Strategy]",1)
        z={
          "rank":idx,"original_taxon":n,
          "original_broad_RAIL_minus_ROAD_q":round(float(r["q"]),10),
          "original_stations_present_n39":int(r["actual_Preslia_station_occupancy_n39"]),
          "original_publication_studies_with_taxon":int(r["original_CzechVeg_source_studies_with_species"]),
          "NCBI_assembly_hit_status":genome["status"],
          "NCBI_assembly_records_n":genome["count"],
          "NCBI_assembly_first_UID":genome["ids"][0] if genome["ids"] else "",
          "NCBI_assembled_genome_detected":("yes" if genome["count"] is not None and genome["count"]>0
                                        else "no_hit_in_this_search" if genome["status"]=="verified"
                                        else "unverified_query_error"),
          "NCBI_TSA_master_hit_status":tsa["status"],
          "NCBI_TSA_master_projects_n":tsa["count"],
          "NCBI_assembled_transcriptome_detected":("yes" if tsa["count"] is not None and tsa["count"]>0
                                        else "no_hit_in_this_search" if tsa["status"]=="verified"
                                        else "unverified_query_error"),
          "NCBI_WGS_SRA_hit_status":wgs["status"],
          "NCBI_WGS_SRA_experiment_records_n":wgs["count"],
          "WGS_SRA_is_NOT_georeferenced_population_sample_n":True,
          "GPS_verified_population_WGS_n":None,
          "GPS_population_genomics_status":"not_verified",
          "original_czech_railway_context_not_adaptation":True
        }
        for tag,x in (("assembly",genome),("TSA",tsa),("WGS",wgs)):
            if x["status"]!="verified":z["query_error_"+tag]=x.get("error")
        results.append(z)
        print(f"screened_{idx}_of_{a.max} species={n} assembly={z['NCBI_assembled_genome_detected']} tsa={z['NCBI_assembled_transcriptome_detected']} wgs_records={wgs['count']}",flush=True)
    out=Path(a.output_csv);out.parent.mkdir(parents=True,exist_ok=True)
    fields=list(results[0]); fields=[x for x in fields if not x.startswith("query_error_")]+["query_error_assembly","query_error_TSA","query_error_WGS"]
    with out.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore");w.writeheader()
        w.writerows(results)
    g=Counter(x["NCBI_assembled_genome_detected"] for x in results)
    t=Counter(x["NCBI_assembled_transcriptome_detected"] for x in results)
    w=Counter("hit" if x["NCBI_WGS_SRA_experiment_records_n"] and x["NCBI_WGS_SRA_experiment_records_n"]>0 else
              "zero" if x["NCBI_WGS_SRA_experiment_records_n"]==0 else "error" for x in results)
    j={"status":"SOURCE_INPUT_VERIFIED_NCBI_METADATA_SCREEN",
       "n_source_original_rail_context_qpositive_eligible":len(pos),
       "n_screened":len(results),
       "genome_assembly_status":dict(g),
       "assembled_TSA_status":dict(t),"NCBI_SRA_WGS_lead_status":dict(w),
       "n_GPS_verified_population_genomes":None,
       "pop_genomes_status":"NOT_AUDITED: do not call WGS experiment counts individual samples",
       "meta_source":"tetsu genuine v17.7 sign-unfiltered 39-station CzechVeg original-source report",
       "NCBI_search_scope":"NCBI Assembly latest + nuccore TSA master + SRA WGS strategy, exact source binomial",
       "limitations":["No explicit assembly accession was resolved: returned first Entrez UID only",
          "A positive WGS record does not imply a georeferenced individual or railway contrast",
          "No NCBI hit is not genome/transcriptome absence across all repositories",
          "BROAD rail-context is a text field, not a surveyed active railway track",
          "Genome/transcriptome availability is not evidence of selection"]}
    Path(a.output_json).write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding="utf-8")
    print("NCBI_INVENTORY_SUMMARY_JSON="+json.dumps(j,ensure_ascii=False,separators=(",",":")),flush=True)
    print("NCBI_INVENTORY_ROWS_JSON="+json.dumps(results,ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()
