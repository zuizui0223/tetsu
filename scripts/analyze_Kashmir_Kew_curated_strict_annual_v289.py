#!/usr/bin/env python3
"""v28.9 evidence-curated high-impact life span and robust partial identification.

Take 104/123 original railway taxon source matches to Poppenwimer 2023 global
life cycles, then add ONLY original-name-verified perennial species with
independent Kew POWO scientific source evidence. Do not auto-fill any disputed
annual/biennial and perennial/annual species. Main effect is fully observed
27-site average under globally consistent remaining unknown life spans;
bootstrap spans are technical original-site resampling, not population PI.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyze_Kashmir_annual_microsites_v285 import (
    original_taxon_headers,life_span_crosswalk
)
from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix
from analyze_Kashmir_global_unknown_trait_partial_id_v288 import (
    analyze,unknown_taxon_priority,common_source_sites,one_site_linear_contrast,
    partial_identification,WEIGHTS,
)

ALLOWED_TRAITS={"annual","biennial","perennial"}
REQUIRED_OVERRIDES={
    "Cynodon dactylon (L.) Pers.":"perennial",
    "Breea arvensis (L.) Less.":"perennial",
}
EXPLICITLY_WITHHELD={
    "Lactuca serriola L.":"original Kew flora life span annual OR biennial",
    "Senecio chrysanthemoides DC.":"source conflict: Kashmir weed flora annual, official IIIM herbarium perennial",
    "Lotus corniculatus L.":"Kew broad perennial OR annual; subtype uncertainty",
}


def curated_merge(original:list[str],global_crosswalk:Path,independent_overrides:Path):
    existing,meta=life_span_crosswalk(global_crosswalk,original)
    if not {"original_taxon","life_span","trait_source"}.issubset(set(meta.keys())):
        # meta is summary, not actual CSV column list; check through source rows below.
        pass
    expected_names={x for x in original}
    source_rows={}
    with global_crosswalk.open(newline="",encoding="utf-8-sig") as f:
        r=csv.DictReader(f)
        for row in r:
            key=row["original_taxon"]
            if key in source_rows: raise ValueError("Duplicate global source original taxon")
            source_rows[key]=dict(row)
    if len(source_rows)!=len(existing):
        raise ValueError("Original global source crosswalk metadata count mismatch")
    details=[]
    with independent_overrides.open(newline="",encoding="utf-8-sig") as f:
        reader=csv.DictReader(f)
        needed={"original_taxon","life_span","accepted_name",
                "synonym_evidence_url","life_span_evidence_url","trait_source"}
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError("Curated source proof columns missing")
        for row in reader:
            name=(row["original_taxon"] or "").strip()
            span=(row["life_span"] or "").strip().casefold()
            if name not in expected_names:
                raise ValueError("Curated taxon not one of original 123 source names")
            if name in source_rows:
                raise ValueError("Never overwrite accepted verified global-source life-span annotation")
            if name in EXPLICITLY_WITHHELD:
                raise ValueError("Disputed taxon was incorrectly inferred in curated source")
            if REQUIRED_OVERRIDES.get(name)!=span:
                raise ValueError("Curated evidence must match pre-recorded exact source names and values")
            if (not row["accepted_name"] or not row["trait_source"] or
                any(not (row[k] or "").startswith("https://powo.science.kew.org/taxon/")
                    for k in ("synonym_evidence_url","life_span_evidence_url"))):
                raise ValueError("Original Kew scientific accepted-name/life-cycle URL incomplete")
            source_rows[name]={"original_taxon":name,"life_span":span,
                "trait_source":row["trait_source"]+" | "+
                    row["synonym_evidence_url"]+" | "+row["life_span_evidence_url"]}
            details.append({"original_taxon":name,"life_span":span,
                            "accepted_name":row["accepted_name"],
                            "synonym_evidence_url":row["synonym_evidence_url"],
                            "life_span_evidence_url":row["life_span_evidence_url"]})
    if set(x["original_taxon"] for x in details)!=set(REQUIRED_OVERRIDES):
        raise ValueError("Exactly two independently verified source taxa must be present")
    return [source_rows[x] for x in original if x in source_rows],{
        "global_original_taxa":len(existing),
        "new_independent_accepted_name_taxa":len(details),
        "total_curated_original_taxa":len(source_rows),
        "source_original_overrides":details,
        "deliberately_unresolved_top_ranking_examples":EXPLICITLY_WITHHELD,
        "independent_curation_file_sha256":hashlib.sha256(independent_overrides.read_bytes()).hexdigest()
    }


def technical_source_bootstrap(plots:dict,classified:dict[int,str],n_taxa:int,
                               original_sites:list[str],zone:str,reference:str,
                               replicate:int=1999,seed:int=289)->dict:
    if replicate<0 or not original_sites: raise ValueError("Invalid bootstrap")
    site_contrasts=[one_site_linear_contrast(
       plots,classified,n_taxa,s,zone,reference,WEIGHTS["two_year_mean"])
       for s in original_sites]
    observed=partial_identification(site_contrasts)
    if not replicate:
        return {"observed":observed,"bootstrap_draws":0}
    rng=random.Random(seed)
    lowers,uppers=[],[]
    for _ in range(replicate):
        drawn=rng.choices(site_contrasts,k=len(site_contrasts))
        inferred=partial_identification(drawn)
        lowers.append(inferred["sharp_lower"])
        uppers.append(inferred["sharp_upper"])
    lowers.sort();uppers.sort()
    def quantile(v,p):
        t=(len(v)-1)*p
        low=int(t); high=min(low+1,len(v)-1)
        return v[low]+(t-low)*(v[high]-v[low])
    return {"observed":observed,
            "bootstrap_draws":replicate,
            "site_resampling_not_field_randomization_CI":True,
            "lower_bound_technical_q025_q975":[quantile(lowers,.025),quantile(lowers,.975)],
            "upper_bound_technical_q025_q975":[quantile(uppers,.025),quantile(uppers,.975)],
            "fraction_bootstrap_sharp_upper_negative":sum(x<0 for x in uppers)/replicate}


def main()->None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rail",type=Path,required=True)
    ap.add_argument("--global-traits",type=Path,required=True)
    ap.add_argument("--curated-traits",type=Path,required=True)
    ap.add_argument("--merged-traits",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--bootstrap",type=int,default=1999)
    args=ap.parse_args()
    original=original_taxon_headers(args.rail)
    plots,source_meta=parse_xlsx_matrix(args.rail)
    merged,merged_meta=curated_merge(original,args.global_traits,args.curated_traits)
    args.merged_traits.parent.mkdir(parents=True,exist_ok=True)
    with args.merged_traits.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=["original_taxon","life_span","trait_source"])
        w.writeheader();w.writerows(merged)
    traits,trait_meta=life_span_crosswalk(args.merged_traits,original)
    if len(traits)!=len(merged):
        raise AssertionError("Merged source crosswalk did not roundtrip original names")
    report=analyze(plots,traits,len(original))
    sites=common_source_sites(plots)
    result={
        "version":"v28.9",
        "source_original_rail_sha256":source_meta["original_file_sha256"],
        "independent_trait_global_and_kew_audited":True,
        "source_curated_meta":merged_meta,
        "merged_crosswalk_sha256":trait_meta["crosswalk_sha256"],
        "sample_partial_identification_not_global_population_causal_effect":True,
        "original_site_taxon_matrix_not_republished":True,
        "partial_id":report,
        "original_site_bootstrap":technical_source_bootstrap(
            plots,traits,len(original),sites,"Track","Interior Plot 1",
            args.bootstrap),
        "remaining_priority":unknown_taxon_priority(plots,traits,original),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "version":result["version"],
        "matched_species":len(traits),
        "source_sites":len(sites),
        "curated_life_spans":merged_meta["source_original_overrides"],
        "reference_Track_minus_Interior1_all27":
           report["comparisons"]["Track__minus__Interior Plot 1"]["all"],
        "site_bootstrap":result["original_site_bootstrap"],
        "remaining_unknown_width_priority":result["remaining_priority"],
        "output":str(args.output)
    },ensure_ascii=False))


if __name__=="__main__":
    main()
