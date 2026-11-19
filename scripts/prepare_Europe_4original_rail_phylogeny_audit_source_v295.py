#!/usr/bin/env python3
"""v29.5 original four EUROPE RAILWAY floras -> PRIVATE CI-only phylogenetic tip QA.

The site-by-species/trait rows exported by this helper are transient within
GitHub Actions runner and MUST NEVER be committed, uploaded or printed.
Original source hashes, 4 publication cohorts, site counts and v24.5 LE DA
global zlog axes are frozen by v28.5 frozen_reference() before any phylo audit.
The three original V.PhyloMaker2 mega-trees are tested against only
source-authentic binomial species with both positive LEDA traits.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from pathlib import Path
from collections import Counter,defaultdict

from audit_Europe_4independent_railstudies_fixed_community_trait_label_null_v285 import (
    STUDIES,N,LEDASHA,publisher_czech,publisher_italy,publisher_hr_ba,
    extract_original_xlsx, frozen_reference
)
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda

def encode_original_source_rows(cohorts:dict,traits:dict)->tuple[list[dict],dict]:
    if set(cohorts)!=set(STUDIES):raise ValueError("Four true original publication cohorts missing")
    rows=[]
    qa={}
    for study in STUDIES:
        sites=cohorts[study]
        if len(sites)!=N[study]:raise ValueError("True original station or ballast sites drift")
        group_ids={g:str(i+1) for i,g in enumerate(sorted({x['group'] for x in sites}))}
        source_species=set()
        incidence=0
        for j,site in enumerate(sites):
            names=site['taxa']
            if len(names)<5:raise ValueError("Original source site has too few dual-trait species")
            for taxon in sorted(names):
                if taxon not in traits[study]:
                    raise ValueError("Original real source species had no paired LEDA canopy+seed trait")
                height,seed=traits[study][taxon]
                if not (math.isfinite(height) and math.isfinite(seed)):
                    raise ValueError("Original source trait coordinate nonfinite")
                if len(taxon.split())!=2:raise ValueError("Source matched binomial not exact")
                rows.append({
                    "study":study,
                    "site_num":j+1,
                    "group_num":group_ids[site['group']],
                    "taxon":taxon,
                    "zlog_canopy_height":height,
                    "zlog_seed_mass":seed,
                })
                incidence+=1
                source_species.add(taxon)
        qa[study]={
            'true_original_field_sites':len(sites),
            'true_original_geographical_groups':len(group_ids),
            'species_exact_LEDA_dual_trait_source_count':len(source_species),
            'source_field_site_by_taxon_incidence':incidence,
            'n_original_all_sites_with_ge5_dual_trait_species':len(sites),
        }
    if sum(q['true_original_field_sites'] for q in qa.values())!=171:
        raise ValueError('All original studies did not preserve exactly 171 original physical site units')
    if len(rows)!=len({(r["study"],r["site_num"],r["taxon"]) for r in rows}):
        raise ValueError("Repeated source study/site/species observation inserted")
    return rows,qa

def run(cz_s1:Path,cz_s2:Path,it_pdf:Path,hr_pdf:Path,
        si_xlsx:Path,si_audit:Path,height:Path,seed:Path,
        private_out:Path,aggregate_out:Path)->dict:
    for file,key in ((height,"height"),(seed,"seed")):
        if hashlib.sha256(file.read_bytes()).hexdigest()!=LEDASHA[key]:
            raise ValueError("Original LE DA source checksum changed")
    h,hqa=read_leda(height,"height_m")
    s,sqa=read_leda(seed,"seed_mass_mg")
    published={
        "CZ_Czechia":publisher_czech(cz_s1,cz_s2),
        "IT_Italy_Parma":publisher_italy(it_pdf),
        "HR_BA_2017_one_study":publisher_hr_ba(hr_pdf)
    }
    studies={k:v[0] for k,v in published.items()}
    original_qa={k:v[1] for k,v in published.items()}
    original_si=extract_original_xlsx(si_xlsx,si_audit)
    cohorts,trait,qa=frozen_reference(
        (studies,{**original_qa,"slovenia":original_si}),
        {"height":h,"seed":s})
    rows,support=encode_original_source_rows(cohorts,trait)
    private_out.parent.mkdir(parents=True,exist_ok=True)
    with private_out.open("w",encoding="utf-8",newline="") as out:
        writer=csv.DictWriter(out,fieldnames=[
            "study","site_num","group_num","taxon",
            "zlog_canopy_height","zlog_seed_mass"])
        writer.writeheader()
        writer.writerows(rows)
    aggregate={
        "version":"v29.5_original_phylo_source_only",
        "original_four_independent_publisher_plant_flora_units":support,
        "n_original_publications":len(STUDIES),
        "n_original_field_sites_total":171,
        "source_original_LEDA_sha":{
            "height":hqa["sha256_actual_download"],
            "seed_mass":sqa["sha256_actual_download"]},
        "source_three_Europe_publication_zlog_axes_frozen":True,
        "source_taxa_per_site_and_local_traits_are_PRIVATE_CI_only_not_released":True,
        "author_original_sites_and_species_taxa_membership_changed":False,
        "phylogenetic_tip_mapping_has_NOT_been_inferred_yet":True,
    }
    aggregate_out.parent.mkdir(parents=True,exist_ok=True)
    aggregate_out.write_text(json.dumps(aggregate,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "study_source_species_counts":{k:v["species_exact_LEDA_dual_trait_source_count"] for k,v in support.items()},
        "study_source_sites":{k:v["true_original_field_sites"] for k,v in support.items()},
        "original_true_site_taxon_incidence_total":len(rows),
        "private_source_taxa_not_printed":True},ensure_ascii=False))
    return aggregate

def main():
    p=argparse.ArgumentParser()
    for field in ("cz-s1","cz-s2","it-pdf","hr-pdf","si-xlsx","si-audit",
                  "height","seed","private-source-csv","aggregate-qa"):
        p.add_argument("--"+field,required=True,type=Path)
    a=p.parse_args()
    run(a.cz_s1,a.cz_s2,a.it_pdf,a.hr_pdf,a.si_xlsx,a.si_audit,
        a.height,a.seed,a.private_source_csv,a.aggregate_qa)
if __name__=="__main__":
    main()
