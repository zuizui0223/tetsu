#!/usr/bin/env python3
"""v29.4: source-verified 4 European railway floras, CLASSIFICATION-SCALE falsification.

A genuine source original study/site/species matrix is FROZEN. Compare
original v28.5 broad 2D trait relabeling, v29.3 within-genus relabeling,
and a new GBIF family-CONSTRAINED bivariate relabeling. GBIF is a source
taxonomy proxy, NOT a calibrated molecular plant phylogeny.

No unmatched genus is guessed: family-unknown taxa remain immobile and
source geography/species composition never changes. Report predeclared
source taxonomy coverage/exchangeability; if insufficient, NOT IDENTIFIABLE.
"""
from __future__ import annotations
import argparse, concurrent.futures, datetime, hashlib, json, math, re, time
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

import numpy as np
from validate_Europe_original_4rail_studies_genus_trait_label_null_v293 import (
    STUDIES, N, METRIC_KEYS, publisher_czech, publisher_italy, publisher_hr_ba,
    extract_original_xlsx, frozen_reference, build_matrix, metric_values,
    source_genus_null, original_species_and_genus_index,
    FROZEN_V285_QUINTILE_HEIGHT, FROZEN_V285_QUINTILE_JOINT, LEDASHA,
)
from audit_Europe_4independent_railstudies_fixed_community_trait_label_null_v285 import (
    source_study_null, MODES
)
from audit_SI_CZ_original_continuous_trait_coverage_v220 import read_leda

GBIF_ENDPOINT="https://api.gbif.org/v1/species/match"
GBIF_USER_AGENT="tetsu-ecology-taxonomic-family-qa/29.4 (research source verifications)"
MIN_ACCEPTED_CONFIDENCE=90
MAX_GENUS_LOOKUPS=1200
MAX_GBIF_WORKERS=4
TAXONOMY_TIMEOUT=25
MIN_FAMILY_ASSIGNED_ORIGINAL_INCIDENCE_FRACTION=.70
MIN_MULTI_GENUS_FAMILIES=8
MIN_FAMILY_EXCHANGEABLE_INCIDENCE_FRACTION=.30
MIN_SUPPORTED_SOURCE_GEOGRAPHY_GROUPS=2
DRAWS_PER_CHAIN=299
SEEDS=(294119,294120)
METRICS=("cwm2D_mean_within_original_geography",
         "cwm_height_absolute_mean_within_original_geography",
         "cwm_seed_mass_absolute_mean_within_original_geography")

def validate_gbif_genus_taxonomy(genus_name:str, data:dict)->dict:
    """No fuzzy/higher-rank or ambiguous matches promoted to confident family."""
    if not isinstance(data,dict):return {"status":"NOT_IDENTIFIABLE_MALFORMED"}
    checks={
       "matchType":data.get("matchType")=="EXACT",
       "rank":data.get("rank")=="GENUS",
       "confidence":isinstance(data.get("confidence"),(float,int)) and
                    data.get("confidence",0)>=MIN_ACCEPTED_CONFIDENCE,
       "kingdom":data.get("kingdom")=="Plantae",
       "genus":str(data.get("genus") or data.get("canonicalName") or "").strip()==genus_name,
       "family":bool(re.fullmatch("[A-Za-z][A-Za-z-]{3,}",str(data.get("family") or ""))),
    }
    if not all(checks.values()):
        return {"status":"NOT_IDENTIFIABLE_TAXONOMY","failed_rules":[k for k,v in checks.items() if not v]}
    return {"status":"EXACT_ACCEPTED_GENUS_FAMILY","family":str(data["family"]),
            "familyKey":data.get("familyKey"),
            "confidence":data["confidence"],"usageKey":data.get("usageKey"),
            "matchType":data.get("matchType")}

def lookup_genus_family(genus_name:str)->dict:
    url=GBIF_ENDPOINT+"?"+urlencode({
        "name":genus_name, "rank":"GENUS", "kingdom":"Plantae", "strict":"true"})
    err=None
    for attempt in range(3):
        try:
            req=Request(url,headers={"Accept":"application/json",
                                      "User-Agent":GBIF_USER_AGENT})
            with urlopen(req,timeout=TAXONOMY_TIMEOUT) as response:
                if response.status!=200:raise ValueError("GBIF genus match HTTP unexpected")
                raw=response.read(60_000)
                if len(raw)>=60_000:raise ValueError("GBIF response oversized")
                obj=json.loads(raw)
            return validate_gbif_genus_taxonomy(genus_name,obj)
        except HTTPError as e:
            err=e
            if e.code not in (429,500,502,503,504):break
            time.sleep((attempt+1)*1.5)
        except (URLError,TimeoutError,ValueError,json.JSONDecodeError) as e:
            err=e
            time.sleep(.5*(attempt+1))
    return {"status":"ERROR_EXTERNAL_TAXONOMY_API","error_class":type(err).__name__}

def lookup_study_genera(genera:set[str],offline:dict|None=None)->dict:
    if not genera or len(genera)>MAX_GENUS_LOOKUPS:
        raise ValueError("Source genuine genus count is empty or exceeds fixed bounded classification lookup")
    if offline is not None:
        # Unit-test mode. The real workflow must pass None and call first-party API.
        if set(offline)!=set(genera):raise ValueError("Synthetic taxonomy map genus support mismatch")
        return dict(offline)
    names=sorted(genera)
    out={}
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_GBIF_WORKERS) as pool:
        future={pool.submit(lookup_genus_family,g):g for g in names}
        for task in concurrent.futures.as_completed(future):
            g=future[task]
            try:out[g]=task.result()
            except Exception as err:
                out[g]={"status":"ERROR_EXTERNAL_TAXONOMY_API",
                        "error_class":type(err).__name__}
    if set(out)!=set(genera):raise ValueError("Incomplete public genus taxonomy lookup")
    return out

def family_index(species:list[str],taxonomy:dict)->dict[str,list[int]]:
    groups=defaultdict(list)
    for i,taxon in enumerate(species):
        g=taxon.split()[0]
        info=taxonomy.get(g,{})
        if info.get("status")=="EXACT_ACCEPTED_GENUS_FAMILY":
            family=info["family"]
            groups["FAMILY:"+family].append(i)
        else:
            # unresolved taxa do not receive invented families or swap.
            groups["UNRESOLVED:"+taxon].append(i)
    if sum(map(len,groups.values()))!=len(species):
        raise ValueError("Original plant species lost from family partition")
    return dict(sorted(groups.items()))

def original_family_source_gate(sites:list,traits:dict,taxonomy:dict)->dict:
    names,genusgroups=original_species_and_genus_index(sites,traits)
    groups=family_index(names,taxonomy)
    known={names[i] for k,idx in groups.items() if k.startswith("FAMILY:") for i in idx}
    multi_generafamilies={
       k for k,idx in groups.items() if k.startswith("FAMILY:")
       and len({names[i].split()[0] for i in idx})>=2}
    movable={names[i] for k in multi_generafamilies for i in groups[k]}
    total=sum(len(s["taxa"]) for s in sites)
    known_inc=sum(len(s["taxa"]&known) for s in sites)
    mobile_inc=sum(len(s["taxa"]&movable) for s in sites)
    by_geo=defaultdict(lambda:[0,0,0])
    for site in sites:
        key=site["group"]
        by_geo[key][0]+=len(site["taxa"]&known)
        by_geo[key][1]+=len(site["taxa"]&movable)
        by_geo[key][2]+=len(site["taxa"])
    source_geo={str(g):{"family_assigned_incidence_fraction":v[0]/v[2],
                        "multi_genus_family_exchangeable_fraction":v[1]/v[2]}
                for g,v in by_geo.items()}
    n_geo=sum(x["multi_genus_family_exchangeable_fraction"]
              >=MIN_FAMILY_EXCHANGEABLE_INCIDENCE_FRACTION
              and x["family_assigned_incidence_fraction"]
              >=MIN_FAMILY_ASSIGNED_ORIGINAL_INCIDENCE_FRACTION
              for x in source_geo.values())
    eligible=(known_inc/total>=MIN_FAMILY_ASSIGNED_ORIGINAL_INCIDENCE_FRACTION and
              len(multi_generafamilies)>=MIN_MULTI_GENUS_FAMILIES and
              mobile_inc/total>=MIN_FAMILY_EXCHANGEABLE_INCIDENCE_FRACTION and
              n_geo>=MIN_SUPPORTED_SOURCE_GEOGRAPHY_GROUPS)
    return {
        "n_original_sites":len(sites),
        "n_original_trait_known_species":len(names),
        "n_original_genera":len(genusgroups),
        "n_assigned_trait_known_taxa":len(known),
        "n_unresolved_trait_known_taxa_frozen":len(names)-len(known),
        "n_supported_source_families":sum(x.startswith("FAMILY:") for x in groups),
        "n_multiple_source_genera_families":len(multi_generafamilies),
        "source_family_assigned_incidence_fraction":known_inc/total,
        "source_exchangeable_multi_genus_family_incidence_fraction":mobile_inc/total,
        "original_source_geography_support":source_geo,
        "n_geographies_meeting_predeclared_support":n_geo,
        "source_family_conditioned_null_eligible":eligible,
        "status":"FAMILY_NULL_IDENTIFIABLE" if eligible
                 else "NOT_IDENTIFIABLE_ORIGINAL_FAMILY_TAXONOMY_SUPPORT"
    }

def shuffle_within_family(T:np.ndarray,groups:dict,rng:np.random.Generator)->np.ndarray:
    out=T.copy()
    if sum(map(len,groups.values()))!=len(T):raise ValueError("Species do not partition")
    for name,idx in groups.items():
        if len(idx)<2 or not name.startswith("FAMILY:"):continue
        perm=rng.permutation(idx)
        out[idx,:]=T[perm,:]
        if sorted(map(tuple,out[idx,:]))!=sorted(map(tuple,T[idx,:])):
            raise ValueError("Paired traits crossed family")
    if sorted(map(tuple,out))!=sorted(map(tuple,T)):
        raise ValueError("Original trait pairs changed across taxonomy null")
    return out

def family_null_one_study(sites,traits,taxonomy,seed_offset=0,draws=DRAWS_PER_CHAIN)->dict:
    gate=original_family_source_gate(sites,traits,taxonomy)
    if not gate["source_family_conditioned_null_eligible"]:
        return {"status":gate["status"],"source_family_support_QA":gate,
                "effects":None,"not_zero_biological_effect":True}
    W,T,geo,qa=build_matrix(sites,traits)
    names,genusgroups=original_species_and_genus_index(sites,traits)
    if len(names)!=T.shape[0]:raise ValueError("Original sorted taxon order mismatch")
    groups=family_index(names,taxonomy)
    obs=metric_values(W,T,qa["i"],qa["j"],qa["within"])
    draws_values={m:[] for m in METRICS}
    for seed in SEEDS:
        rng=np.random.default_rng(seed+seed_offset)
        for _ in range(draws):
            p=shuffle_within_family(T,groups,rng)
            y=metric_values(W,p,qa["i"],qa["j"],qa["within"])
            for m in METRICS:draws_values[m].append(y[m])
    return {"status":"ORIGINAL_SOURCE_FAMILY_NULL_ESTIMATED","source_family_support_QA":gate,
            "source_species_site_incidence_AND_geography_frozen":True,
            "paired_LEDA_trait_vector_always_in_source_family":True,
            "technical_taxon_label_shuffles_per_null":draws*len(SEEDS),
            "effects":{m:{
                "original_observed":obs[m],
                "source_family_null_mean":mean(values),
                "observed_minus_family_trait_null":obs[m]-mean(values),
                "technical_upper_rank_not_field_p":
                    (1+sum(v>=obs[m] for v in values))/(len(values)+1)
            } for m,values in draws_values.items()}}

def run(cz_s1,cz_s2,it_pdf,hr_pdf,si_xlsx,si_audit,height,seed,
        prev_json:Path,draws:int=DRAWS_PER_CHAIN)->dict:
    for file,k in ((height,"height"),(seed,"seed")):
        if hashlib.sha256(file.read_bytes()).hexdigest()!=LEDASHA[k]:
            raise ValueError("Author-authenticated LEDA file SHA drift")
    h,_=read_leda(height,"height_m");s,_=read_leda(seed,"seed_mass_mg")
    published={"CZ_Czechia":publisher_czech(cz_s1,cz_s2),
               "IT_Italy_Parma":publisher_italy(it_pdf),
               "HR_BA_2017_one_study":publisher_hr_ba(hr_pdf)}
    studies={k:v[0] for k,v in published.items()}
    sourceqa={k:v[1] for k,v in published.items()}
    original_si=extract_original_xlsx(si_xlsx,si_audit)
    cohorts,traits,gate=frozen_reference(
        (studies,{**sourceqa,"slovenia":original_si}),{"height":h,"seed":s})
    if set(cohorts)!=set(STUDIES) or sum(len(x) for x in cohorts.values())!=171:
        raise ValueError("Original four European field studies changed")
    prev=json.loads(prev_json.read_text())
    if prev.get("version")!="v29.3" or len(prev.get("studies",[]))!=4:
        raise ValueError("Original source v29.3 reference aggregate missing")
    prior={s["id"]:s for s in prev["studies"]}
    for label in STUDIES:
        if label not in prior:raise ValueError("Source previous study identifier drift")
        old=source_genus_null(cohorts[label],traits[label],label,draws=299)
        if old["status"]!="SOURCE_GENUS_TRAIT_NULL_ESTIMABLE":
            raise ValueError("Previously source-authenticated v29.3 genus null no longer supported")
        h_old=old["source_effects"][METRICS[1]]["observed_minus_genus_label_null"]
        j_old=old["source_effects"][METRICS[0]]["observed_minus_genus_label_null"]
        if (abs(h_old-prior[label]["height_v293_genus_excess"])>1e-9 or
            abs(j_old-prior[label]["joint_v293_excess"])>1e-9):
            raise ValueError("Prior exact original source genus null drift: "+label)
    all_genera={g for label in STUDIES for site in cohorts[label]
                for t in site["taxa"] for g in [t.split()[0]]}
    taxonomy=lookup_study_genera(all_genera)
    error=sum(v["status"]=="ERROR_EXTERNAL_TAXONOMY_API" for v in taxonomy.values())
    exact=sum(v["status"]=="EXACT_ACCEPTED_GENUS_FAMILY" for v in taxonomy.values())
    if error:
        # A network fault is source-acquisition failure, not evidence of no
        # taxonomic family relationship.
        raise RuntimeError(f"Taxonomy reference unavailable for {error} of {len(taxonomy)} original genera")
    digest=hashlib.sha256(json.dumps(taxonomy,sort_keys=True,
                     ensure_ascii=False,separators=(",",":")).encode()).hexdigest()
    result={}
    for label in STUDIES:
        result[label]=family_null_one_study(
            cohorts[label],traits[label],taxonomy,
            seed_offset=17813*STUDIES.index(label),draws=draws)
    valid={k:v for k,v in result.items()
           if v["status"]=="ORIGINAL_SOURCE_FAMILY_NULL_ESTIMATED"}
    summary=None
    if len(valid)>=3:
        summary={}
        for metric in METRICS:
            d={k:v["effects"][metric]["observed_minus_family_trait_null"]
               for k,v in valid.items()}
            loo=[mean(x for k,x in d.items() if k!=name) for name in d]
            summary[metric]={
                "study_equal_mean_observed_minus_family_null":mean(d.values()),
                "n_positive_original_publications":sum(x>0 for x in d.values()),
                "n_original_independent_study_publications":len(d),
                "study_leave_one_out_mean_range":[min(loo),max(loo)],
                "by_study_source_effect":d,
            }
    return {
        "version":"v29.4",
        "source_scope":"FOUR original railway field publications; MIREN roads not included",
        "author_source_real_study_counts":{k:len(cohorts[k]) for k in STUDIES},
        "original_3study_trait_axis_and_source_cohort_QA":gate,
        "source_taxonomy_reference":{
          "provider":"GBIF legacy v1 species name match, Plantae genus exact only",
          "api":GBIF_ENDPOINT,
          "observation_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat(),
          "taxonomy_snapshot_digest_sha256":digest,
          "unique_source_genera":len(taxonomy),
          "exact_GBIF_genus_family_matched":exact,
          "other_invalid_or_ambiguous_without_family":len(taxonomy)-exact,
          "raw_original_species_per_site_and_taxa_NOT_exported":True,
          "not_a_versioned_phylogenetic_backbone":True,
        },
        "per_original_publication_family_null":result,
        "family_conditioned_meta_descriptive_if_supported":summary,
        "comparison_scales":"v28.5 broad across-genus, v29.3 genus-only, v29.4 family-only; source sites and z-log axes fixed, but taxonomic null exchangeability differs",
        "no_cross_model_difference_interpreted_as_phylogenetic_variance_explained":True,
        "family_taxonomy_is_not_molecular_phylogeny":True,
        "not_railway_specific_adaptation_or_causal_infrastructure_treatment":True,
    }
def main():
    p=argparse.ArgumentParser()
    for key in ("cz-s1","cz-s2","it-pdf","hr-pdf","si-xlsx",
                "si-audit","height","seed","previous-v293","output"):
        p.add_argument("--"+key,required=True,type=Path)
    p.add_argument("--draws",type=int,default=DRAWS_PER_CHAIN)
    args=p.parse_args()
    out=run(args.cz_s1,args.cz_s2,args.it_pdf,args.hr_pdf,args.si_xlsx,
            args.si_audit,args.height,args.seed,args.previous_v293,args.draws)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"v":out["version"],"taxonomy":out["source_taxonomy_reference"],
       "studies":{k:{"status":v["status"],"support":v["source_family_support_QA"],
                     "effects":v["effects"]} for k,v in out["per_original_publication_family_null"].items()},
       "summary":out["family_conditioned_meta_descriptive_if_supported"]},indent=2))
if __name__=="__main__":main()
