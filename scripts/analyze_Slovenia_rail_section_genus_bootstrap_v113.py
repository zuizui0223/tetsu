#!/usr/bin/env python3
"""v11.3: 15-rail-section and source taxonomic clustering sensitivity.

Original published Slovenia 2026 railway botanical field source, Šabić,
Grašič & Rovanšek DOI 10.1007/s44473-026-00155-z. 89 distinct actual
field locations, 173 true same-side BS-BT original pairs, 15 true railway
sections. Strict wind label from original Diaspore.strategy. No genuine
freight seed/wind/mowing mechanism measured.

The v11.0 source result has SAME 261 original taxa with opposite
species-equal vs occurrence-opportunity weighted wind guild contrasts.
This script checks if the sign inversion survives:
  (1) 15-source-section Bayesian positive-weight bootstrap, preserving
      ALL original taxa in every draw (no disappearance of rare species);
  (2) positive-weight resampling of genuine ORIGINAL genera;
  (3) leave-one-genus-out original plant trait-label sensitivity.
The intervals describe within-published-source heterogeneity, NOT
geographically representative global confidence intervals or a causal
inference test.
"""
from __future__ import annotations
import argparse,json,math,random,re
from collections import defaultdict,Counter
from pathlib import Path
from statistics import mean

from analyze_Slovenia_rail_microhabitats_v103 import (
    SOURCE_SHA256,SOURCE_DOI,extract_original_xlsx,validate_original_keys
)
from analyze_Slovenia_rail_species_occupancy_permutation_v109 import (
    original_location_matched_sets,source_guild
)

DRAWS=2500
GENUS_DRAWS=2500
SEED=20261114
STRATA=("strict_annual","other_known_lifespan")


def percentile(x,p):
    x=sorted(x)
    if not x:return None
    z=(len(x)-1)*p; lo=int(math.floor(z)); hi=int(math.ceil(z))
    return x[lo] if lo==hi else x[lo]+(x[hi]-x[lo])*(z-lo)


def weighted_taxon_records(sites,traits,section_weights=None):
    """Every original rail section must have strictly positive mass."""
    sections=sorted({v["section_number"] for v in sites.values()})
    sw=section_weights if section_weights is not None else {k:1.0 for k in sections}
    if set(sw)!=set(sections) or any(not math.isfinite(v) or v<=0 for v in sw.values()):
        raise ValueError("Each original railway section requires a positive finite weight")
    acc=defaultdict(lambda:[0.,0.,0.])
    for loc in sites.values():
        weight=sw[loc["section_number"]]
        bs=loc["BS"];bt=loc["BT"]
        for code in sorted(bs|bt):
            vals=acc[code]
            if code in bs and code in bt:vals[2]+=weight
            elif code in bs:vals[0]+=weight
            else:vals[1]+=weight
    rows=[]
    for code in sorted(acc):
        lab=source_guild(traits[code])
        if lab is None:continue
        a,b,c=acc[code];total=a+b+c
        if not(total>0):raise ValueError("Original observed taxon lost in positive-mass bootstrap")
        rows.append({"taxon_code":code,"wind":lab["wind"],"stratum":lab["stratum"],
                     "occupancy_preference_BS_minus_BT":(a-b)/total,
                     "n_original_locations_present_either_zone":total})
    return rows


def taxon_name_map(original,trait_keys):
    """Only actual original published scientific binomials, never guess EPPO genus."""
    candidates={}
    for sheet in ("species_codebook","traits"):
        for row in original[sheet]:
            code=row.get("EPPO_code")
            if code not in trait_keys:continue
            for key,v in row.items():
                if key=="EPPO_code" or not isinstance(v,str):continue
                val=v.strip()
                if not re.match(r"^[A-Z][a-z]+\s+(?:[a-z]+|sect\.|sp\.)\b",val):
                    continue
                if re.search(r"family|growth|habitat|life|strategy|native",key,re.I):
                    continue
                genus=val.split()[0]
                if len(genus)>=3:
                    candidates.setdefault(code,[]).append((sheet,key,genus,val))
    genus={}
    mismatches=[]
    for code,choices in candidates.items():
        unique={x[2] for x in choices}
        if len(unique)!=1:
            mismatches.append({"EPPO_code":code,"orig_recorded_genera":sorted(unique)})
            continue
        genus[code]=next(iter(unique))
    missing=sorted(set(trait_keys)-set(genus))
    return genus,{
        "name_source_keys_found":sorted(set(k for rec in candidates.values()
                                            for _,k,_,_ in rec)),
        "n_trait_coded_taxa_with_genus":len(genus),
        "n_source_trait_codes_no_unambiguous_binomial":len(missing),
        "n_source_trait_codes_disagree_genus":len(mismatches),
        "disagreement_taxon_codes":mismatches[:12],
        "missing_taxon_codes":missing[:15],
        "source_codebook_lothe_doyph_mismatch_not_silently_resolved":True,
    }


def guild_contrast(rows,taxon_weights=None,require_min=3):
    """Exact v11.0 contrast at positive weights, with lifespan strata kept."""
    acc={s:{"wind":[],"other":[]} for s in STRATA}
    if taxon_weights is None:taxon_weights={}
    for t in rows:
        wt=taxon_weights.get(t["taxon_code"],1.0)
        if not math.isfinite(wt) or wt<0:raise ValueError("Negative or nonfinite taxon weight")
        if wt==0:continue
        acc[t["stratum"]]["wind" if t["wind"] else "other"].append((t,wt))
    if any(min(len(acc[s]["wind"]),len(acc[s]["other"]))<require_min for s in STRATA):
        return None
    out={}
    for mode in ("species","opportunities"):
        total=0.;value=0.
        for s in STRATA:
            guild=acc[s]
            sums={}
            for name,items in guild.items():
                denom=sum(w*(1.0 if mode=="species" else t["n_original_locations_present_either_zone"])
                          for t,w in items)
                numer=sum(t["occupancy_preference_BS_minus_BT"]*
                          w*(1.0 if mode=="species" else t["n_original_locations_present_either_zone"])
                          for t,w in items)
                sums[name]=(denom,numer/denom)
            ns=sums["wind"][0]+sums["other"][0]
            value+=ns*(sums["wind"][1]-sums["other"][1])
            total+=ns
        out[mode]=value/total
    out["difference_species_minus_occupancy"]=out["species"]-out["opportunities"]
    out["sign_reversal"]=out["species"]>0 and out["opportunities"]<0
    out["source_taxa_count"]=sum(len(v) for groups in acc.values() for v in groups.values())
    return out


def positive_dirichlet_weights(keys,rng):
    vals=[rng.expovariate(1.0) for _ in keys]
    scale=len(vals)/sum(vals)
    return {k:val*scale for k,val in zip(keys,vals)}


def draws_summary(draws):
    if not draws:return {"n_draws":0,"not_identifiable":True}
    return {
        "n_draws":len(draws),
        "species_equal_median":round(percentile([r["species"] for r in draws],.5),6),
        "species_equal_95_source_cluster_interval":[round(percentile([r["species"] for r in draws],p),6)
                                                   for p in (.025,.975)],
        "opportunity_weighted_median":round(percentile([r["opportunities"] for r in draws],.5),6),
        "opportunity_weighted_95_source_cluster_interval":[
            round(percentile([r["opportunities"] for r in draws],p),6) for p in (.025,.975)],
        "n_draws_retain_species_plus_occupancy_minus":sum(r["sign_reversal"] for r in draws),
        "fraction_draws_retain_sign_reversal":round(mean(r["sign_reversal"] for r in draws),6),
        "n_draws_both_positive":sum(r["species"]>=0 and r["opportunities"]>=0 for r in draws),
        "n_draws_both_negative":sum(r["species"]<=0 and r["opportunities"]<=0 for r in draws),
        "no_claim_probability_global_RAIL_sign_reversal":True,
    }


def run(original,verified110):
    if (verified110.get("complete") is not True
        or verified110.get("version")!="v11.0"
        or verified110.get("original_source_XLSX_sha256")!=SOURCE_SHA256):
        raise ValueError("Must provide genuine source verified original v11.0 receipt")
    sites,traits,integrity=original_location_matched_sets(original)
    base=weighted_taxon_records(sites,traits)
    if len(base)!=261:
        raise ValueError("Original 261 source taxa were not reproducible")
    direct=guild_contrast(base)
    if direct is None:raise ValueError("Original source taxon guilds not identified")
    for key,target in (
        ("species","full_1plus_species_equal_wind_minus_other_preference"),
        ("opportunities","full_1plus_occupancy_weighted_wind_minus_other_preference")):
        if not math.isclose(direct[key],verified110[target],abs_tol=1e-12):
            raise ValueError("Original v11.0 real source contrast not reproduced "+key)
    original_sections=sorted({v["section_number"] for v in sites.values()})
    if len(original_sections)!=15 or len(sites)!=89:
        raise ValueError("Original 15-section/89-site field sample changed")
    rng=random.Random(SEED)
    sectiondraws=[]
    for _ in range(DRAWS):
        w=positive_dirichlet_weights(original_sections,rng)
        rows=weighted_taxon_records(sites,traits,w)
        if len(rows)!=261:raise ValueError("Positive-mass section resampling lost source taxa")
        d=guild_contrast(rows)
        if d is None:raise ValueError("Source species label groups lost during section bootstrap")
        sectiondraws.append(d)
    genusmap,genusaudit=taxon_name_map(original,{t["taxon_code"] for t in base})
    # Source taxonomy must be sufficient for genus-cluster analysis.
    known=[r for r in base if r["taxon_code"] in genusmap]
    distinct_genera=sorted(set(genusmap[x["taxon_code"]] for x in known))
    if len(known)<220 or len(distinct_genera)<120:
        raise ValueError("Original scientific names cannot support genus clustering")
    taxonomic_base=guild_contrast(known)
    if taxonomic_base is None:raise ValueError("Missing strict wind/other in source taxonomic strata")
    genusdraws=[]
    for _ in range(GENUS_DRAWS):
        w=positive_dirichlet_weights(distinct_genera,rng)
        tw={k:w[genusmap[k]] for k in genusmap}
        d=guild_contrast(known,tw)
        if d is not None:genusdraws.append(d)
    genusloo=[]
    for g in distinct_genera:
        tw={t["taxon_code"]:float(genusmap[t["taxon_code"]]!=g) for t in known}
        d=guild_contrast(known,tw)
        genusloo.append({"omitted_original_genus":g,"source_taxa_removed":
                         sum(genusmap[t["taxon_code"]]==g for t in known),
                         "species_equal":None if d is None else round(d["species"],6),
                         "opportunity_weighted":None if d is None else round(d["opportunities"],6),
                         "retains_original_opposite_signs":False if d is None else d["sign_reversal"]})
    return {
        "complete":True,"version":"v11.3",
        "original_source_DOI":SOURCE_DOI,"original_supplement_sha256":SOURCE_SHA256,
        "original_field_locations":len(sites),"original_railway_sections":len(original_sections),
        "original_matched_track_side_pairs":sum(v["n_matched_sides"] for v in sites.values()),
        "original_source_taxa_both_lifespan_and_diaspore":len(base),
        "source_original_contrasts_reproduced":direct,
        "positive_weight_section_bootstrap_all_261_taxa":draws_summary(sectiondraws),
        "source_taxonomic_name_audit":genusaudit,
        "source_taxa_with_genus_for_sensitivity":len(known),
        "n_original_genera_with_source_names":len(distinct_genera),
        "original_taxonomy_known_taxa_only":taxonomic_base,
        "original_genus_resampling":draws_summary(genusdraws),
        "n_leave_one_genus_out_identified":sum(x["species_equal"] is not None for x in genusloo),
        "n_leave_one_genus_out_retain_opposite_signs":sum(x["retains_original_opposite_signs"]
                                                        for x in genusloo),
        "leave_one_genus_out_examples_not_all_species_records":genusloo,
        "major_caveat":(
            "Positive-weight Bayesian section/genus bootstrap has NO phylogeny model and does "
            "NOT deliver globally representative confidence intervals. All observed species "
            "have strictly positive sampling mass so original rare taxa never disappear. "
            "Taxonomic genus is a correlated-species proxy, NOT independent multi-country "
            "replication. Cannot prove railway-mediated propagule movement, wind mortality, "
            "maintenance selection, phenotypic plasticity or adaptation."
        ),
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("original-xlsx","source-audit102","verified110","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    source=extract_original_xlsx(a.original_xlsx,a.source_audit102)
    receipt=json.loads(Path(a.verified110).read_text(encoding="utf-8"))
    result=run(source,receipt)
    f=Path(a.output);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({
        "complete":True,"original_source_taxa":result["original_source_taxa_both_lifespan_and_diaspore"],
        "original":result["source_original_contrasts_reproduced"],
        "sections":result["positive_weight_section_bootstrap_all_261_taxa"],
        "taxonomy":result["source_taxonomic_name_audit"],
        "genera":result["n_original_genera_with_source_names"],
        "genus_bootstrap":result["original_genus_resampling"],
        "genus_loo_retains":result["n_leave_one_genus_out_retain_opposite_signs"]
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()
