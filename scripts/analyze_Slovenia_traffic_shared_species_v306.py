#!/usr/bin/env python3
"""v306: WHICH source-identified plants account for Slovenia freight/rail similarity?

Exploratory AFTER observing v305. Read exact publisher Slovenian 2022-23
BS/BT flora and genuine 2020 operator/Eurostat section-matched traffic.
No pseudoreplication of nine traffic segments as 36 independent pairs.
No unmeasured plant traits assigned: 'Bidens' is an author's strategy
archetype, NOT observed physical attachment to a train.
"""
from __future__ import annotations
import argparse, hashlib, json, math, random
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx, validate_original_keys
from analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 import to_original_sites, jaccard_similarity
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import haversine_km
from analyze_Slovenia_historical_train_frequency_local_species_coupling_v275 import mapped_segment_data
from validate_SI_actual_2020_freight_vegetation_beta_holdout_v233 import exact_within_line_labels

SEEDS=(306202,306203)
DRAWS=199
BURN=8000
SPACING=400
ZONES=("BS","BT","union")
WIND=frozenset(("Epilobium","Phragmites"))
TOL=1e-9

def zoneset(site, zone):
    if zone=="union":return set(site["BS"])|set(site["BT"])
    return set(site[zone])

def site_cohort(original,traffic):
    original_sites,siteqa=to_original_sites(original)
    sites=sorted((x for area in original_sites.values() for x in area),
                 key=lambda r:r["source_id"])
    assert len(sites)==89 and len({s["source_id"] for s in sites})==89
    by={r["source_id"]:i for i,r in enumerate(sites)}
    segments=defaultdict(list)
    for t in traffic:
        sid=str(t["sid"])
        if sid not in by or sites[by[sid]]["section"]!=str(t["section"]):
            raise ValueError("Official 2020 traffic to botanical source ID/section drift")
        segments[t["source_segment"]].append((by[sid],t))
    keys=sorted(segments,key=str)
    if len(keys)!=9 or sum(len(v) for v in segments.values())!=27 or len({v[0][1]["section"] for v in segments.values()})!=5:
        raise ValueError("Nine/27/five genuine historical source field support drift")
    for key in keys:
        if len({(t["freight"],t["passenger"],t["section"]) for _,t in segments[key]})!=1:
            raise ValueError("Real physical source traffic segment assignments disagree")
    return sites,keys,segments,siteqa

def traits_and_names(original):
    l,z,sp,tr,cb,sourceqa=validate_original_keys(original)
    return tr,cb,sourceqa

def site_pair_dyads(sites,keys,groups):
    """Fixed source geography/exposure: never selected using species results."""
    dyads=[]
    for i,ka in enumerate(keys):
        for j in range(i+1,len(keys)):
            kb=keys[j]
            ga,gb=groups[ka],groups[kb]
            members=[(a,b) for a,_ in ga for b,_ in gb]
            geo=[haversine_km(sites[a]["lat"],sites[a]["lon"],
                              sites[b]["lat"],sites[b]["lon"])
                 for a,b in members]
            ta,tb=ga[0][1],gb[0][1]
            dyads.append({
                "left":i,"right":j,"members":members,
                "log_distance":float(np.log1p(np.mean(geo))),
                "freight":float((np.log1p(ta["freight"])+np.log1p(tb["freight"]))/2),
                "passenger":float((np.log1p(ta["passenger"])+np.log1p(tb["passenger"]))/2),
            })
    if len(dyads)!=36:raise ValueError("Original nine segment dyads not recovered")
    return dyads

def metric_and_contributions(source_sets,dyads,taxa):
    """Mean site-pair Jaccard exactly decomposed into source taxon columns."""
    ix={k:j for j,k in enumerate(taxa)}
    vals=np.zeros((len(dyads),len(taxa)))
    for i,d in enumerate(dyads):
        for a,b in d["members"]:
            x,y=source_sets[a],source_sets[b]
            both=x&y; either=x|y
            if not either:raise ValueError("Original BS/BT source pair has empty union")
            amount=1/(len(either)*len(d["members"]))
            for taxon in both:
                vals[i,ix[taxon]]+=amount
    y=vals.sum(axis=1)
    if not np.allclose(y,[np.mean([jaccard_similarity(source_sets[a],source_sets[b])
                                  for a,b in d["members"]]) for d in dyads],
                       atol=1e-12):
        raise ValueError("Species contributions do not sum to genuine original mean Jaccard")
    return y,vals

def exposure_design(dyads,labels=None):
    x=[]
    for d in dyads:
        freight=d["freight"] if labels is None else (labels[d["left"]]+labels[d["right"]])/2
        x.append([1.,freight,d["log_distance"]])
    return np.asarray(x,float)

def fit_slope(dyads,y,labels=None):
    X=exposure_design(dyads,labels)
    if np.linalg.matrix_rank(X)!=3:raise ValueError("True real-dyad regression not identifiable")
    return float(np.linalg.lstsq(X,np.asarray(y,float),rcond=None)[0][1])

def contribution_analysis(sites,dyads,traits,codebook,zone,baseline,year=2020):
    full=[zoneset(r,zone) for r in sites]
    taxa=sorted(set().union(*full))
    y,mat=metric_and_contributions(full,dyads,taxa)
    X=exposure_design(dyads)
    c=np.linalg.lstsq(X,mat,rcond=None)[0][1]
    raw=float(np.linalg.lstsq(X,y,rcond=None)[0][1])
    if abs(sum(c)-raw)>1e-10:raise ValueError("Decomposed species freight slopes do not sum")
    if baseline is not None and abs(raw-baseline)>TOL:
        raise ValueError(f"Frozen first-party v305 original flora/source no-urban coefficient changed: {zone} {raw} != {baseline}")
    def taxoninfo(code,beta):
        t=traits.get(code,{})
        name=t.get("latin_name") or codebook.get(code,{}).get("latin_name")
        status=t.get("Nativeness.status")
        life=t.get("Lifespan.categories")
        strategy=t.get("Diaspore.strategy")
        return {
            "taxon":code,"source_latin_name":name or "codebook-unmatched",
            "contribution_to_no_urban_adjusted_freight_slope":float(beta),
            "status_explicit_alien":bool(isinstance(status,str) and "alien" in status.lower()),
            "original_nativeness_status":status,
            "source_strict_annual":life=="annual",
            "original_lifespan":life,
            "original_diaspore_strategy":strategy,
            "author_wind_archetype":strategy in WIND,
            "author_Bidens_archetype_not_real_train_attachment":strategy=="Bidens",
        }
    infos=[taxoninfo(k,v) for k,v in zip(taxa,c)]
    def groupsum(fn):
        yes=[t for t in infos if fn(t)]
        return {"n_source_species":len(yes),
                "signed_contribution_sum":float(sum(x["contribution_to_no_urban_adjusted_freight_slope"] for x in yes)),
                "positive_contribution_sum":float(sum(max(x["contribution_to_no_urban_adjusted_freight_slope"],0) for x in yes)),
                "negative_contribution_sum":float(sum(min(x["contribution_to_no_urban_adjusted_freight_slope"],0) for x in yes))}
    statusgroups={
        "explicit_alien":groupsum(lambda t:t["status_explicit_alien"]),
        "not_explicitly_alien_INCLUDES_native_unknown_crop":groupsum(lambda t:not t["status_explicit_alien"]),
        "strict_annual":groupsum(lambda t:t["source_strict_annual"]),
        "not_strict_annual_INCLUDES_mixed_unknown":groupsum(lambda t:not t["source_strict_annual"]),
        "author_wind_Epilobium_Phragmites":groupsum(lambda t:t["author_wind_archetype"]),
        "author_Bidens_not_measured_rail_attachment":groupsum(lambda t:t["author_Bidens_archetype_not_real_train_attachment"]),
        "source_missing_nativeness":groupsum(lambda t:t["original_nativeness_status"] is None),
        "source_missing_diaspore_strategy":groupsum(lambda t:t["original_diaspore_strategy"] is None),
    }
    strategy={x:groupsum(lambda t:t["original_diaspore_strategy"]==x)
              for x in sorted({t["original_diaspore_strategy"] for t in infos
                               if t["original_diaspore_strategy"] is not None})}
    high=sorted(range(9),key=lambda k:-np.mean([d["freight"] for d in dyads if k in (d["left"],d["right"])]))[:3]
    # rank actual segments by their true log freight inferred from the same pairs
    high_pairs=[i for i,d in enumerate(dyads) if d["left"] in high and d["right"] in high]
    high_species=mat[high_pairs].mean(axis=0) if high_pairs else np.zeros(len(taxa))
    info_by_code={t["taxon"]:t for t in infos}
    top_high=[{**info_by_code[k], "mean_contribution_to_high3_segment_pair_jaccard":float(high_species[i])}
              for i,k in sorted(enumerate(taxa),key=lambda pair:-high_species[pair[0]])[:15]]
    return {
        "zone":zone,"field_taxa_in_89_original_sites":len(taxa),
        "real_original_mean_dyad_Jaccard":float(y.mean()),
        "no_urban_observed_freight_slope":raw,
        "sum_taxon_contributions_reproduces_v305":True,
        "mean_sitepair_denominator_not_pooled_segment_union":True,
        "top_positive_taxon_contributors":sorted(infos,key=lambda t:-t["contribution_to_no_urban_adjusted_freight_slope"])[:15],
        "top_negative_taxon_contributors":sorted(infos,key=lambda t:t["contribution_to_no_urban_adjusted_freight_slope"])[:10],
        "all_taxon_beta_positive_sum":float(sum(max(float(v),0) for v in c)),
        "all_taxon_beta_negative_sum":float(sum(min(float(v),0) for v in c)),
        "source_defined_guilds":statusgroups,
        "source_original_diaspore_strategy_categories":strategy,
        "three_highest_freight_segments_pair_count":len(high_pairs),
        "top_shared_taxa_in_highest_three_segments":top_high,
    }

def swap_communities(sets,rng,proposals,strata):
    if not strata:raise ValueError("Need source site strata")
    good=[list(g) for g in strata if len(g)>=2]
    accepted=0
    for _ in range(proposals):
        group=good[rng.randrange(len(good))]
        i,j=rng.sample(group,2)
        a=sets[i]-sets[j]
        b=sets[j]-sets[i]
        if not a or not b:continue
        u=rng.choice(tuple(a));v=rng.choice(tuple(b))
        sets[i].remove(u);sets[i].add(v)
        sets[j].remove(v);sets[j].add(u)
        accepted+=1
    return accepted

def frequency_margin_null(sites,dyads,zone,observed,mode="global",draws=DRAWS):
    if draws!=DRAWS:raise ValueError("Predeclared 2x199 technical null only")
    source=[zoneset(r,zone) for r in sites]
    if mode=="global":strata=[list(range(len(sites)))]
    elif mode=="source_year_track_layout":
        d=defaultdict(list)
        for i,s in enumerate(sites):d[(s["year"],s["layout"])].append(i)
        strata=[s for s in d.values() if len(s)>=2]
    else:raise ValueError("Unregistered null mode")
    orig_sizes=[len(x) for x in source]
    orig_prevalence=Counter(t for row in source for t in row)
    results=[];chains=[]
    for chain,seed in enumerate(SEEDS):
        rng=random.Random(seed+101*ZONES.index(zone)+(1 if mode!="global" else 0))
        m=[set(r) for r in source]
        accept=swap_communities(m,rng,BURN,strata)
        values=[]
        for _ in range(draws):
            accept+=swap_communities(m,rng,SPACING,strata)
            if [len(x) for x in m]!=orig_sizes or Counter(x for r in m for x in r)!=orig_prevalence:
                raise ValueError("Genuine species frequency/site richness source margins not preserved")
            y=np.array([np.mean([jaccard_similarity(m[i],m[j]) for i,j in row["members"]])
                        for row in dyads])
            value=fit_slope(dyads,y)
            results.append(value);values.append(value)
        chains.append({"chain":chain+1,"acceptance_fraction":accept/(BURN+SPACING*draws),
                       "mean_null_freight_slope":float(np.mean(values))})
    if len(results)!=2*draws or min(x["acceptance_fraction"] for x in chains)<=0:
        raise ValueError("Source fixed margin null not actually generated")
    result=np.array(results)
    return {
        "mode":mode,"source_full_89_site_taxon_incidence_margins_FROZEN":True,
        "n_technical_correlated_draws_NOT_field_replicates":len(results),
        "observed_no_urban_freight_slope":float(observed),
        "null_mean_freight_slope":float(np.mean(result)),
        "observed_minus_null_mean":float(observed-np.mean(result)),
        "technical_upper_tail_rank_NOT_field_p":float((1+sum(x>=observed for x in results))/(1+len(results))),
        "technical_lower_tail_rank_NOT_field_p":float((1+sum(x<=observed for x in results))/(1+len(results))),
        "technical_null_2p5_97p5_NOT_population_CI":list(map(float,np.quantile(result,[.025,.975]))),
        "chain_diagnostics":chains,
        "source_real_traffic_year_and_geography_UNCHANGED":True,
    }

def finite_original_within_line_traffic_labels(keys,groups,dyads,y):
    meta={i:{"national_line":groups[keys[i]][0][1]["section"],
             "freight_2020":groups[keys[i]][0][1]["freight"]}
          for i in range(len(keys))}
    exact_total=0;coefs=[]
    for n,mapping in exact_within_line_labels(meta):
        if mapping is None:exact_total=n;continue
        coefs.append(fit_slope(dyads,y,mapping))
    if exact_total!=48 or len(coefs)!=48:
        raise ValueError("Exact original 48 source-traffic within-line labels changed")
    actual=fit_slope(dyads,y)
    return {
        "source_original_true_national_lines":len({x["national_line"] for x in meta.values()}),
        "n_permuted_traffic_assignments_within_true_lines":len(coefs),
        "observed_adjusted_for_geographic_distance_slope":actual,
        "technical_two_sided_absolute_rank_NOT_randomization_field_p":
            float(sum(abs(v)>=abs(actual)-1e-12 for v in coefs)/len(coefs)),
        "technical_one_sided_positive_rank_NOT_field_p":
            float(sum(v>=actual-1e-12 for v in coefs)/len(coefs)),
        "source_null_minmax":list(map(float,(min(coefs),max(coefs)))),
        "not_a_2015_2020_independent_study_or_causal_traffic_experiment":True,
    }

def main():
    p=argparse.ArgumentParser()
    for name in ("source-xlsx","source-audit","station-pdf","v305-frozen","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    original=extract_original_xlsx(a.source_xlsx,a.source_audit)
    traffic,trafficqa=mapped_segment_data(original,Path(a.station_pdf))
    sites,keys,groups,siteqa=site_cohort(original,traffic)
    t,cb,tqa=traits_and_names(original)
    dyads=site_pair_dyads(sites,keys,groups)
    v305=json.loads(Path(a.v305_frozen).read_text(encoding="utf8"))
    saved=v305["buffers"]["r1000_inner0"]["traffic_pairs"]["models"]
    result={
      "version":"v306",
      "source":"original Slovenia Šabić et al. 2026 railway plants, official operator 2023 PDF and original Eurostat 2020",
      "original_publisher_supplement_SHA256":hashlib.sha256(Path(a.source_xlsx).read_bytes()).hexdigest(),
      "n_true_original_botanical_sites":len(sites),
      "n_mapped_traffic_sites":sum(map(len,groups.values())),
      "n_true_traffic_segments":len(keys),
      "n_DEpendent_segment_pairs_not_independent":len(dyads),
      "n_original_study_programmes":1,
      "original_source_traffic_crosswalk_QA":trafficqa,
      "status":"ORIGINAL_SOURCE_REPLAY_AND_SOURCE_TAXON_IDENTITIES",
      "original_field_site_species_GPS_not_exported":True,
      "zones":{}
    }
    for zone in ZONES:
        ref=saved[f"{zone}/freight/urban=False"]
        if ref["status"]!="descriptive":
            raise ValueError("Source v305 original no-urban reference not evaluable")
        val=contribution_analysis(sites,dyads,t,cb,zone,ref["coefficients"]["freight"])
        source_sets=[zoneset(s,zone) for s in sites]
        y=np.array([np.mean([jaccard_similarity(source_sets[i],source_sets[j])
                             for i,j in d["members"]]) for d in dyads])
        val["source_full_89_site_frequency_fixed_null"]={
            mode:frequency_margin_null(sites,dyads,zone,val["no_urban_observed_freight_slope"],mode)
            for mode in ("global","source_year_track_layout")
        }
        val["true_national_line_48_freight_label_assignments"]=finite_original_within_line_traffic_labels(
            keys,groups,dyads,y)
        val["previous_v305_urban_adjusted_reference_NOT_new_taxon_decomposition"] = (
            v305["buffers"]["r1000_inner0"]["traffic_pairs"]["models"][
                f"{zone}/freight/urban=True"]["coefficients"]["freight"])
        result["zones"][zone]=val
    output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf8")
    print(json.dumps({
        "version":result["version"],"sites":result["n_true_original_botanical_sites"],
        "mapped_sites":result["n_mapped_traffic_sites"],"segments":result["n_true_traffic_segments"],
        "zones":{z:{
          "slope":x["no_urban_observed_freight_slope"],
          "global_margin_null_excess":x["source_full_89_site_frequency_fixed_null"]["global"]["observed_minus_null_mean"],
          "global_margin_null_technical_upper_rank":x["source_full_89_site_frequency_fixed_null"]["global"]["technical_upper_tail_rank_NOT_field_p"],
          "within_line_traffic_relabel_abs_rank":x["true_national_line_48_freight_label_assignments"]["technical_two_sided_absolute_rank_NOT_randomization_field_p"],
          "alien_contribution":x["source_defined_guilds"]["explicit_alien"]["signed_contribution_sum"],
          "annual_contribution":x["source_defined_guilds"]["strict_annual"]["signed_contribution_sum"],
          "wind_strategy_contribution":x["source_defined_guilds"]["author_wind_Epilobium_Phragmites"]["signed_contribution_sum"],
        } for z,x in result["zones"].items()}
    },ensure_ascii=False,indent=2))
if __name__=="__main__":main()
