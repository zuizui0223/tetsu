#!/usr/bin/env python3
"""v11.9: source-verified original CzechVeg railway vs roadside plant trait JOIN.

Original independent CzechVeg-Open v1.1 117739 vegetation plot headers,
original source vegetation plot species table, CC BY 4.0
10.5281/zenodo.19231152. Join exact published species names to freely
downloadable Pladias 2018 Sádlo Diaspore.strategy archetypes and
Raunkiaer life forms (public usage rules apply).

Contrasts describe RAIL TEXT CANDIDATE plots vs ROAD TEXT CONTROL plots,
NOT original rail-track BS-versus-BT exact area matched pairs. Each plot
is a genuine original phytosociological relevé and its own source SURF_AREA.
No source label is imputed based solely on distance to railway, reference
title or presence of rare species. A source plot's true rail right-of-way
membership and survey-area treatment are not independently geovalidated.

No causal train-transported seeds, gust/wind mortality, herbicide mowing,
flower timing, biomass, plasticity or genetically evolved local adaptation.
No raw Pladias trait source files or full CzechVeg source re-exported.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re,unicodedata
from collections import Counter,defaultdict
from pathlib import Path
from statistics import mean,median

from audit_CzechVeg_rail_habitat_original_v117 import (
  DOI,EXPECTED_PLOTS,HEADERS_NAME,original_zenodo_entry,
  STRONG_RAIL,ROAD,CONFOUND_BUS
)
from audit_Pladias_public_diaspore_lifeform_v118 import (
  xlsx_rows,DISP_CODES
)
from analyze_Slovenia_rail_species_occupancy_weighting_v110 import group_estimator

SPECIES_NAME="CzechVeg_Open_species_v1-1_2026-03-26.csv"
N_SRC_RAIL_HABITAT_EXPECTED=200 # initial scan found ~200 candidate mentions
MIN_PLOT_SIZE_M2=1
MAX_PLOT_SIZE_M2=1000


def verify_file(meta,name,source_file):
    file=original_zenodo_entry(meta,name)
    size=Path(source_file).stat().st_size
    if size!=file["source_bytes"]:
        raise ValueError("Original CzechVeg file bytes differs from Zenodo published source")
    digest=hashlib.new(file["checksum_algorithm"])
    with open(source_file,"rb") as f:
        for chunk in iter(lambda:f.read(2**20),b""):digest.update(chunk)
    if digest.hexdigest()!=file["checksum"]:
        raise ValueError("Original CzechVeg vegetation file md5/sha mismatch")
    return {"source_key":file["source_key"],"checksum":file["checksum"],
            "checksum_algorithm":file["checksum_algorithm"],
            "source_bytes":size,
            "original_publisher_file_url":file["official_url"]}


def normalize_taxon(x):
    x=unicodedata.normalize("NFKC",str(x or "")).strip()
    x=re.sub(r"\s+"," ",x).replace("×","x ")
    x=re.sub(r"\s+"," ",x).strip()
    return x


def parse_traits(diaspore_file,life_file):
    rawd=Path(diaspore_file).read_bytes()
    rawl=Path(life_file).read_bytes()
    _,fd,rd=xlsx_rows(rawd)
    _,fl,rl=xlsx_rows(rawl)
    if len(rd)!=3147 or len(rl)!=4599:
        raise ValueError("Original Pladias 3147 strategy and 4599 life form source rows changed")
    if "Taxon" not in fd or "Value [enum]" not in fd or "Taxon" not in fl:
        raise ValueError("Published Pladias source trait taxonomy columns not recognized")
    source_d=defaultdict(set)
    source_l=defaultdict(lambda:{"all":set(),"dominant":set()})
    for row in rd:
        name=normalize_taxon(row["Taxon"])
        value=str(row.get("Value [enum]") or "").strip()
        label=value.split(" ",1)[0]
        if name and label in DISP_CODES:
            source_d[name].add(label)
    for row in rl:
        name=normalize_taxon(row["Taxon"])
        value=str(row.get("Value [enum]") or "").strip().lower()
        dominance=str(row.get("Dominance [boolean]") or "").strip().lower()
        if name and value:
            source_l[name]["all"].add(value)
            if dominance in ("pravda","true","1","ano","yes"):
                source_l[name]["dominant"].add(value)
    good_d={k:next(iter(v)) for k,v in source_d.items() if len(v)==1}
    conflict_d=set(k for k,v in source_d.items() if len(v)>1)
    good_l={}
    ambiguous_l=[]
    for name,d in source_l.items():
        candidate=d["dominant"] if d["dominant"] else d["all"]
        if len(candidate)==1:good_l[name]=next(iter(candidate))
        elif len(candidate)>1:ambiguous_l.append(name)
    return good_d,good_l,{
      "source_Pladias_original_diaspore_rows":len(rd),
      "source_Pladias_original_lifeform_rows":len(rl),
      "source_strategy_SHA256":hashlib.sha256(rawd).hexdigest(),
      "source_lifeform_SHA256":hashlib.sha256(rawl).hexdigest(),
      "n_unique_strategy_taxon_names":len(source_d),
      "n_unique_strict_unambiguous_strategy_names":len(good_d),
      "n_source_strategy_conflicting_taxa_excluded":len(conflict_d),
      "n_source_lifeform_unambiguous_names":len(good_l),
      "n_source_lifeform_ambiguous_taxa_excluded":len(ambiguous_l),
      "source_dispersal_code_counts":dict(Counter(good_d.values())),
      "source_rights_require_Pl adias_board_contact_for_substantial_research".replace(" ",""):
          True,
    }


def positive_area(v):
    try:x=float(v)
    except (ValueError,TypeError):return None
    if not math.isfinite(x) or x<=0:return None
    return x


def plot_labels(header_file):
    plots={}
    stats=Counter()
    geographic=Counter()
    original_ids=set()
    with open(header_file,encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        cols=set(reader.fieldnames or [])
        expected={"RELEVE_NR","SURF_AREA","HABITAT","LOCALITY",
                  "SKALICKY1988","YEAR"}
        if expected-cols:raise ValueError("Original full CzechVeg plot header schema not recognized")
        for row in reader:
            stats["original_header_rows"]+=1
            plot=str(row.get("RELEVE_NR") or "").strip()
            if not plot or plot in original_ids:
                raise ValueError("Original CzechVeg plot identifier duplicate or missing")
            original_ids.add(plot)
            a=positive_area(row.get("SURF_AREA"))
            h=str(row.get("HABITAT") or "")
            l=str(row.get("LOCALITY") or "")
            rail_h=bool(STRONG_RAIL.search(h)) and not bool(CONFOUND_BUS.search(h))
            rail_s=bool(STRONG_RAIL.search(h+" "+l)) and not bool(CONFOUND_BUS.search(h+" "+l))
            road_h=bool(ROAD.search(h))
            road_s=bool(ROAD.search(h+" "+l))
            # Treat multi-vector site description as neither exclusive treatment.
            strict="rail" if rail_h and not road_h else (
                "road" if road_h and not rail_h else None)
            broad="rail" if rail_s and not road_s else (
                "road" if road_s and not rail_s else None)
            if strict:stats["n_strict_"+strict]+=1
            if broad:stats["n_broad_"+broad]+=1
            if a is None:stats["n_plot_source_area_missing"]+=1
            if a is None or not(MIN_PLOT_SIZE_M2<=a<=MAX_PLOT_SIZE_M2):
                if strict or broad:stats["n_labeled_sites_excluded_out_of_area_scope"]+=1
                continue
            source_year=str(row.get("YEAR") or "").strip()
            region=str(row.get("SKALICKY1988") or "").strip()
            plots[plot]={"strict":strict,"broad":broad,
                         "area":a,"year":source_year,"region":region}
            if strict:stats["n_eligible_strict_"+strict]+=1
            if broad:stats["n_eligible_broad_"+broad]+=1
            if strict=="rail":geographic[region]+=1
    if stats["original_header_rows"]!=EXPECTED_PLOTS:
        raise ValueError("Unmodified CzechVeg full original 117739 plots missing")
    return plots,{"counts":dict(stats),"n_original_valid_area_plot_id_records":len(plots),
                  "strict_rail_original_region_counts":dict(geographic),
                  "n_original_geographical_regions_strict_rail":sum(x>0 for x in geographic.values())}


def plot_taxa(source_species_csv,eligible,disp,life):
    """Stream full original CSV. Count species once per original plot.
    Do not store full source original dataset in output or GitHub."""
    retained=defaultdict(set)
    scanned=0
    matched_raw=Counter()
    with open(source_species_csv,encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        fields=set(reader.fieldnames or [])
        need={"RELEVE_NR","KAPLAN_PLUS_NAME_SIMPLE","PLANT_GROUP"}
        if need-fields:
            raise ValueError("Original CzechVeg species CSV lacks mandated source species plot keys")
        for r in reader:
            scanned+=1
            plot=str(r.get("RELEVE_NR") or "").strip()
            if plot not in eligible:continue
            meta=eligible[plot]
            if not meta["strict"] and not meta["broad"]:continue
            if str(r.get("PLANT_GROUP") or "").strip().lower()!="vascular":continue
            matched_raw["selected_original_source_vascular_plant_rows"]+=1
            name=normalize_taxon(r.get("KAPLAN_PLUS_NAME_SIMPLE"))
            if not name or name.lower() in ("na","nan","null"):continue
            retained[plot].add(name)
    pairs=set((plot,name) for plot,names in retained.items() for name in names)
    if not pairs:raise ValueError("Source original source plot names from actual CzechVeg species table unavailable")
    source_by_class=Counter()
    source_by_type=Counter()
    source_by_trait=Counter()
    source_taxa={}
    for plot,name in sorted(pairs):
        attr=eligible[plot]
        d=disp.get(name)
        l=life.get(name)
        for typ in ("strict","broad"):
            grp=attr[typ]
            if grp:
                source_by_class[(typ,grp)]+=1
                if d is not None:source_by_type[(typ,grp)]+=1
                if d is not None and l is not None:source_by_trait[(typ,grp)]+=1
        source_taxa.setdefault(plot,[]).append((name,d,l))
    return source_taxa,{
        "n_original_source_species_csv_records_scanned":scanned,
        "n_original_source_vascular_rows_at_rail_road_labeled_plots":
            matched_raw["selected_original_source_vascular_plant_rows"],
        "n_unique_plot_by_taxon_after_de_dupe":len(pairs),
        "source_unique_plot_taxon_counts":{f"{a}_{b}":n for (a,b),n in source_by_class.items()},
        "source_exact_Pladias_strategy_covered_plot_taxon_counts":
            {f"{a}_{b}":n for (a,b),n in source_by_type.items()},
        "source_exact_Pladias_strategy_plus_lifeform_plot_taxon_counts":
            {f"{a}_{b}":n for (a,b),n in source_by_trait.items()},
        "no_fuzzy_name_matches_or_genus_inferred_traits":True,
    }


def guild_panel(eligible,plot_taxon,typ,subset=None):
    if subset is None:subset=set(eligible)
    candidate={k:v for k,v in eligible.items() if k in subset and v[typ] in ("rail","road")}
    nrail=sum(v[typ]=="rail" for v in candidate.values())
    nroad=sum(v[typ]=="road" for v in candidate.values())
    info={"n_original_rail_candidate_plots":nrail,"n_original_road_control_plots":nroad}
    if min(nrail,nroad)<10:
        return {**info,"identified":False,"cause":"fewer than 10 source plot opportunities in one original class"}
    taxa=defaultdict(lambda:{"rail":set(),"road":set(),"strategy":None,"life":None})
    for plot,meta in candidate.items():
        for name,d,l in plot_taxon.get(plot,[]):
            if d is None or l is None:continue
            x=taxa[name]
            if x["strategy"] not in (None,d) or x["life"] not in (None,l):
                raise ValueError("Taxon source strategy or lifeform drifted within one name")
            x["strategy"]=d;x["life"]=l
            x[meta[typ]].add(plot)
    rows=[]
    for name,rec in sorted(taxa.items()):
        cr=len(rec["rail"]);cd=len(rec["road"])
        if cr+cd==0:continue
        rows.append({
            "taxon_code":name,
            "wind":rec["strategy"] in ("Epilobium","Phragmites"),
            # Raunkiaer therophyte is a PROXY for strict annual only.
            "stratum":"strict_annual" if rec["life"]=="therophyte" else "other_known_lifespan",
            "n_original_locations_present_either_zone":cr+cd,
            "occupancy_preference_BS_minus_BT":cr/nrail-cd/nroad,
        })
    e=group_estimator(rows,1,"species")
    w=group_estimator(rows,1,"opportunities")
    r={**info,"n_source_unique_joint_trait_taxa_in_original_rail_or_road_vegetation":len(rows),
       "n_source_joint_trait_taxa_by_wind_label":dict(Counter(
           "Epilobium+Phragmites" if v["wind"] else "other" for v in rows)),
       "exact_name_life_and_diaspore_coverage":True,
       "identified":bool(e["identified"] and w["identified"]),
       "critical_scope":"Rail/text vs road/text candidate plot occurrence prevalence, NOT matched BS/BT within-railway zones. Therophyte life form is not identical to source strict annual lifespan."
    }
    if r["identified"]:
        r["strict_wind_minus_other_species_equal_rail_vs_road_prevalence_gap"]=e[
            "wind_minus_other_preference_lifespan_adjusted"]
        r["strict_wind_minus_other_occupancy_weighted_rail_vs_road_prevalence_gap"]=w[
            "wind_minus_other_preference_lifespan_adjusted"]
        r["opposite_weighting_signs"]=(
            r["strict_wind_minus_other_species_equal_rail_vs_road_prevalence_gap"]*
            r["strict_wind_minus_other_occupancy_weighted_rail_vs_road_prevalence_gap"]<0)
        r["original_joint_trait_source_lifespan_classes"]=e["by_stratum"]
    return r


def region_balanced_checks(plotmeta,plot_taxon,typ):
    # Region is from ORIGINAL botanical phyto-geographic source.
    regions=defaultdict(lambda:{"rail":[],"road":[]})
    for plot,meta in plotmeta.items():
        label=meta[typ]
        reg=meta["region"]
        if label in ("rail","road") and reg not in ("","NA"):
            regions[reg][label].append(plot)
    regions_eligible={code:vals for code,vals in regions.items()
                      if len(vals["rail"])>=5 and len(vals["road"])>=5}
    scores=[]
    for region,vals in sorted(regions_eligible.items()):
        group=set(vals["rail"]+vals["road"])
        test=guild_panel(plotmeta,plot_taxon,typ,group)
        scores.append({"source_phytogeographic_region":region,
                       "nrail":len(vals["rail"]),"nroad":len(vals["road"]),
                       "identified":test["identified"],
                       "source_equal":test.get("strict_wind_minus_other_species_equal_rail_vs_road_prevalence_gap"),
                       "source_occupancy":test.get("strict_wind_minus_other_occupancy_weighted_rail_vs_road_prevalence_gap"),
                       "sign_reversal":test.get("opposite_weighting_signs")})
    valid=[x for x in scores if x["identified"]]
    return {
        "n_original_phyto_regions_with_at_least_5_rail_and_5_road_candidate_plots":len(scores),
        "n_original_regions_with_guild_comparison_identified":len(valid),
        "n_original_regions_with_sign_reversal":sum(x["sign_reversal"] for x in valid),
        "n_regions_wind_equal_species_more_rail_than_other":sum(x["source_equal"]>0 for x in valid),
        "n_regions_wind_occupancy_more_rail_than_other":sum(x["source_occupancy"]>0 for x in valid),
        "source_region_comparisons_not_pairwise_causal_rail_vs_road":scores,
    }


def run(zenodo_meta,header_csv,species_csv,disp,life,receipt117,receipt118):
    if not(receipt117.get("complete") is True and receipt117.get("version")=="v11.7"
           and receipt118.get("complete") is True and receipt118.get("version")=="v11.8"):
        raise ValueError("Require original independent CzechVeg and Pladias source receipts")
    head=verify_file(zenodo_meta,HEADERS_NAME,header_csv)
    origsp=verify_file(zenodo_meta,SPECIES_NAME,species_csv)
    strategies,lifestrata,traits_audit=parse_traits(disp,life)
    for key,expected in (
        ("source_strategy_SHA256",receipt118["strategy_original_source_metadata"]["source_file_sha256"]),
        ("source_lifeform_SHA256",receipt118["lifeform_original_source_metadata"]["source_file_sha256"])):
        if traits_audit[key]!=expected:
            raise ValueError("Original public Pladias trait source SHA differs from verified v11.8")
    plotmeta,plot_audit=plot_labels(header_csv)
    taxa,species_audit=plot_taxa(species_csv,plotmeta,strategies,lifestrata)
    primary=guild_panel(plotmeta,taxa,"strict")
    broad=guild_panel(plotmeta,taxa,"broad")
    phyto_strict=region_balanced_checks(plotmeta,taxa,"strict")
    phyto_broad=region_balanced_checks(plotmeta,taxa,"broad")
    return {
        "complete":True,"version":"v11.9",
        "original_vegetation_DOI":"10.5281/zenodo.19231152",
        "original_Pladias_DOI":"10.23855/preslia.2021.001",
        "original_diaspore_taxonomy_source":"Sadlo et al. (2018) Preslia 90:1-22",
        "original_source_Czechveg_header_provenance":head,
        "original_source_Czechveg_species_provenance":origsp,
        "original_source_Pladias_trait_quality":traits_audit,
        "original_source_site_label_and_area_audit":plot_audit,
        "original_source_species_plot_and_trait_join_quality":species_audit,
        "high_precision_original_HABITAT_text_label_RAIL_vs_ROAD":primary,
        "broader_source_site_HABITAT_or_LOCALITY_label_RAIL_vs_ROAD":broad,
        "region_stratified_high_precision_source_comparison":phyto_strict,
        "region_stratified_broader_source_comparison":phyto_broad,
        "important_limits":[
            "This is NOT the same paired BS/BT design as Slovenia. Railway candidate plots and road candidate plots were collected for different source studies and not allocated or observed randomly.",
            "No claim any text-labeled candidate is definitively inside active rail right of way; original habitat label needs manual source/GIS audit.",
            "Area distributions, region-year source effort, site vegetation type and habitat differences are not fully matched in this exploratory screen.",
            "Exact Pladias original Latin binomial matches only; taxonomic unknowns and nonmatches are excluded, NOT assigned nonwind.",
            "Original Raunkiaer therophyte classification is a proxy and NOT an exact source strict annual lifespan measurement.",
            "Source occurrence count is not biomass or direct propagated seed arrival.",
            "Pladias substantial research-use governance requires Board contact before manuscript, and original published CzechVeg source authors must be cited.",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    for n in ("zenodo-record","header-csv","species-csv","dispersal-xlsx",
              "lifeform-xlsx","verified117","verified118","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    meta=json.loads(Path(a.zenodo_record).read_text(encoding="utf-8"))
    r117=json.loads(Path(a.verified117).read_text(encoding="utf-8"))
    r118=json.loads(Path(a.verified118).read_text(encoding="utf-8"))
    result=run(meta,a.header_csv,a.species_csv,a.dispersal_xlsx,
               a.lifeform_xlsx,r117,r118)
    dst=Path(a.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({
        "complete":True,"version":"v11.9",
        "header":result["original_source_Czechveg_header_provenance"],
        "species":result["original_source_Czechveg_species_provenance"],
        "traits":result["original_source_Pladias_trait_quality"],
        "plot_audit":result["original_source_site_label_and_area_audit"],
        "source_taxon_match":result["original_source_species_plot_and_trait_join_quality"],
        "strict":result["high_precision_original_HABITAT_text_label_RAIL_vs_ROAD"],
        "broad":result["broader_source_site_HABITAT_or_LOCALITY_label_RAIL_vs_ROAD"],
        "strict_regions":result["region_stratified_high_precision_source_comparison"],
        "broad_regions":result["region_stratified_broader_source_comparison"]
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()
