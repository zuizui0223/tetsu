#!/usr/bin/env python3
"""v19.8 TRUE original-source publication-holdout railway plant community prediction.

Do guild/trait features transfer more reliably between botanical publications
than taxon identity? Compare effort/trait-completeness-only (negative control),
exact Pladias Raunkiaer life form + Sádlo diaspore-strategy composition,
and full source taxon membership. Train on 3/4 of independent *publication
titles*, score on untouched 1/4, preserving original same-publication x
phytoregion x decade x source measured plot-area rail/road matched BLOCKS.
Each heldout publication has equal weight in the final result; within each
heldout block evaluate the rank concordance between source RAIL-text plots
and source ROAD-text plots (block AUC). No environmental causal inference.
Full source CzechVeg MD5 and Pladias original trait-row gates are required.

A "rail" category means only original source HABITAT/LOCALITY textual clues,
not legally verified active rail right-of-way, train-mediated propagules or
adaptive phenotype. Original source botanical plots may have survey bias;
4 folds and source-paper AUC are NOT replicated randomized interventions.
"""
from __future__ import annotations
import argparse,json,math,random,statistics,sys
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyze_CzechVeg_Pladias_rail_road_diaspore_v119 import (
    verify_file,plot_labels,plot_taxa,parse_traits,
    HEADERS_NAME,SPECIES_NAME
)
from analyze_CzechVeg_rail_road_region_year_area_overlap_v121 import header_source_context
from analyze_CzechVeg_rail_road_same_study_trait_v123 import (
    source_study_meta,common_original_study_blocks,SPECS
)
from validate_CzechVeg_publication_holdout_rail_species_v197 import split_source_studies

SEED=20261116
PANEL_THRESHOLDS=("all_vascular_min5","both_trait50_min5")
MODELS=("effort","guild","taxon")
PERMUTE_GUILD_TRAIN_N=24
MAX_ITER=650
C=0.7
SOURCE_N_PLOTS=117739
SOURCE_N_PLOT_TAXON_ROWS=2510340
MATCHING=next(x for x in SPECS if x[0]=="S2_study_region_decade_area4_min2")

def source_blocks_and_traits(zenodo,header_csv,species_csv,diaspore_xlsx,life_xlsx):
    head=verify_file(zenodo,HEADERS_NAME,header_csv)
    taxonsrc=verify_file(zenodo,SPECIES_NAME,species_csv)
    dispersal,life,traitq=parse_traits(diaspore_xlsx,life_xlsx)
    if traitq["source_Pladias_original_diaspore_rows"]!=3147 or traitq[
        "source_Pladias_original_lifeform_rows"]!=4599:
        raise ValueError("True original Pladias trait row census changed")
    base,labelqa=plot_labels(header_csv)
    year,yearqa=header_source_context(header_csv,base)
    meta,studyqa=source_study_meta(header_csv,year)
    plots,plantqa=plot_taxa(species_csv,meta,dispersal,life)
    if plantqa["n_original_source_species_csv_records_scanned"]!=SOURCE_N_PLOT_TAXON_ROWS:
        raise ValueError("Original 2.51M botanical species rows missing")
    blocks,blockqa=common_original_study_blocks(meta,plots,"broad",MATCHING)
    if len({b[0] for b in blocks})!=36 or len(blocks)!=44:
        raise ValueError("True source original v19.7 44 block/36 publication support changed")
    return blocks,plots,{
        "source_md5_headers":head,
        "source_md5_species":taxonsrc,
        "source_pladias_semantic_qa":traitq,
        "original_matched_frame_44block_36study":blockqa,
        "original_plant_taxon_join":plantqa,
        "original_source_metadata":{"plot":labelqa["counts"],
                                    "source_published_study_title_support":studyqa,
                                    "year_region":yearqa},
    }

def site_characteristics(plant_tuple):
    # Source CSV deduped as distinct taxon×plot in v11.9.
    taxa={x[0] for x in plant_tuple}
    d={x[0]:x[1] for x in plant_tuple if x[1] is not None}
    l={x[0]:x[2] for x in plant_tuple if x[2] is not None}
    n=len(taxa)
    return {"taxa":taxa,"life":l,"disp":d,"n":n,
            "n_life":len(l),"n_disp":len(d),
            "frac_life":len(l)/n if n else 0,
            "frac_disp":len(d)/n if n else 0}

def site_pass(x,panel):
    if x["n"]<5:return False
    if panel=="all_vascular_min5":return True
    if panel=="both_trait50_min5":
        return x["n_life"]>=3 and x["n_disp"]>=3 and x["frac_life"]>=.5 and x["frac_disp"]>=.5
    raise ValueError("Unknown original site source-coverage panel")

def frame_records(blocks,plants,panel):
    characteristics={k:site_characteristics(x) for k,x in plants.items()}
    selected={}
    for key,arms in blocks.items():
        a=[i for i in arms["rail"] if site_pass(characteristics[i],panel)]
        b=[i for i in arms["road"] if site_pass(characteristics[i],panel)]
        if len(a)>=2 and len(b)>=2:
            selected[key]={"rail":a,"road":b}
    site_data={}
    orig={}
    for key,arms in selected.items():
        for name,label in (("rail",1),("road",0)):
            for pid in arms[name]:
                if pid in site_data:raise AssertionError("Source plot used in multiple original matched blocks")
                site_data[pid]=characteristics[pid]
                orig[pid]={"label":label,"study":key[0],"block":key}
    pub=sorted({k[0] for k in selected})
    return selected,site_data,orig,{
        "original_blocks_retained":len(selected),
        "original_publications_retained":len(pub),
        "original_rail_text_plots":sum(len(v["rail"]) for v in selected.values()),
        "original_road_text_plots":sum(len(v["road"]) for v in selected.values()),
        "source_taxon_eligibility_does_not_look_at_group_effect_direction":True,
        "original_site_coverage_panel":panel,
        "heldout_publications_never_overlap_training":True,
        "loss_from_original_44_study_blocks":44-len(selected),
    }

def features(site,mode):
    n=site["n"]
    if n<=0:raise ValueError("Original plot with no vascular taxon")
    d={
        "effort:log_n_taxa":math.log1p(n),
        "coverage:life":site["frac_life"],
        "coverage:diaspore":site["frac_disp"],
    }
    if mode=="effort":return d
    if mode=="guild":
        lc=Counter(site["life"].values())
        dc=Counter(site["disp"].values())
        if lc:
            for label,count in lc.items():d["life:"+label]=count/sum(lc.values())
        if dc:
            for label,count in dc.items():d["disp:"+label]=count/sum(dc.values())
    elif mode=="taxon":
        for t in site["taxa"]:
            d["species:"+t]=1.0
    else:
        raise ValueError("Unknown model")
    return d

def study_balanced_train_weights(keys,y,original):
    # Equal total observation mass per source publication, half source-rail,
    # half road. Scale to total train plot count for stable penalization.
    paper_labels=Counter((original[pid]["study"],int(label)) for pid,label in zip(keys,y))
    nstudies=len({original[pid]["study"] for pid in keys})
    if nstudies<2:raise ValueError("Need multiple original source papers")
    weight_scale=len(keys)/nstudies
    return [weight_scale*.5/paper_labels[(original[pid]["study"],int(label))]
            for pid,label in zip(keys,y)]

def fit_model(keys,original,site_data,mode,override_labels=None):
    from sklearn.feature_extraction import DictVectorizer
    from sklearn.linear_model import LogisticRegression
    if not keys:raise ValueError("No training source plots")
    y=[int(original[k]["label"]) for k in keys] if override_labels is None else [
        int(override_labels[k]) for k in keys]
    if len(set(y))!=2:raise ValueError("Both source arms required for heldout model")
    vectorizer=DictVectorizer(sparse=True,sort=True)
    x=vectorizer.fit_transform([features(site_data[k],mode) for k in keys])
    model=LogisticRegression(C=C,penalty="l2",solver="liblinear",max_iter=MAX_ITER,
                             random_state=SEED,class_weight=None)
    weights=study_balanced_train_weights(keys,y,original)
    model.fit(x,y,sample_weight=weights)
    return model,vectorizer,{"n_training_plots":len(keys),
                              "n_training_nonzero_features":x.shape[1],
                              "n_training_original_papers":len(set(original[k]["study"] for k in keys)),
                              "source_study_balanced_and_arm_balanced_weights":True,
                              "feature_dictionary_fit_only_on_training_source_studies":True}

def source_block_auc(positives,negatives):
    # Identical to within-original-publisher block ROC AUC (no cross-study
    # case/control score calibration); each original study aggregated equally.
    if not positives or not negatives:raise ValueError("Absent railway or road group in original block")
    win=sum((1 if p>q else .5 if p==q else 0) for p in positives for q in negatives)
    return win/(len(positives)*len(negatives))

def evaluate_within_original_papers(model,vectorizer,mode,test_blocks,site_data):
    test_ids=sorted({p for v in test_blocks.values() for p in v["rail"]+v["road"]})
    score=model.decision_function(vectorizer.transform([features(site_data[k],mode)
                                                         for k in test_ids]))
    byid=dict(zip(test_ids,map(float,score)))
    bypub=defaultdict(list)
    for block,v in test_blocks.items():
        auc=source_block_auc([byid[p] for p in v["rail"]],
                             [byid[p] for p in v["road"]])
        bypub[block[0]].append(auc)
    studies={paper:statistics.mean(vals) for paper,vals in sorted(bypub.items())}
    return studies,{
        "n_original_publications_heldout":len(studies),
        "n_matched_original_blocks_heldout":len(test_blocks),
        "n_actual_original_site_records_heldout":len(test_ids),
        "one_paper_one_vote_regardless_of_how_many_study_plots":True,
    }

def permutation_of_training_labels(train_blocks,original,rng):
    # Permute rails/roads within each pre-existing study/region/decade/area
    # block. Original source-paper and within-block case counts conserved.
    label={}
    for key,arms in sorted(train_blocks.items()):
        siteids=sorted(arms["rail"]+arms["road"])
        signs=[1]*len(arms["rail"])+[0]*len(arms["road"])
        rng.shuffle(signs)
        label.update(zip(siteids,signs))
        assert sum(label[p] for p in siteids)==len(arms["rail"])
    return label

def run_fold(blocks,sites,original,folds,fold_idx,permutations=PERMUTE_GUILD_TRAIN_N):
    hold=folds[fold_idx]
    train={k:v for k,v in blocks.items() if k[0] not in hold}
    test={k:v for k,v in blocks.items() if k[0] in hold}
    assert {k[0] for k in train}.isdisjoint({k[0] for k in test})
    if not train or not test:
        return {"identified":False,"reason":"No study-overlap source blocks available in fold"}
    keys=sorted({k for v in train.values() for k in v["rail"]+v["road"]})
    recorded={}
    perpaper={}
    for mode in MODELS:
        model,vec,trainmeta=fit_model(keys,original,sites,mode)
        scores,testmeta=evaluate_within_original_papers(model,vec,mode,test,sites)
        perpaper[mode]=scores
        recorded[mode]={"study_balanced_auc":statistics.mean(scores.values()),
                        "n_heldout_study":len(scores),
                        "n_heldout_study_auc_gt_half":sum(a>.5 for a in scores.values()),
                        "train_qa":trainmeta,"test_qa":testmeta}
    paired={
        "guild_minus_effort_auc":statistics.mean(
            perpaper["guild"][study]-perpaper["effort"][study]
            for study in perpaper["guild"]),
        "taxon_minus_effort_auc":statistics.mean(
            perpaper["taxon"][study]-perpaper["effort"][study]
            for study in perpaper["taxon"]),
        "guild_minus_taxon_auc":statistics.mean(
            perpaper["guild"][study]-perpaper["taxon"][study]
            for study in perpaper["guild"])
    }
    # A label-exchange null for interpretability of model overfitting, not
    # a randomized railway treatment or a valid paper-level significance test.
    random_null=[]
    rng=random.Random(SEED+fold_idx*877)
    for _ in range(permutations):
        y=permutation_of_training_labels(train,original,rng)
        clf,vec,_=fit_model(keys,original,sites,"guild",override_labels=y)
        pp,_=evaluate_within_original_papers(clf,vec,"guild",test,sites)
        random_null.append(statistics.mean(pp.values()))
    return {
        "identified":True,"fold":fold_idx+1,
        "heldout_paper_n":len(perpaper["effort"]),
        "heldout_source_publication_hashes":{
          str(study):__import__("hashlib").sha256(study.encode()).hexdigest()[:12]
          for study in sorted(perpaper["effort"])},
        "source_publication_holdout_model_AUC":recorded,
        "within_identical_study_paired_model_AUC_differences":paired,
        "synthetic_within_source_block_train_label_permutations_n":permutations,
        "technical_train_label_permutation_guild_auc":{
          "mean":statistics.mean(random_null),"minimum":min(random_null),
          "maximum":max(random_null),
          "fraction_permuted_guild_auc_ge_actual":sum(x>=recorded["guild"]["study_balanced_auc"]
                                                     for x in random_null)/len(random_null)},
        "heldout_per_paper_AUC_nonindependent_comparisons":{
          __import__("hashlib").sha256(k.encode()).hexdigest()[:12]:
              {m:round(perpaper[m][k],6) for m in MODELS}
          for k in sorted(perpaper["effort"])},
        "heldout_source_papers_unequal_survey_effort":True,
        "all_three_models_same_exact_heldout_source_sites":True,
        "no_field_p_value_or_verified_RAIL_site_natural_habitat":True,
    }

def run(zenodo,header_csv,species_csv,diaspore_xlsx,life_xlsx):
    blocks,plants,provenance=source_blocks_and_traits(
        zenodo,header_csv,species_csv,diaspore_xlsx,life_xlsx)
    # Fold assignment remains fixed at the original 44-block study source
    # opportunity, regardless of later trait/site coverage gate.
    originalfolds,source_fold_masses=split_source_studies(blocks)
    panels={}
    for panel in PANEL_THRESHOLDS:
        chosen,sites,original,frame=frame_records(blocks,plants,panel)
        foldresults=[run_fold(chosen,sites,original,originalfolds,i)
                     for i in range(len(originalfolds))]
        valid=[x for x in foldresults if x["identified"]]
        study_rows={}
        for f in valid:
            for key,vals in f["heldout_per_paper_AUC_nonindependent_comparisons"].items():
                if key in study_rows:raise AssertionError("Original publication appears in more than one test fold")
                study_rows[key]=vals
        pooled={
            m:statistics.mean(x[m] for x in study_rows.values())
            for m in MODELS} if study_rows else {}
        paired={
            k:statistics.mean(x[k] for x in (
                 {"guild_minus_effort":v["guild"]-v["effort"],
                  "taxon_minus_effort":v["taxon"]-v["effort"],
                  "guild_minus_taxon":v["guild"]-v["taxon"]}
                 for v in study_rows.values()))
            for k in ("guild_minus_effort","taxon_minus_effort","guild_minus_taxon")
        } if study_rows else {}
        panels[panel]={
          "retained_source_frame":frame,
          "n_heldout_original_publications_with_pairwise_AUC":len(study_rows),
          "total_original_blocks_with_source_traits":len(chosen),
          "folds":foldresults,
          "study_equal_mean_heldout_AUC":pooled,
          "paired_heldout_AUC_differences":paired,
          "n_heldout_source_publications_guild_AUC_greater_than_effort":
            sum(v["guild"]>v["effort"] for v in study_rows.values()),
          "n_heldout_source_publications_taxon_AUC_greater_than_effort":
            sum(v["taxon"]>v["effort"] for v in study_rows.values()),
          "train_and_test_original_paper_split_opportunity_only":True,
        }
    return {
        "version":"v19.8_source_paper_holdout_functional_vs_species",
        "original_publisher":"CzechVeg-Open Chytry et al 2026; original Pladias Sadlo et al 2018",
        "source_zenodo_record":"10.5281/zenodo.19231152",
        "source":provenance,
        "study_folds_opportunity_balanced_no_source_taxon_leak":{
            "study_n":sum(len(x) for x in originalfolds),
            "source_study_opportunity_fold_masses":source_fold_masses,
            "original_study_fold_sizes":[len(x) for x in originalfolds]},
        "original_study_heldout_predictability":panels,
        "source_rail_road_habitat_is_TEXT_label_only":True,
        "publishing_studies_are_not_ran_random_sample_of_Europe":True,
        "MCMC_or_training_permutation_samples_not_field_p_values":True,
        "original_vegetation_and_trait_raw_source_not_redistributed":True,
        "substantial_Pl adias_use_before_publication_contact_board".replace(" ",""):True
    }

def main():
    p=argparse.ArgumentParser()
    for fld in ("zenodo-record","source-headers","source-species",
                "source-diaspore-xlsx","source-lifeform-xlsx","output"):
        p.add_argument("--"+fld,required=True)
    a=p.parse_args()
    res=run(json.loads(Path(a.zenodo_record).read_text(encoding="utf8")),
            a.source_headers,a.source_species,
            a.source_diaspore_xlsx,a.source_lifeform_xlsx)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(res,ensure_ascii=False,indent=2),encoding="utf8")
    short={}
    for panel,v in res["original_study_heldout_predictability"].items():
        short[panel]={
            "source_study_n":v["n_heldout_original_publications_with_pairwise_AUC"],
            "source_frame":v["retained_source_frame"],
            "source_mean_paper_AUC":v["study_equal_mean_heldout_AUC"],
            "model_pair_differences":v["paired_heldout_AUC_differences"],
            "source_studies_positive_guild_increment":
               v["n_heldout_source_publications_guild_AUC_greater_than_effort"],
            "source_studies_positive_species_increment":
               v["n_heldout_source_publications_taxon_AUC_greater_than_effort"],
            "source_folds":[{"fold":f.get("fold"),
                "identified":f["identified"],
                "original_test_papers":f.get("heldout_paper_n"),
                "study_auc":{k:vv["study_balanced_auc"]
                        for k,vv in f.get("source_publication_holdout_model_AUC",{}).items()},
                "label_permuted_guild_null":
                    f.get("technical_train_label_permutation_guild_auc")} for f in v["folds"]]}
    print(json.dumps({"version":res["version"],"source_full_CzechVeg":
      {k:v for k,v in provenance_summary(res["source"]).items()},
      "fixed_original_source_study_fold":res["study_folds_opportunity_balanced_no_source_taxon_leak"],
      "heldout_report":short},ensure_ascii=False,indent=2),flush=True)

def provenance_summary(provenance):
    return {"CzechVeg_header_MD5":provenance["source_md5_headers"]["checksum"],
            "CzechVeg_species_MD5":provenance["source_md5_species"]["checksum"],
            "Pladias_lifeform_original_rows":provenance[
                "source_pladias_semantic_qa"]["source_Pladias_original_lifeform_rows"],
            "Pladias_diaspore_original_rows":provenance[
                "source_pladias_semantic_qa"]["source_Pladias_original_diaspore_rows"],
            "source_original_44_matched_blocks":provenance[
                "original_matched_frame_44block_36study"][
                    "n_original_study_blocks_with_rail_and_road_common_support"]}

if __name__=="__main__":
    main()
