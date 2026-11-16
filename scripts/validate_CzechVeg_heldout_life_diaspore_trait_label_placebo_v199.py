#!/usr/bin/env python3
"""v19.9: which externally observed plant traits transfer to unseen railway publications?

Continuation of source-verified CzechVeg/Pladias study-holdout v19.8.
Break down the v19.8 heldout functional signal into lifeform-only and
diaspore-only features, both using exactly the same actual source survey
sites, original study folds, outcome labels, baseline controls and training
model. Compare REAL trait mapping against matched synthetic LABEL SHUFFLES.

The shuffled "placebo" holds actual source taxon incidence, identity,
frequency, missing-trait masks, global marginal trait counts and
joint lifeform/dispersal correlation for dual-known taxa fixed; it only
shuffles which taxon receives each trait combination, either globally or
within source-occupancy-frequency strata. This is an internal
phenotypic-annotation mapping null, NOT a field randomized treatment.
Rail/road are text descriptions, not confirmed right-of-way or managed
trackbed. No inferred viable train seed transport or evolutionary selection.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random,statistics,sys
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from validate_CzechVeg_source_publication_holdout_guild_vs_taxon_v198 import (
    source_blocks_and_traits,frame_records,study_balanced_train_weights,
    source_block_auc,split_source_studies,PANEL_THRESHOLDS,SEED,C,MAX_ITER
)
N_PERMUTATIONS=32
MODELS=("effort","lifeform","dispersal","combined")
SHUFFLE_DESIGNS=("global","occupancy_stratified")

def feature(site,mode,trait_map=None):
    x={
        "effort:log_n_taxa":math.log1p(site["n"]),
        "coverage:life":site["frac_life"],
        "coverage:diaspore":site["frac_disp"],
    }
    if mode=="effort":return x
    if mode not in MODELS:raise ValueError("Unknown source trait ablation")
    l_map=trait_map["life"] if trait_map else None
    d_map=trait_map["disp"] if trait_map else None
    if mode in ("lifeform","combined"):
        count=Counter((l_map[t] if l_map is not None else label)
                      for t,label in site["life"].items())
        if count:
            for name,n in count.items():
                x["life:"+name]=n/sum(count.values())
    if mode in ("dispersal","combined"):
        count=Counter((d_map[t] if d_map is not None else label)
                      for t,label in site["disp"].items())
        if count:
            for name,n in count.items():
                x["disp:"+name]=n/sum(count.values())
    return x

def source_heldout_auc(blocks,sites,original,study_folds,mode,trait_map=None):
    from sklearn.feature_extraction import DictVectorizer
    from sklearn.linear_model import LogisticRegression
    published={}
    fold_reports=[]
    for i,unseen in enumerate(study_folds):
        train={key:value for key,value in blocks.items() if key[0] not in unseen}
        test={key:value for key,value in blocks.items() if key[0] in unseen}
        if not train or not test:
            fold_reports.append({"fold":i+1,"identified":False})
            continue
        train_ids=sorted({p for b in train.values() for p in b["rail"]+b["road"]})
        test_ids=sorted({p for b in test.values() for p in b["rail"]+b["road"]})
        assert set(train_ids).isdisjoint(test_ids)
        assert {original[p]["study"] for p in train_ids}.isdisjoint(
            {original[p]["study"] for p in test_ids})
        y=[original[p]["label"] for p in train_ids]
        vectorizer=DictVectorizer(sparse=True,sort=True)
        training_x=vectorizer.fit_transform([feature(sites[p],mode,trait_map) for p in train_ids])
        model=LogisticRegression(C=C,penalty="l2",solver="liblinear",
                                 max_iter=MAX_ITER,random_state=SEED,
                                 class_weight=None)
        weights=study_balanced_train_weights(train_ids,y,original)
        model.fit(training_x,y,sample_weight=weights)
        score=model.decision_function(
            vectorizer.transform([feature(sites[p],mode,trait_map) for p in test_ids]))
        pscore=dict(zip(test_ids,map(float,score)))
        paper_auc=defaultdict(list)
        for key,b in test.items():
            paper_auc[key[0]].append(source_block_auc(
                [pscore[p] for p in b["rail"]],
                [pscore[p] for p in b["road"]]
            ))
        auc_per_paper={
            hashlib.sha256(paper.encode()).hexdigest()[:12]:statistics.mean(vals)
            for paper,vals in paper_auc.items()}
        if set(auc_per_paper)&set(published):raise AssertionError("Paper was heldout twice")
        published.update(auc_per_paper)
        fold_reports.append({
            "fold":i+1,"identified":True,
            "n_heldout_papers":len(auc_per_paper),
            "mean_original_paper_equal_auc":statistics.mean(auc_per_paper.values()),
            "n_source_train_only_model_features":training_x.shape[1],
            "n_original_heldout_source_site_records":len(test_ids)})
    if not published:
        return {"identified":False,"reason":"No evaluable heldout original publication"}
    return {"identified":True,"n_source_original_heldout_publications":len(published),
            "paper_balanced_heldout_AUC":statistics.mean(published.values()),
            "per_fold_original_source":fold_reports,
            "original_source_paper_AUC_by_hash":published,
            "original_source_studies_not_derived_from_independent_country_samples":True}

def original_trait_maps(sites):
    life={};disp={};frequency=Counter()
    for site in sites.values():
        frequency.update(site["taxa"])
        for taxon,v in site["life"].items():
            if taxon in life and life[taxon]!=v:
                raise ValueError("Pladias lifeform changed for original full taxon")
            life[taxon]=v
        for taxon,v in site["disp"].items():
            if taxon in disp and disp[taxon]!=v:
                raise ValueError("Pladias diaspore strategy changed for original full taxon")
            disp[taxon]=v
    return life,disp,frequency

def occupancy_bin(n):
    if n<=1:return "one"
    if n<=3:return "2-3"
    if n<=7:return "4-7"
    if n<=15:return "8-15"
    return "16+"

def shuffled_taxon_mapping(life,disp,frequency,rng,design):
    if design not in SHUFFLE_DESIGNS:raise ValueError("Invalid source counterfactual")
    taxa=sorted(set(life)|set(disp))
    strata=defaultdict(list)
    for t in taxa:
        status=("both" if t in life and t in disp else
                "life" if t in life else "disp")
        bin_=occupancy_bin(frequency[t]) if design=="occupancy_stratified" else "any"
        strata[(status,bin_)].append(t)
    lcopy={};dcopy={}
    changed=0
    for key,names in sorted(strata.items()):
        # Shuffle correlated pair of lifeform and dispersal together when
        # both known; keep exact known/missing status of every source taxon.
        status,_=key
        if status=="both":
            vals=[(life[t],disp[t]) for t in names]
            rng.shuffle(vals)
            for t,(a,b) in zip(names,vals):
                lcopy[t]=a;dcopy[t]=b
                changed+=int((life[t],disp[t])!=(a,b))
        elif status=="life":
            vals=[life[t] for t in names];rng.shuffle(vals)
            for t,a in zip(names,vals):
                lcopy[t]=a;changed+=int(life[t]!=a)
        else:
            vals=[disp[t] for t in names];rng.shuffle(vals)
            for t,a in zip(names,vals):
                dcopy[t]=a;changed+=int(disp[t]!=a)
    if set(lcopy)!=set(life) or set(dcopy)!=set(disp):
        raise AssertionError("Pladias source trait-missingness coverage changed under null")
    if Counter((lcopy[t],dcopy[t]) for t in set(lcopy)&set(dcopy)) != Counter(
          (life[t],disp[t]) for t in set(life)&set(disp)):
        raise AssertionError("Paired Pladias lifeform/diaspore distribution altered by shuffling")
    if Counter(lcopy.values())!=Counter(life.values()) or Counter(dcopy.values())!=Counter(disp.values()):
        raise AssertionError("Published original traits marginal category counts altered by null")
    return {"life":lcopy,"disp":dcopy},{"n_taxa_with_nontrivial_trait_label_change":changed,
                                      "original_taxa_count_known_either_trait":len(taxa),
                                      "original_occupancy_bin_strata_count":len(strata),
                                      "pladias_trait_coverage_unchanged":True}

def analyze_panel(blocks,plants,study_folds,panel,permutations=N_PERMUTATIONS):
    selected,sites,original,source_q=frame_records(blocks,plants,panel)
    original_auc={mode:source_heldout_auc(
        selected,sites,original,study_folds,mode) for mode in MODELS}
    if any(not r["identified"] for r in original_auc.values()):
        return {"identified":False,"original_source":source_q,
                "reason":"Not enough source publication coverage"}
    paper_set={frozenset(r["original_source_paper_AUC_by_hash"]) for r in original_auc.values()}
    if len(paper_set)!=1:raise AssertionError("Guild and baseline evaluated on different source publications")
    perstudy={}
    for paper in sorted(next(iter(original_auc.values()))["original_source_paper_AUC_by_hash"]):
        perstudy[paper]={
            mode:original_auc[mode]["original_source_paper_AUC_by_hash"][paper]
            for mode in MODELS}
    paired={
        a+"_minus_"+b:statistics.mean(v[a]-v[b] for v in perstudy.values())
        for a,b in (("combined","effort"),("lifeform","effort"),
                    ("dispersal","effort"),("combined","lifeform"),
                    ("combined","dispersal"),("lifeform","dispersal"))
    }
    l,d,freq=original_trait_maps(sites)
    placebo={}
    for design in SHUFFLE_DESIGNS:
        rng=random.Random(SEED+1729*SHUFFLE_DESIGNS.index(design)+
                          len(sites)*7+len(l))
        results=[];changed=[]
        for j in range(permutations):
            mapping,info=shuffled_taxon_mapping(l,d,freq,rng,design)
            shuffled=source_heldout_auc(selected,sites,original,study_folds,"combined",mapping)
            results.append(shuffled["paper_balanced_heldout_AUC"])
            changed.append(info["n_taxa_with_nontrivial_trait_label_change"])
        vals=sorted(results)
        observed=original_auc["combined"]["paper_balanced_heldout_AUC"]
        placebo[design]={
            "n_synthetic_global_taxon_trait_label_reassignments":permutations,
            "placebo_mean_source_study_AUC":statistics.mean(vals),
            "placebo_q025_technical_not_field_CI":vals[int(.025*len(vals))],
            "placebo_q975_technical_not_field_CI":vals[min(len(vals)-1,int(.975*len(vals)))],
            "placebo_min":vals[0],"placebo_max":vals[-1],
            "observed_guild_AUC":observed,
            "observed_minus_placebo_mean":observed-statistics.mean(vals),
            "fraction_placebo_AUC_ge_observed_not_p_value":
                sum(v>=observed for v in vals)/len(vals),
            "average_source_taxa_truly_reassigned_trait_labels":statistics.mean(changed),
            "trait_missingness_mask_originally_preserved":True,
            "original_trait_marginal_and_joint_distribution_preserved":True,
            "phylogenetic_relationships_and_site_management_not_preserved_or_identified":True
        }
    return {
        "identified":True,
        "original_source_panel":source_q,
        "same_source_publications_and_plot_ids_for_all_four_models":True,
        "n_original_source_heldout_publications":len(perstudy),
        "original_study_equal_holdout_model_auc":{
            k:v["paper_balanced_heldout_AUC"] for k,v in original_auc.items()},
        "study_equal_paired_AUC_differences":paired,
        "source_paper_AUC_positive_increment_counts":{
           a+"_better_than_"+b:sum(v[a]>v[b] for v in perstudy.values())
           for a,b in (("combined","effort"),("lifeform","effort"),
                       ("dispersal","effort"),("lifeform","dispersal"))
        },
        "original_study_AUC_by_hash":perstudy,
        "trait_source_taxa_counts":{
            "lifeform":len(l),"dispersal":len(d),
            "joint_known":len(set(l)&set(d)),
            "n_unique_taxa_in_retained_source_plots":len(freq)},
        "study_folds_original_publication":original_auc["combined"]["per_fold_original_source"],
        "actual_trait_mapping_vs_Pladias_label_shuffled_nulls":placebo,
        "no_field_p_value_or_crosscountry_independent_test":True,
        "source_habitat_is_original_site_TEXT_not_verified_true_ballast":True,
    }

def run(zenodo,header,species,diaspore,life_xlsx):
    blocks,plants,source=source_blocks_and_traits(
        zenodo,header,species,diaspore,life_xlsx)
    folds,opportunity=split_source_studies(blocks)
    panels={k:analyze_panel(blocks,plants,folds,k) for k in PANEL_THRESHOLDS}
    return {
        "version":"v19.9_original_pladias_functional_ablation_and_label_null",
        "source_full_czechveg_and_Pladias_QA":source,
        "original_44_source_blocks_and_36_study_folds":{
            "n_source_studies":sum(len(x) for x in folds),
            "n_source_folds":len(folds),
            "study_opportunity_mass":opportunity},
        "original_study_out_of_source_guild_vs_taxon_mapping":panels,
        "source_trait_labels_SHUFFLED_within_known_census_not_raw_records":True,
        "source_railroad_source_text_not_confirmed_active_ballast":True,
        "no_phenotype_genetic_adaptation_train_seed_transport_or_causal_management":True
    }

def main():
    p=argparse.ArgumentParser()
    for fld in ("zenodo-record","headers","species","diaspore-xlsx","lifeform-xlsx","output"):
        p.add_argument("--"+fld,required=True)
    a=p.parse_args()
    r=run(json.loads(Path(a.zenodo_record).read_text(encoding="utf8")),
          a.headers,a.species,a.diaspore_xlsx,a.lifeform_xlsx)
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding="utf8")
    short={}
    for name,x in r["original_study_out_of_source_guild_vs_taxon_mapping"].items():
        short[name]={
            "identified":x["identified"],
            "n_source_studies":x.get("n_original_source_heldout_publications"),
            "source_original_plot_frame":x.get("original_source_panel"),
            "heldout_AUC":x.get("original_study_equal_holdout_model_auc"),
            "heldout_AUC_ablation_differences":x.get("study_equal_paired_AUC_differences"),
            "source_paper_increment_counts":x.get("source_paper_AUC_positive_increment_counts"),
            "real_vs_shuffled_traits":x.get("actual_trait_mapping_vs_Pladias_label_shuffled_nulls")
        }
    print(json.dumps({
        "version":r["version"],
        "source_44blocks":r["original_44_source_blocks_and_36_study_folds"],
        "panels":short},ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":main()
