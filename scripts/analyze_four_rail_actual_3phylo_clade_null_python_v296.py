#!/usr/bin/env python3
"""v29.6 AUTHOR-ORIGINAL 4 railway plant communities x genuine 3 V.PhyloMaker2
branch-length trees. Pure Python, NO R or apt-get; no grafted phylogenetic tips.

The original site×taxon incidence is fixed, original matched LE DA height/seed
vectors are paired and only reallocated within observed original tree clades.
Study fields are independent only at FOUR ORIGINAL PUBLICATIONS, not 171 sites
or three taxonomic tree alternatives. All taxa/site source rows stay in runner
private TEMP folder, never in GitHub artifacts.
"""
from __future__ import annotations
import argparse,csv,json,math,re
from collections import Counter,defaultdict,deque
from pathlib import Path
from statistics import mean
import numpy as np
from audit_VPhyloMaker2_RData_python_no_R_environment_v296 import (
    inspect_tree_rdata,TREES,COMMIT
)
from prepare_Europe_4original_rail_phylogeny_audit_source_v295 import (
    encode_original_source_rows
)

STUDIES=("CZ_Czechia","IT_Italy_Parma","HR_BA_2017_one_study","SI_Slovenia_100m_track")
SITES=(39,28,15,89)
MODELS=("original_source_tip_occupancy_quintile",
        "source_tree_clade_rootdepth_fraction_0.25",
        "source_tree_clade_rootdepth_fraction_0.50",
        "source_tree_clade_rootdepth_fraction_0.75")
FRACTIONS=(None,.25,.5,.75)
N_PER_CHAIN=199
SEEDS=(296119,296120)
MIN_MULTITAXON_CLADES=5
MIN_MOVABLE_INCIDENCE=.25
MIN_SOURCE_GEOS=2
EXPECTED_TIP_COUNTS={
  "LCVP":(433,245,163,195),
  "WP":(416,242,161,189),
  "TPL":(444,248,170,202)}
# From source original 4x3 V.PhyloMaker2 audit run #44517543086.
# These are original CWM height distances ON the DIRECT TIP ONLY
# physical study cohort with >=5 matched taxa per actual site.
EXPECTED_SOURCE_DIRECT_TIP_HEIGHT={
  "LCVP":(.115896747380225,.136848676217088,.1928370852884,.273269737385033),
  "WP":(.111559522478825,.132856099063241,.202736623392164,.251815943995816),
  "TPL":(.117758595809626,.134910314590972,.200609743227737,.251253493136368)
}
TOOL_VERSION="v29.6"
TREE_NO_NEW_TIPS=True

def original_root_depths(edge:np.ndarray,lengths:np.ndarray,n_tips:int,n_nodes:int)->tuple[np.ndarray,np.ndarray,dict]:
    """Find real tree root and compute actual cumulative phylogenetic branch lengths."""
    if edge.ndim!=2 or edge.shape[1]!=2 or len(edge)!=len(lengths):
        raise ValueError("Original source phylo edge and edge.length shapes inconsistent")
    e=np.asarray(edge,dtype=np.int64)
    L=np.asarray(lengths,dtype=np.float64)
    size=n_tips+n_nodes+1
    if e.min()<1 or e.max()>=size or len(e)<n_tips:
        raise ValueError("Original R phylogenetic node ID outside fixed 1-based phylo universe")
    if not np.all(np.isfinite(L)) or np.any(L<0):
        raise ValueError("Original phylogeny has missing, negative or invalid actual branch lengths")
    parent=np.zeros(size,dtype=np.int64)
    children=[[] for _ in range(size)]
    if len(np.unique(e[:,1]))!=len(e):
        raise ValueError("Original source tree has a child with multiple ancestors")
    for i,(p,c) in enumerate(e):
        parent[c]=p
        children[p].append((c,float(L[i])))
    allparents=set(e[:,0].tolist())
    roots=allparents-set(e[:,1].tolist())
    if len(roots)!=1:raise ValueError("Original phylogenetic source root not uniquely determined")
    root=next(iter(roots))
    depth=np.full(size,np.nan,dtype=np.float64)
    depth[root]=0.
    queue=deque([root])
    reached=0
    while queue:
        p=queue.popleft()
        reached+=1
        for child,weight in children[p]:
            if np.isfinite(depth[child]):
                raise ValueError("Original tree cycle or duplicated original child")
            depth[child]=depth[p]+weight
            queue.append(child)
    if reached!=len(allparents|set(e[:,1].tolist())):
        raise ValueError("Disconnected original phylogenetic source tree components")
    dt=depth[1:n_tips+1]
    if not np.all(np.isfinite(dt)) or not np.all(dt>0):
        raise ValueError("Original tree tips have missing/zero root-to-tip path length")
    source={"original_root_node":root,
            "original_total_root_to_tip_distance_min":float(dt.min()),
            "original_total_root_to_tip_distance_median":float(np.median(dt)),
            "original_total_root_to_tip_distance_max":float(dt.max()),
            "root_to_tip_branchlength_CV":float(dt.std(ddof=1)/dt.mean()),
            "all_original_edge_lengths_nonnegative":True,
            "phylogenetic_cut_fraction_of_median_branchlength_not_absolute_ma":True}
    return parent,depth,source

def original_true_clade_groups(tip_positions_one_based:list[int],parent:np.ndarray,
                               depth:np.ndarray,cut:float)->list[np.ndarray]:
    if not math.isfinite(cut) or cut<0 or len(tip_positions_one_based)<2:
        raise ValueError("Original tree clade cut or exact source tip positions invalid")
    by=defaultdict(list)
    for i,tip in enumerate(tip_positions_one_based):
        if tip<1 or tip>=len(parent) or not math.isfinite(depth[tip]):
            raise ValueError("Real phylogenetic tip outside backbone")
        clade=tip
        guard=0
        while depth[clade]>cut:
            nextnode=int(parent[clade])
            if not nextnode:raise ValueError("True source ancestor unavailable at requested branch cut")
            clade=nextnode
            guard+=1
            if guard>len(parent):raise ValueError("Original phylogenetic node path has cycle")
        by[clade].append(i)
    groups=[np.asarray(v,dtype=int) for _,v in sorted(by.items())]
    if sorted(sum((x.tolist() for x in groups),[]))!=list(range(len(tip_positions_one_based))):
        raise ValueError("Original source phylogeny clades do not partition exact tip species")
    return groups

def original_source_private_csv(path:Path)->dict:
    """Input rows are original actual biological source and must remain runner-private."""
    records={}
    with path.open(encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        if reader.fieldnames != ["study","site_num","group_num","taxon",
                                "zlog_canopy_height","zlog_seed_mass"]:
            raise ValueError("Original true botanical source source schema drifted")
        for line in reader:
            study=line["study"]
            if study not in STUDIES:raise ValueError("Not one of 4 original author field publications")
            k=(int(line["site_num"]),line["taxon"])
            if not re.fullmatch(r"[A-Z][A-Za-z-]+ [a-z][a-z-]+",line["taxon"]):
                raise ValueError("Original author botanical source is not strict binomial")
            r=records.setdefault(study,{})
            if k in r:raise ValueError("Original author source duplicate site-taxonomy field record")
            h=float(line["zlog_canopy_height"])
            s=float(line["zlog_seed_mass"])
            if not math.isfinite(h) or not math.isfinite(s):
                raise ValueError("Invalid original matched Oldenburg LE DA height and seed")
            r[k]=(line["group_num"],h,s)
    if set(records)!=set(STUDIES):raise ValueError("Original four authentic publication cohorts missing")
    return records

def source_study_matrix(source_records:dict,study:str,backbone_tip_lookup:dict[str,int])->tuple[np.ndarray,np.ndarray,list[str],np.ndarray,dict]:
    entries=source_records[study]
    nsite=SITES[STUDIES.index(study)]
    source_taxa={t for _,t in entries}
    exact=sorted(source_taxa & set(backbone_tip_lookup))
    if len(exact)<20:raise ValueError("Too few unaugmented original source phylo tips to make ecological comparisons")
    all_sites=defaultdict(dict)
    groups={}
    trait_by_taxon={}
    for (site,tax), (geo,h,s) in entries.items():
        if not 1<=site<=nsite:raise ValueError("Original field study site index drift")
        if site in groups and groups[site]!=geo:raise ValueError("Original site has multiple original geographical fields")
        groups[site]=geo
        if tax in backbone_tip_lookup:
            all_sites[site][tax]=(h,s)
            if tax in trait_by_taxon and trait_by_taxon[tax]!=(h,s):
                raise ValueError("True original species LEDA trait values conflict across original sites")
            trait_by_taxon[tax]=(h,s)
    if set(groups)!=set(range(1,nsite+1)):
        raise ValueError("Original physical sample sites changed")
    tip_counts={site:len(all_sites[site]) for site in groups}
    original_kept_sites=[site for site in range(1,nsite+1) if tip_counts[site]>=5]
    if len(original_kept_sites)<math.ceil(.8*nsite):
        raise ValueError("v29.5 original direct-tip study site support gate failed")
    remaining=sorted({tax for site in original_kept_sites for tax in all_sites[site]})
    idx={t:i for i,t in enumerate(remaining)}
    P=np.zeros((len(original_kept_sites),len(remaining)),dtype=np.int8)
    g=[]
    for row,site in enumerate(original_kept_sites):
        g.append(groups[site])
        for taxon in all_sites[site]:
            P[row,idx[taxon]]=1
    n=P.sum(axis=1)
    if np.any(n<5):raise ValueError("Retained original field plot still has fewer than 5 real backbone tip taxa")
    g=np.asarray(g,dtype=str)
    uniq=sorted(set(g))
    valid_groups=sum(np.sum(g==label)>=2 for label in uniq)
    if valid_groups<2:raise ValueError("Original field geography strata lost under exact-tip matching")
    h=np.asarray([trait_by_taxon[t] for t in remaining],dtype=float)
    pairs=np.triu_indices(len(g),1)
    within=g[pairs[0]]==g[pairs[1]]
    ii=pairs[0][within];jj=pairs[1][within]
    if len(ii)<2:raise ValueError("Original within-publication original geographical site pairs missing")
    W=P/n[:,None]
    original=source_CWM_pairwise_beta(W,h,ii,jj)
    qa={
      "original_study":study,
      "source_original_field_sites":nsite,
      "source_original_species_with_both_LEDA_traits":len(source_taxa),
      "source_original_exact_megatree_tip_species_n":len(exact),
      "source_original_field_sites_with_5_or_more_direct_tips":len(original_kept_sites),
      "source_original_field_site_median_fraction_of_direct_tip_species":
         float(np.median([
            tip_counts[s]/sum(1 for ss,t in entries if ss==s) for s in range(1,nsite+1)])),
      "n_original_matched_taxa_after_site_eligibility":len(remaining),
      "source_original_geographical_groups_with_2_or_more_retained_sites":valid_groups,
      "source_retained_within_original_geography_site_pairs_not_independent_field_replicates":len(ii),
      "source_CWM_height_observed_on_original_direct_tips":original["height"],
      "original_site_taxon_membership_and_LEDA_trait_per_species_unchanged":True,
    }
    if not (0<=qa["source_original_field_site_median_fraction_of_direct_tip_species"]<=1):
        raise ValueError("Original author site denominator error, correction not applied")
    return P,h,remaining,g,{"qa":qa,"W":W,"ii":ii,"jj":jj}

def source_CWM_pairwise_beta(W:np.ndarray,T:np.ndarray,ii:np.ndarray,jj:np.ndarray)->dict:
    if T.ndim!=2 or T.shape[1]!=2 or W.ndim!=2 or W.shape[1]!=T.shape[0]:
        raise ValueError("Original source community×paired traits matrix dimensions invalid")
    C=W@T
    delta=C[ii]-C[jj]
    return {"height":float(np.abs(delta[:,0]).mean()),
            "seed":float(np.abs(delta[:,1]).mean()),
            "joint":float(np.linalg.norm(delta,axis=1).mean())}

def source_occupancy_quintiles(P:np.ndarray)->list[np.ndarray]:
    order=np.argsort(P.sum(axis=0),kind="stable")
    return [a.astype(int) for a in np.array_split(order,5)]
def phylo_exchangeability(P:np.ndarray,groups:list[np.ndarray],g:np.ndarray)->dict:
    mobile=np.unique(np.concatenate([grp for grp in groups if len(grp)>=2])) if any(len(a)>=2 for a in groups) else np.array([],dtype=int)
    if len(mobile)>0:
        fraction=float(P[:,mobile].sum()/P.sum())
    else:fraction=0.
    groupfracs={}
    for geo in sorted(set(g)):
        q=P[g==geo]
        groupfracs[str(geo)]=float(q[:,mobile].sum()/q.sum()) if len(mobile)>0 else 0.
    multi=sum(len(a)>=2 for a in groups)
    group_support=sum(p>=MIN_MOVABLE_INCIDENCE for p in groupfracs.values())
    allowed=multi>=MIN_MULTITAXON_CLADES and fraction>=MIN_MOVABLE_INCIDENCE and group_support>=MIN_SOURCE_GEOS
    return {"source_phylogenetic_clades_in_original_study":len(groups),
            "n_source_clades_with_2plus_real_botanical_species":multi,
            "n_original_trait_known_species_in_exchangeable_clades":int(len(mobile)),
            "source_site_species_incidence_exchangeable_fraction":fraction,
            "source_original_geographic_groups_exchangeable_25percent_n":group_support,
            "source_original_geographical_group_exchangeable_frac":groupfracs,
            "predeclared_real_phylo_exchangeability_gate":allowed}
def permute_paired_traits(T:np.ndarray,groups:list[np.ndarray],rng:np.random.Generator)->np.ndarray:
    out=T.copy()
    for grp in groups:
        if len(grp)>1:out[grp]=T[rng.permutation(grp)]
    if sorted(map(tuple,out))!=sorted(map(tuple,T)):
        raise ValueError("Shuffled bivariate original LE DA trait values not conserved")
    return out
def source_null(P,T,groups,g,ii,jj,study_idx:int,tree_idx:int,mode_idx:int,draws:int)->dict:
    support=phylo_exchangeability(P,groups,g)
    if mode_idx==0:
        support["predeclared_real_phylo_exchangeability_gate"]=True
        status="SOURCE_EXACT_TIP_OCCUPANCY_QUINTILE_REFERENCE"
    else:
        status=("SOURCE_REAL_PHYLOGENETIC_CLADE_DEPTH_ESTIMABLE"
                if support["predeclared_real_phylo_exchangeability_gate"]
                else "NOT_IDENTIFIABLE_PHYLOGENETIC_CLADE_EXCHANGEABILITY")
    result={"status":status,"original_species_site_matrices_FROZEN":True,
            "paired_LEDA_bivariate_vectors_conserved":True,"source_support":support,"metrics":None}
    if status=="NOT_IDENTIFIABLE_PHYLOGENETIC_CLADE_EXCHANGEABILITY":
        return result
    if draws<3:raise ValueError("Insufficient technical label randomizations")
    W=P/P.sum(axis=1)[:,None]
    observed=source_CWM_pairwise_beta(W,T,ii,jj)
    sums={k:[] for k in ("height","seed","joint")}
    for chain,seed in enumerate(SEEDS):
        rng=np.random.default_rng(seed+tree_idx*17011+study_idx*3803+mode_idx*101)
        for rep in range(draws):
            p=permute_paired_traits(T,groups,rng)
            test=source_CWM_pairwise_beta(W,p,ii,jj)
            for key in sums:sums[key].append(test[key])
    result["metrics"]={key:{"observed":observed[key],
            "true_original_clade_label_null_mean":mean(vals),
            "original_observed_minus_original_clade_label_null":observed[key]-mean(vals),
            "technical_upper_rank_NOT_field_p":(1+sum(x>=observed[key] for x in vals))/(len(vals)+1)}
        for key,vals in sums.items()}
    result["n_technical_label_assignments_not_field_studies"]=2*draws
    return result

def run(private:Path,tree_dir:Path,output:Path,draws:int=N_PER_CHAIN)->dict:
    source_records=original_source_private_csv(private)
    alltrees={}
    for tree_ix,name in enumerate(TREES):
        obj,qa=inspect_tree_rdata(tree_dir/("GBOTB.extended."+name+".rda"),name)
        tips=[str(a).replace("_"," ") for a in np.asarray(obj["tip.label"]).reshape(-1)]
        if len(set(tips))!=len(tips):raise ValueError("Original phylogeny normalized source tip names duplicate")
        lookup={name:i+1 for i,name in enumerate(tips)}
        par,depth,tree_depth=original_root_depths(
            np.asarray(obj["edge"]),np.asarray(obj["edge.length"]),
            len(tips),int(np.asarray(obj["Nnode"]).reshape(-1)[0]))
        alltrees[name]={
            "source_authentic_original_tree_quality":qa,
            "source_ancestry_depth":tree_depth,
            "original_studies":{}
        }
        for study_ix,study in enumerate(STUDIES):
            P,T,remaining,geo,r=source_study_matrix(source_records,study,lookup)
            qaq=r["qa"]
            if qaq["source_original_exact_megatree_tip_species_n"]!=EXPECTED_TIP_COUNTS[name][study_ix]:
                raise ValueError("v29.5 true original author×backbone source exact-tip count drift")
            if abs(qaq["source_CWM_height_observed_on_original_direct_tips"]-
                   EXPECTED_SOURCE_DIRECT_TIP_HEIGHT[name][study_ix])>1e-8:
                raise ValueError("v29.5 true original matched height beta drift")
            modes={}
            positions=[lookup[s] for s in remaining]
            median_depth=tree_depth["original_total_root_to_tip_distance_median"]
            for mode_ix,mode in enumerate(MODES):
                frac=FRACTIONS[mode_ix]
                groups=(source_occupancy_quintiles(P) if mode_ix==0 else
                        original_true_clade_groups(positions,par,depth,frac*median_depth))
                calc=source_null(P,T,groups,geo,r["ii"],r["jj"],study_ix,tree_ix,mode_ix,draws)
                calc["original_tree_cut_fraction_of_source_median_root_to_tip_branch_length"]=frac
                calc["original_tree_actual_cut_branchlength"]=frac*median_depth if frac is not None else None
                modes[mode]=calc
            alltrees[name]["original_studies"][study]={"source_original_community_and_true_tip_support_QA":qaq,
                                                        "four_original_phylo_label_nulls":modes}
    summary={}
    for tree,dat in alltrees.items():
        summary[tree]={}
        for mode in MODES:
            supported={}
            for study,record in dat["original_studies"].items():
                d=record["four_original_phylo_label_nulls"][mode]
                if d["metrics"] is not None:
                    supported[study]=d["metrics"]["height"]["original_observed_minus_original_clade_label_null"]
            loo=([mean(v for k,v in supported.items() if k!=removed) for removed in supported]
                 if len(supported)>=3 else [])
            summary[tree][mode]={
                "original_field_publications_n_source_support":len(supported),
                "original_study_level_height_effects":supported,
                "descriptive_original_publication_equal_height_excess":
                   mean(supported.values()) if len(supported)>=3 else None,
                "descriptive_n_original_studies_positive":
                   sum(v>0 for v in supported.values()) if len(supported)>=3 else None,
                "original_publication_LOO_mean_range_not_population_CI":
                   [min(loo),max(loo)] if loo else None,
            }
    result={
        "version":"v29.6",
        "original_phylogeny_publisher_reference":"jinyizju/V.PhyloMaker2 original pinned 2022 commit "+COMMIT,
        "data_provenance":"four genuine original field floristic publication datasets + Oldenburg pinned canopy height and seed mass; precomputed source-author site CSV only in runner temp",
        "original_independent_publisher_studies_n":4,
        "original_study_physical_field_sites_n":171,
        "original_tree_nomenclatures_are_alternative_not_independent_field_papers":list(TREES),
        "source_studies_original_phylo_branch_length_clade_analyses":alltrees,
        "source_12_study_by_tree_original_cohort_ancestry_QA":True,
        "original_four_field_paper_equal_descriptive_summaries":summary,
        "source_branch_cut_fractions_of_root_to_tip_median":[.25,.5,.75],
        "two_technical_label_chains":2,"technical_permutations_per_chain":draws,
        "no_artificial_phylogeny_tip_grafting":True,
        "original_RData_tree_tip_presence_not_equivalent_to_species_having_DNA_sequence":True,
        "original_site_by_species_incidence_and_plant_tips_never_exported":True,
        "source_original_plant_leaves_and_gps_NEVER_released":True,
        "technical_permutation_ranks_NOT_biological_study_p_values":True,
        "no_rail_vs_open_nonrail_matched_causal_field_study":True}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"original_true_phylogenetic_null_status":"SOURCE_AUTHENTICATED","results":{
        tree:{mode:{"studies":x["original_field_publications_n_source_support"],
                     "height_excess":x["descriptive_original_publication_equal_height_excess"],
                     "source_studies_positive":x["descriptive_n_original_studies_positive"]}
              for mode,x in analysis.items()}
        for tree,analysis in summary.items()}},indent=2))
    return result

def _self_test():
    edge=np.array([[5,6],[5,7],[6,1],[6,2],[7,3],[7,4]],dtype=int)
    parent,depth,qa=original_root_depths(edge,np.ones(6),4,3)
    assert qa["original_root_node"]==5
    assert all(abs(depth[1:5]-2)<1e-12)
    groups=original_true_clade_groups([1,2,3,4],parent,depth,1.)
    assert sorted(map(len,groups))==[2,2]
    taxa=np.array([[1.,2.],[3.,4.],[5.,6.],[7.,8.]])
    shuffled=permute_paired_traits(taxa,groups,np.random.default_rng(22))
    assert sorted(map(tuple,shuffled))==sorted(map(tuple,taxa))
    P=np.array([[1,0,1,0],[0,1,1,0],[1,1,0,1]])
    assert sum(map(len,source_occupancy_quintiles(P)))==4
    print("v29.6 ORIGINAL PHYLOGENETIC CLADES Python synthetic root-to-tip, exact phylo clade partition, paired trait test PASS")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--private-original-source",type=Path)
    p.add_argument("--tree-dir",type=Path)
    p.add_argument("--output",type=Path)
    p.add_argument("--draws",type=int,default=N_PER_CHAIN)
    a=p.parse_args()
    if a.self_test:
        _self_test()
    else:
        if not all((a.private_original_source,a.tree_dir,a.output)):
            p.error("Original author source, pinned tree directory and aggregate output required")
        run(a.private_original_source,a.tree_dir,a.output,a.draws)
