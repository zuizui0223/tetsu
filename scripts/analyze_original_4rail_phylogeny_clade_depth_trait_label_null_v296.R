#!/usr/bin/env Rscript
# v29.6 Genuine 4 original railway flora surveys x 3 unaugmented V.PhyloMaker2
# time-scaled backbones. Original species-site incidences, site richness, field
# grouping and real continuous original paired height/seed trait values FROZEN.
# Permute bivariate trait VALUES within true tree clades cut at predeclared
# relative root-to-tip branch-length depths 0.25, 0.50, 0.75; compare an
# occupancy-quintile species trait null IN THE EXACT SAME BACKBONE TIP SUBSET.
# This is an observational spatial trait structure diagnostic, not adaptation.
# All original private taxon lists remain in Actions runner; emit aggregates.

suppressWarnings(options(stringsAsFactors=FALSE))
args <- commandArgs(trailingOnly=TRUE)
SEEDS <- c(296119L,296120L)
DRAWS_EACH <- 199L
CUT_FRACTIONS <- c(0.25,0.5,0.75)
MIN_MOVABLE_CLADES <- 5L
MIN_MOVEABLE_INCIDENCE <- 0.25
MIN_SUPPORTED_SOURCE_GEOGRAPHY <- 2L
STUDIES <- c("CZ_Czechia","IT_Italy_Parma","HR_BA_2017_one_study",
             "SI_Slovenia_100m_track")
BACKBONES <- c("LCVP","WP","TPL")
EXPECTED_TIPS <- matrix(c(
   433L,245L,163L,195L,
   416L,242L,161L,189L,
   444L,248L,170L,202L),
   nrow=4L,ncol=3L,dimnames=list(STUDIES,BACKBONES))
EXPECTED_MATCHED_HEIGHT <- matrix(c(
   0.115896747380225,0.136848676217088,0.1928370852884,0.273269737385033,
   0.111559522478825,0.132856099063241,0.202736623392164,0.251815943995816,
   0.117758595809626,0.134910314590972,0.200609743227737,0.251253493136368),
   nrow=4L,ncol=3L,dimnames=list(STUDIES,BACKBONES))
EXPECTED_N_SITES <- c(39L,28L,15L,89L)
names(EXPECTED_N_SITES) <- STUDIES

tree_depths <- function(edge,edge_length,n_tip,n_node) {
  if (!is.matrix(edge)||ncol(edge)!=2L ||
       length(edge_length)!=nrow(edge) || any(!is.finite(edge_length)) ||
       any(edge_length<0)) stop("Actual original phylo edge/branch lengths invalid")
  nmax <- n_tip+n_node
  if(any(edge>nmax)||any(edge<1L)||anyDuplicated(edge[,2]))
     stop("Phylogeny edge node IDs invalid")
  roots <- setdiff(unique(edge[,1]),edge[,2])
  if(length(roots)!=1L) stop("Phylogeny must have exactly one root")
  parent <- integer(nmax)
  parent[edge[,2]] <- edge[,1]
  distances <- rep(NA_real_,nmax)
  distances[roots] <- 0
  left <- seq_len(nrow(edge))
  iterations <- 0L
  while(length(left)>0L) {
     p <- edge[left,1]
     c <- edge[left,2]
     good <- is.finite(distances[p])
     if(!any(good)) stop("Cannot topologically resolve original root-to-tip branches")
     distances[c[good]] <- distances[p[good]]+edge_length[left[good]]
     left <- left[!good]
     iterations <- iterations+1L
     if(iterations>nmax)stop("Original tree edge ancestry loop")
  }
  if(any(!is.finite(distances[seq_len(n_tip)])) ||
     any(distances[seq_len(n_tip)]<=0)) stop("Some original phylogenetic tips have no measured depth")
  list(parent=parent,root=roots,depth=distances,
       source_tip_age_min=min(distances[seq_len(n_tip)]),
       source_tip_age_median=median(distances[seq_len(n_tip)]),
       source_tip_age_max=max(distances[seq_len(n_tip)]),
       source_tip_root_depth_cv=sd(distances[seq_len(n_tip)])/mean(distances[seq_len(n_tip)]))
}
cut_true_original_tree <- function(matched_tip_node_indices,parent,depth,cut_depth) {
  if(anyNA(matched_tip_node_indices) || any(matched_tip_node_indices<1L) ||
     any(!is.finite(depth[matched_tip_node_indices])) || !is.finite(cut_depth))
       stop("Tree clade source tip support or cutoff invalid")
  target <- as.integer(matched_tip_node_indices)
  tries <- 0L
  while(TRUE) {
    above <- which(depth[target] > cut_depth)
    if(!length(above)) break
    if(any(parent[target[above]]<1L))stop("Phylogenetic clade cut reached above the original root")
    target[above] <- parent[target[above]]
    tries <- tries+1L
    if(tries>length(parent))stop("Phylo cut ancestry failed")
  }
  # As cut grows deeper, the set of disjoint clades refines; every label is
  # a genuine original internal phylogenetic node, not a genus heuristic.
  groups <- split(seq_along(target),as.character(target))
  if(!identical(sort(unlist(groups,use.names=FALSE)),seq_along(target)))
    stop("Original tips not exactly partitioned into phylogenetic clades")
  unname(groups)
}
standard_occupancy_quintile_groups <- function(P) {
  occupancy <- colSums(P)
  order <- order(occupancy,seq_along(occupancy))
  n <- length(order)
  if(n<20L)stop("Original direct-tip source flora has fewer than twenty joint trait taxa")
  out <- split(order,cut(seq_len(n),breaks=5L,labels=FALSE))
  if(length(out)!=5L||any(lengths(out)<2L))stop("Frozen source occupancy quintile support failed")
  unname(out)
}
assess_clade_exchangeability <- function(P,groups,geo,required=TRUE) {
  if(!is.matrix(P)||ncol(P)<20L||nrow(P)<5L ||
       anyNA(P)||!all(P %in% c(0,1)))stop("Invalid true source site-incidence matrix")
  mobile <- sort(unique(unlist(groups[lengths(groups)>=2L],use.names=FALSE)))
  n_multiples <- sum(lengths(groups)>=2L)
  valid <- length(mobile)>0L
  source_total <- sum(P)
  prop <- if(valid)sum(P[,mobile,drop=FALSE])/source_total else 0.
  geo_fractions <- vapply(unique(geo),function(g) {
    sub <- P[geo==g,,drop=FALSE]
    if(nrow(sub)<1L) return(0.)
    if(!length(mobile))return(0.)
    sum(sub[,mobile,drop=FALSE])/sum(sub)
  },numeric(1))
  n_geo <- sum(geo_fractions>=MIN_MOVEABLE_INCIDENCE)
  supported <- valid && n_multiples>=MIN_MOVABLE_CLADES &&
      prop>=MIN_MOVEABLE_INCIDENCE &&
      n_geo>=MIN_SUPPORTED_SOURCE_GEOGRAPHY
  list(n_source_original_phylo_groups=length(groups),
       n_source_phylo_groups_at_least_two_study_taxa=n_multiples,
       n_study_taxa_in_exchangeable_clades=length(mobile),
       source_trait_known_site_taxon_incidence_fraction_exchangeable=prop,
       n_original_geography_groups_with_25pct_exchangeable=n_geo,
       original_geography_fraction_exchangeable=as.list(geo_fractions),
       predeclared_source_min_exchangeable_clades=MIN_MOVABLE_CLADES,
       predeclared_source_min_exchangeable_incidence=.25,
       predeclared_min_supported_geographies=MIN_SUPPORTED_SOURCE_GEOGRAPHY,
       eligible_branch_clade_null=supported,
       status=if(supported)"IDENTIFIABLE_DIRECT_BACKBONE_CLADE_LABEL_NULL"
              else "NOT_IDENTIFIABLE_PHYLOGENETIC_CLADE_EXCHANGEABILITY")
}
matrix_for_original_study <- function(orig, tree_tip_names) {
  if(anyDuplicated(tree_tip_names))stop("Original megatree normalized tips duplicate")
  names <- sort(intersect(unique(orig$taxon),tree_tip_names))
  if(length(names)<20L) stop("Too few exact direct backbone species in study")
  site_ids <- sort(unique(orig$site_num))
  split_sites <- split(orig,factor(orig$site_num,levels=site_ids),drop=TRUE)
  groups <- vapply(split_sites,function(a){
     lab <- unique(a$group_num)
     if(length(lab)!=1L)stop("Site belongs to multiple original publication geography groups")
     as.character(lab)
  },character(1))
  P <- matrix(0L,nrow=length(site_ids),ncol=length(names))
  H <- matrix(NA_real_,nrow=length(names),ncol=2L)
  for(i in seq_along(names)) {
    sub <- orig[orig$taxon==names[[i]],c("zlog_canopy_height","zlog_seed_mass"),drop=FALSE]
    if(any(!is.finite(as.matrix(sub))))stop("Original source height and mass missing")
    if(any(abs(sub$zlog_canopy_height-sub$zlog_canopy_height[[1]])>1e-10) ||
       any(abs(sub$zlog_seed_mass-sub$zlog_seed_mass[[1]])>1e-10))
       stop("Original source taxon had inconsistent 2D LEDA median among sites")
    H[i,] <- as.numeric(sub[1,])
  }
  for(i in seq_along(split_sites)) {
    r <- split_sites[[i]]
    ids <- match(r$taxon,names)
    ids <- ids[!is.na(ids)]
    P[i,ids] <- 1L
  }
  n <- rowSums(P)
  if(any(n<5L))stop("The original author field plot lacks 5 directly placed tip species")
  W <- sweep(P,1,n,"/")
  if(abs(sum(P)-sum(rowSums(P)))>1e-10)stop("Source community matrix drift")
  pair <- which(upper.tri(matrix(0,nrow(P),nrow(P))) &
                outer(groups,groups,"=="),arr.ind=TRUE)
  if(nrow(pair)<1L)stop("No original same-group independent sites for beta")
  list(names=names,P=P,W=W,H=H,geo=groups,pair=pair)
}
trait_cwm_beta <- function(W,H,pairs) {
  cent <- W %*% H
  dis <- cent[pairs[,1],,drop=FALSE]-cent[pairs[,2],,drop=FALSE]
  c(height=mean(abs(dis[,1])),
    seed=mean(abs(dis[,2])),
    joint=mean(sqrt(rowSums(dis*dis))))
}
permute_paired_traits <- function(H,groups) {
  out <- H
  for(g in groups) {
    if(length(g)<2L)next
    donors <- sample(g,length(g),replace=FALSE)
    out[g,] <- H[donors,,drop=FALSE]
  }
  if(!identical(sort(apply(out,1,paste,collapse="~")),
                sort(apply(H,1,paste,collapse="~"))))
     stop("True original 2D trait vector pairs not preserved by phylo relabel")
  out
}
source_null <- function(S,phylo_group,model_id,study_id,mode_idx) {
  original <- trait_cwm_beta(S$W,S$H,S$pair)
  results <- matrix(NA_real_,nrow=2L*DRAWS_EACH,ncol=3L,
                    dimnames=list(NULL,names(original)))
  q <- 1L
  for(k in seq_along(SEEDS)){
    set.seed(SEEDS[[k]]+model_id*17011L+study_id*3803L+mode_idx*101L)
    for(iter in seq_len(DRAWS_EACH)) {
      p <- permute_paired_traits(S$H,phylo_group)
      results[q,] <- trait_cwm_beta(S$W,p,S$pair)
      q <- q+1L
    }
  }
  nmean <- colMeans(results)
  list(
    observed=original,
    null_mean=nmean,
    excess=original-nmean,
    technical_upper_rank_NOT_field_p=vapply(seq_along(original),function(i)
       (1+sum(results[,i]>=original[[i]]))/(nrow(results)+1),numeric(1)),
    total_technical_label_shuffles_NOT_independent_field_studies=nrow(results))
}
self_test <- function(){
  tree <- matrix(c(5,6,5,7,6,1,6,2,7,3,7,4),
                 byrow=TRUE,ncol=2)
  phy <- tree_depths(tree,rep(1,nrow(tree)),4L,3L)
  stopifnot(phy$root==5L,
            all(abs(phy$depth[1:4]-2)<1e-9))
  gs <- cut_true_original_tree(1:4,phy$parent,phy$depth,1)
  stopifnot(length(gs)==2L,identical(sort(lengths(gs)),c(2L,2L)))
  trait <- matrix(c(1,3,4,5,2,6,7,8),ncol=2)
  set.seed(22)
  shuffled <- permute_paired_traits(trait,gs)
  stopifnot(identical(sort(apply(trait,1,paste,collapse="~")),
                       sort(apply(shuffled,1,paste,collapse="~"))))
  p <- matrix(c(1,0,1,0,0,1,0,1,1,1,0,0,0,0,1,1),
              nrow=4L,byrow=TRUE)
  out <- standard_occupancy_quintile_groups(cbind(p,p,p,p,p))
  stopifnot(length(out)==5L,
            sum(lengths(out))==ncol(cbind(p,p,p,p,p)))
  cat("v29.6 synthetic tree-depth, exact clade partition, trait vector and occupancy group self-test PASS\n")
}
if (length(args)==1L && args[[1]]=="--self-test") {
  self_test()
  quit(save="no",status=0)
}
if(length(args)!=4L) stop("Usage: Rscript phylo_v296.R private-original-source.csv tree-dir output.csv sourceQA.json")
self_test()
private_csv <- args[[1]]
tree_folder <- args[[2]]
out_csv <- args[[3]]
out_qa <- args[[4]]
orig <- read.csv(private_csv,stringsAsFactors=FALSE,check.names=FALSE)
if(!identical(names(orig),c("study","site_num","group_num","taxon","zlog_canopy_height","zlog_seed_mass")))
    stop("Original private source field scheme changed")
if(nrow(orig)<2000L || !setequal(unique(orig$study),STUDIES) ||
   length(unique(paste(orig$study,orig$site_num,orig$taxon)))!=nrow(orig))
  stop("Original four independent publication species×site rows invalid")
study_records <- split(orig,factor(orig$study,levels=STUDIES),drop=TRUE)
if(length(study_records)!=4L || any(vapply(study_records,function(d)
   length(unique(d$site_num)),integer(1))!=EXPECTED_N_SITES))
   stop("Authentic original 4-source geographic site units drift")
out <- list()
tree_qa <- list()
for(ti in seq_along(BACKBONES)) {
  tree_name <- BACKBONES[[ti]]
  original_file <- file.path(tree_folder,paste0("GBOTB.extended.",tree_name,".rda"))
  if(!file.exists(original_file))stop("Source megatree not downloaded")
  env <- new.env(parent=emptyenv())
  rdata <- load(original_file,envir=env)
  obj_name <- paste0("GBOTB.extended.",tree_name)
  if(!obj_name %in% rdata)stop("Author megatree original R object name drift")
  tree <- env[[obj_name]]
  tips <- gsub("_"," ",tree$tip.label,fixed=TRUE)
  if(anyDuplicated(tips))stop("Original source tree has duplicate tip names")
  node <- tree_depths(tree$edge,tree$edge.length,length(tips),tree$Nnode)
  # Cutoffs relative to the source author's median total root-to-tip branch
  # length; NEVER label them as absolute divergence times in Ma.
  scale_depth <- node$source_tip_age_median
  tree_qa[[tree_name]] <- list(
    n_source_megatree_tips=length(tips),
    source_tree_edge_count=nrow(tree$edge),
    root_to_tip_length_min=node$source_tip_age_min,
    root_to_tip_length_median=node$source_tip_age_median,
    root_to_tip_length_max=node$source_tip_age_max,
    root_to_tip_length_cv=node$source_tip_root_depth_cv,
    cutoff_basis="fraction of original published tree median root-to-tip total branch length, NOT claim of absolute Ma",
    no_species_added_to_source_tree=TRUE)
  for(si in seq_along(STUDIES)) {
    study <- STUDIES[[si]]
    rec <- study_records[[study]]
    S <- matrix_for_original_study(rec,tips)
    if(length(S$names)!=EXPECTED_TIPS[study,tree_name])
       stop("v29.5 true source direct-tip count drift: perhaps taxonomy or author data changed")
    original_raw <- trait_cwm_beta(S$W,S$H,S$pair)
    if(abs(original_raw["height"]-EXPECTED_MATCHED_HEIGHT[study,tree_name])>1e-11)
       stop("v29.5 source exact original true direct-tip stature beta drift")
    exact_tip_idx <- match(S$names,tips)
    if(anyNA(exact_tip_idx))stop("Direct tip identity correspondence changed")
    base_groups <- standard_occupancy_quintile_groups(S$P)
    phylo_modes <- list(base_groups)
    labels <- "original_species_occupancy_quintile"
    depths <- NA_real_
    fractions <- NA_real_
    for(fraction in CUT_FRACTIONS) {
      real_clades <- cut_true_original_tree(exact_tip_idx,node$parent,node$depth,
                                           scale_depth*fraction)
      phylo_modes[[length(phylo_modes)+1L]] <- real_clades
      labels <- c(labels,paste0("true_tree_clade_depth_fraction_",sprintf("%.2f",fraction)))
      fractions <- c(fractions,fraction)
      depths <- c(depths,scale_depth*fraction)
    }
    for(mi in seq_along(phylo_modes)) {
      m <- phylo_modes[[mi]]
      support <- assess_clade_exchangeability(S$P,m,S$geo)
      if(mi==1L) {
        # A large occupancy-quintile null is always source-eligible, but not a
        # tree-clade permutation; label scope honestly.
        support$eligible_branch_clade_null <- TRUE
        support$status <- "ORIGINAL_TIP_SUBSET_OCCUPANCY_NULL_REFERENCE"
      }
      values <- if(support$eligible_branch_clade_null)
                  source_null(S,m,ti,si,mi)
                else NULL
      out[[length(out)+1L]] <- data.frame(
        backbone=tree_name,
        study=study,
        null_model=labels[[mi]],
        relative_cut_fraction=fractions[[mi]],
        original_tree_branch_length_depth_cut=depths[[mi]],
        n_original_physical_field_sites=nrow(S$P),
        n_original_trait_known_exact_tree_tips=ncol(S$P),
        n_orig_source_study_plant_geographical_groups=length(unique(S$geo)),
        source_original_within_geography_site_pairs=nrow(S$pair),
        n_true_tree_clusters_at_cut=support$n_source_original_phylo_groups,
        n_true_exchangeable_clusters=support$n_source_phylo_groups_at_least_two_study_taxa,
        n_species_in_exchangeable_source_clades=support$n_study_taxa_in_exchangeable_clades,
        fraction_true_site_taxon_incidence_exchangeable=support$source_trait_known_site_taxon_incidence_fraction_exchangeable,
        n_source_geography_groups_with_25pct_movable=support$n_original_geography_groups_with_25pct_exchangeable,
        prespecified_support_status=support$status,
        n_technical_trait_label_shuffles_if_identifiable=
          if(is.null(values))0L else 2L*DRAWS_EACH,
        original_observed_exact_backbone_CWM_height=original_raw[["height"]],
        source_height_null_mean=if(is.null(values))NA_real_ else values$null_mean[["height"]],
        observed_minus_clade_null_height=if(is.null(values))NA_real_ else values$excess[["height"]],
        source_height_technical_upper_rank_NOT_field_p=if(is.null(values))NA_real_ else values$technical_upper_rank_NOT_field_p[[1]],
        original_observed_exact_backbone_CWM_seed=original_raw[["seed"]],
        observed_minus_clade_null_seed=if(is.null(values))NA_real_ else values$excess[["seed"]],
        original_observed_exact_backbone_CWM_joint2D=original_raw[["joint"]],
        observed_minus_clade_null_joint2D=if(is.null(values))NA_real_ else values$excess[["joint"]],
        original_site_species_geography_FIXED=TRUE,
        2D_LEDA_height_and_seed_mass_paired_vector_preserved=TRUE,
        tree_tip_NOT_guaranteed_to_have_own_sequenced_DNA=TRUE,
        stringsAsFactors=FALSE)
    }
  }
}
summary <- do.call(rbind,out)
if(nrow(summary)!=length(BACKBONES)*length(STUDIES)*(1L+length(CUT_FRACTIONS)))
  stop("Four genuine railway study field units x three exact-tip trees x four source nulls incomplete")
dir.create(dirname(out_csv),recursive=TRUE,showWarnings=FALSE)
write.csv(summary,out_csv,row.names=FALSE,na="")
library_meta <- list(
 version="v29.6",
 scope="4 independent European railway floristic source papers x 3 author original tree nomenclatures x 3 fixed branch-length clade cuts; not 12 independent studies",
 original_VPhyloMaker2_github_commit="7af3fb5152f691af2e4ec9d5e2e467d1b50505e9",
 DRAWS_each_source_chain=DRAWS_EACH,
 source_independent_technical_label_chains=length(SEEDS),
 source_original_plant_site_rows_PRIVATE_in_runner=TRUE,
 source_true_tree_properties=tree_qa,
 original_inclusion_policy="Only species already on original backbone with original LE DA height and seed mass; no graft or aliases",
 source_phylo_cut_fractions=as.list(CUT_FRACTIONS),
 predeclared_exchangeability_minimums=list(
   source_multispecies_tree_clades=MIN_MOVABLE_CLADES,
   source_phylo_matched_site_taxon_incidence=MIN_MOVEABLE_INCIDENCE,
   source_original_geography_groups=MIN_SUPPORTED_SOURCE_GEOGRAPHY),
 n_original_source_studies=length(STUDIES),
 n_backbones=length(BACKBONES),
 source_outcome_estimable_counts_by_tree_mode=as.list(
   table(summary$backbone,summary$null_model,
        summary$prespecified_support_status=="IDENTIFIABLE_DIRECT_BACKBONE_CLADE_LABEL_NULL")),
 source_original_tip_names_and_field_locations_NOT_exported=TRUE,
 not_independent_field_studies_or_phylogenetic_meta_p_values=TRUE,
 branch_length_tree_backbone_does_not_prove_every_tip_was_originally_DNA_sequenced=TRUE)
# JSON package may not be installed in R default: output concise DCF/CSV metadata.
# main Python action will assemble provenance JSON without raw source plant rows.
writeLines(c(
 "version: v29.6",
 paste("n_original_publications:",length(STUDIES)),
 paste("n_independent_backbone_nomenclatures:",length(BACKBONES)),
 paste("source_tree_root_to_tip_LCVP:",tree_qa$LCVP$root_to_tip_length_median),
 paste("source_tree_root_to_tip_WP:",tree_qa$WP$root_to_tip_length_median),
 paste("source_tree_root_to_tip_TPL:",tree_qa$TPL$root_to_tip_length_median),
 paste("source_tree_depth_CV_LCVP:",tree_qa$LCVP$root_to_tip_length_cv),
 paste("source_tree_depth_CV_WP:",tree_qa$WP$root_to_tip_length_cv),
 paste("source_tree_depth_CV_TPL:",tree_qa$TPL$root_to_tip_length_cv),
 "original_phylo_backbone_tip_taxa_graft_count: 0",
 "n_source_null_modes_per_tree_and_study: 4",
 "original_raw_species_site_geo_info_in_public_files: false"),
 out_qa)
cat("v29.6 SOURCE authenticated clade null summary, all 3 real mega-trees:\n")
for(tree in BACKBONES) {
  z <- summary[summary$backbone==tree,]
  for(mode in unique(z$null_model)) {
    sub <- z[z$null_model==mode,]
    passed <- sum(!is.na(sub$observed_minus_clade_null_height))
    v <- sub$observed_minus_clade_null_height
    cat(tree,mode, " estimable=",passed,"/",length(v),
        " source study positive=",sum(v>0,na.rm=TRUE),"\n",sep="")
  }
}
