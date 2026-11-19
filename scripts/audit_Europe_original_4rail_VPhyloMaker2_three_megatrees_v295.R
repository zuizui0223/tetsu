#!/usr/bin/env Rscript
# v29.5: read three real V.PhyloMaker2 backbone phylo R-data objects.
# INPUT: ORIGINAL SOURCE plant species, plot and trait rows on CI RUNNER ONLY.
# OUTPUT: study-level coverage, study-level CWM sensitivity and tree geometry.
# Do NOT export species identities, phylo original tips or original plot rows.
# No tip-grafting is allowed; these are EXACT-BACKBONE-PLACEMENT coverage stats.

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript audit_...v295.R original_source_CI_ONLY.csv tree_dir aggregate.csv")
}
private_source <- args[[1]]
tree_dir <- args[[2]]
output_path <- args[[3]]
if (!file.exists(private_source)) stop("Original author source transient CSV missing")
df <- read.csv(private_source, stringsAsFactors=FALSE, check.names=FALSE,
               colClasses=c(study="character",site_num="integer",
                            group_num="character",taxon="character",
                            zlog_canopy_height="numeric",zlog_seed_mass="numeric"))
required <- c("study","site_num","group_num","taxon",
              "zlog_canopy_height","zlog_seed_mass")
if (!identical(names(df),required)) stop("Original source header/column drift")
studies <- c("CZ_Czechia","IT_Italy_Parma","HR_BA_2017_one_study",
             "SI_Slovenia_100m_track")
knownsites <- c(39L,28L,15L,89L)
names(knownsites) <- studies
known_height_beta <- c(
  CZ_Czechia=0.11789869761794586,
  IT_Italy_Parma=0.13961365832072362,
  HR_BA_2017_one_study=0.18874831982750337,
  SI_Slovenia_100m_track=0.24102203205485576
)
if (nrow(df)<2000L || length(unique(df$study))!=4L ||
    !setequal(unique(df$study), studies)) stop("Source-authentic field studies drift")
if (any(!is.finite(df$zlog_canopy_height)) ||
    any(!is.finite(df$zlog_seed_mass)) ||
    any(df$site_num < 1L) || any(is.na(df$taxon)) ||
    any(!grepl("^[A-Z][A-Za-z-]+ [a-z][a-z-]+$", df$taxon))) {
  stop("Source contains invalid original plant taxon or scientific trait")
}
k <- paste(df$study,df$site_num,df$taxon,sep=":::")
if (anyDuplicated(k)) stop("Original site-by-species uniqueness lost")

original_study <- split(df, factor(df$study,levels=studies),drop=TRUE)
for (study in studies) {
  d <- original_study[[study]]
  if (length(unique(d$site_num))!=knownsites[[study]])
    stop("Original railway study true field-site count changed")
  if (min(table(d$site_num))<5L) stop("Original field plot fails >=5 true joint trait species")
  g <- tapply(d$group_num,d$site_num,unique)
  if (any(lengths(g)!=1L)) stop("One original plot belongs to multiple published source geographic groups")
}

cmw_beta <- function(heights,groups,eligible=rep(TRUE,length(heights))) {
  ii <- which(eligible & is.finite(heights))
  if(length(ii)<2L)return(NA_real_)
  h <- heights[ii]
  g <- groups[ii]
  mask <- upper.tri(outer(h,h,"-")) & outer(g,g,"==")
  if (!any(mask))return(NA_real_)
  mean(abs(outer(h,h,"-")[mask]))
}
measured_species <- lapply(original_study,function(d)unique(d$taxon))
models <- c("LCVP","WP","TPL")
expected_tips <- c(LCVP=73420L,WP=72570L,TPL=74529L)
output <- list()
for (model in models) {
  file <- file.path(tree_dir,paste0("GBOTB.extended.",model,".rda"))
  if (!file.exists(file)) stop(paste("Missing author ORIGINAL backbone",model))
  env <- new.env(parent=emptyenv())
  objects <- load(file,envir=env)
  oname <- paste0("GBOTB.extended.",model)
  if(!oname %in% objects)stop("Author megatree R-data object name drift")
  tree <- env[[oname]]
  if(is.null(tree$tip.label)||is.null(tree$edge.length)||
     is.null(tree$edge)||is.null(tree$Nnode))stop("Original backbone lacks phylo or branches")
  n_tip <- length(tree$tip.label)
  if (abs(n_tip-expected_tips[[model]])>50L ||
      nrow(tree$edge)<n_tip || length(tree$edge.length)!=nrow(tree$edge)) {
    stop("Original backbone tip/edge support does not match published mega-tree")
  }
  if (any(!is.finite(tree$edge.length)) ||
      any(tree$edge.length<0)) stop("Original backbone branch lengths absent or invalid")
  tree_tips <- gsub("_"," ",tree$tip.label,fixed=TRUE)
  if(anyDuplicated(tree_tips))stop("Original mega-tree nonunique tip labels after underscore conversion")
  n_positive <- sum(tree$edge.length>0)
  n_zero <- sum(tree$edge.length==0)
  # A backbone tip denotes its inclusion in a constructed mega-tree,
  # NOT a guarantee of direct molecular sequencing for that taxon.
  terminal_edges <- match(seq_len(n_tip),tree$edge[,2])
  if(anyNA(terminal_edges))stop("Not every original mega-tree tip has a terminal edge")
  terminal_len <- tree$edge.length[terminal_edges]
  n_taxa_global <- length(tree_tips)
  tipset <- new.env(parent=emptyenv(),hash=TRUE,size=round(n_tip*1.3))
  for (txt in tree_tips) assign(txt,TRUE,envir=tipset)
  tree_genera <- unique(sub(" .*","",tree_tips))
  genusset <- new.env(parent=emptyenv(),hash=TRUE,size=round(length(tree_genera)*1.3))
  for (txt in tree_genera) assign(txt,TRUE,envir=genusset)
  for (study in studies) {
    d <- original_study[[study]]
    is_exact <- vapply(d$taxon,exists,logical(1),envir=tipset,inherits=FALSE)
    total_spec <- unique(d$taxon)
    exact_species <- unique(d$taxon[is_exact])
    has_backbone_genus <- vapply(
        sub(" .*","",d$taxon),exists,logical(1),envir=genusset,inherits=FALSE)
    site_chunks <- split(d, d$site_num)
    groups <- vapply(site_chunks,function(x)unique(x$group_num),character(1))
    original_site_heights <- vapply(site_chunks,function(x)mean(x$zlog_canopy_height),numeric(1))
    matches_per_site <- vapply(site_chunks,
      function(x)sum(vapply(x$taxon,exists,logical(1),envir=tipset,inherits=FALSE)),integer(1))
    frac_per_site <- matches_per_site / lengths(site_chunks)
    all_raw_height <- cmw_beta(original_site_heights,groups)
    if(!is.finite(all_raw_height) ||
       abs(all_raw_height-known_height_beta[[study]])>1e-10) {
      stop("Original 4-study source height CWM baseline drifted before phylogeny")
    }
    site_matched_h <- vapply(site_chunks,function(x){
      ok <- vapply(x$taxon,exists,logical(1),envir=tipset,inherits=FALSE)
      if (!any(ok))return(NA_real_)
      mean(x$zlog_canopy_height[ok])
    },numeric(1))
    strict_sites <- matches_per_site>=5L
    valid_geos <- tapply(strict_sites,groups,sum)
    geos_with_2 <- sum(valid_geos>=2L)
    strict_sample <- sum(strict_sites)
    # Source-prespecified "DIRECT TIP" inference gate.
    # No grafted/ghost branches under this gate.
    strict_candidate <- length(exact_species)/length(total_spec)>=.6 &&
        mean(is_exact)>=.6 && mean(strict_sites)>=.8 && geos_with_2>=2L
    eligible_study_comparator <- cmw_beta(
        original_site_heights,groups,eligible=strict_sites)
    matching_study_comparator <- cmw_beta(
        site_matched_h,groups,eligible=strict_sites)
    exact_idx <- match(exact_species,tree_tips)
    if (anyNA(exact_idx)) stop("Tree tip lookup gave index mismatch")
    matched_terminal <- terminal_len[exact_idx]
    output[[length(output)+1L]] <- data.frame(
      phylogeny=model,
      original_study=study,
      backbone_source_tip_n=n_tip,
      backbone_source_edge_n=nrow(tree$edge),
      backbone_source_positive_edge_fraction=n_positive/length(tree$edge.length),
      backbone_source_zero_edge_n=n_zero,
      source_original_station_or_ballast_site_n=length(site_chunks),
      source_actual_original_geographical_group_n=length(unique(groups)),
      source_taxon_dual_LEDA_n=length(total_spec),
      exact_original_backbone_tip_taxon_n=length(exact_species),
      exact_original_backbone_tip_taxon_fraction=length(exact_species)/length(total_spec),
      source_site_by_taxon_incidence_n=nrow(d),
      exact_backbone_tip_incidence_n=sum(is_exact),
      exact_backbone_tip_incidence_fraction=mean(is_exact),
      original_taxa_backbone_genus_incidence_fraction=mean(has_backbone_genus),
      direct_tip_species_with_positive_terminal_edge=sum(matched_terminal>0),
      direct_tip_species_with_zero_terminal_edge=sum(matched_terminal==0),
      original_field_site_ge5_direct_tip_n=strict_sample,
      original_field_site_ge5_direct_tip_fraction=mean(strict_sites),
      original_field_site_ge50pct_tip_fraction_n=sum(frac_per_site>=.5),
      original_field_site_tip_match_fraction_median=median(frac_per_site),
      original_field_site_matched_species_min=min(matches_per_site),
      original_field_site_matched_species_median=median(matches_per_site),
      original_field_site_matched_species_max=max(matches_per_site),
      original_field_geography_groups_ge2_supported_sites_n=geos_with_2,
      source_original_full_cohort_CWM_height_between_site_distance=all_raw_height,
      same_study_retained_sites_original_taxa_CWM_height_between_site_distance=eligible_study_comparator,
      same_study_retained_sites_direct_tip_taxa_CWM_height_between_site_distance=matching_study_comparator,
      prespecified_exact_tip_no_grafting_inference_gate_status=if(strict_candidate)
            "DIRECT_BACKBONE_TIP_SOURCE_SUPPORT_PASS" else "DIRECT_BACKBONE_TIP_SOURCE_SUPPORT_FAIL",
      original_sites_and_species_LEDA_values_UNCHANGED=TRUE,
      no_augmented_tip_is_claimed_as_observed_branch=TRUE,
      no_claim_of_molecular_sequencing_for_backbone_tip=TRUE,
      stringsAsFactors=FALSE
    )
  }
}
out <- do.call(rbind, output)
if(nrow(out)!=12L ||
   any(table(out$phylogeny)!=4L)||
   any(out$source_original_station_or_ballast_site_n<10L)||
   any(is.na(out$exact_original_backbone_tip_taxon_fraction))) {
  stop("Three original mega-trees vs four independent rail studies yielded incomplete QA")
}
dir.create(dirname(output_path),recursive=TRUE,showWarnings=FALSE)
write.csv(out,output_path,row.names=FALSE,na="")
# Do not print original species names or raw plant observations.
cat("v29.5 verified original source: ",
    length(models),"fixed-backbone nomenclatures,",
    length(studies),"true original railway field publications,",
    nrow(out),"source QA rows; no original taxon identities released\n",sep=" ")
for (tree in models) {
  z <- out[out$phylogeny==tree,]
  cat(tree,": original direct-tip support pass=",
      sum(z$prespecified_exact_tip_no_grafting_inference_gate_status==
          "DIRECT_BACKBONE_TIP_SOURCE_SUPPORT_PASS"),"/4,",
      " source species tip fraction range ",
      sprintf("%.3f",min(z$exact_original_backbone_tip_taxon_fraction)),"-",
      sprintf("%.3f",max(z$exact_original_backbone_tip_taxon_fraction)),"\n",sep="")
}
