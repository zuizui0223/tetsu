# v345 frozen after v344: strict HABITAT versus broad context on the SAME source papers / blocks

## Exact unresolved confound
v343 broad historical RAIL/ROAD source context mean +0.049381 across 32 original publications; v344 strict HABITAT-only is −0.052983 across only 5 source publications / 6 blocks. Because source paper and vegetation samples change with label definition, **a difference in means is NOT an effect of label reliability**. The v345 job will distinguish source-paper selection from original RAIL/ROAD context definition.

## Predeclared analysis, prior to v345 results
- **Repeat identical authoritative CzechVeg-Open v1.1 original source** (Zenodo 19231152, 117739 full plot headers MD5 `0bc920462d1f9eac473e7b2a231e6782`, 2510340 plot×plant records MD5 `b296e5e19f03fb3bad0a8347c317f8a3`) and EIVE1.0 original source MD5 `0f0ce19cd2a781eeed027d8b0acadbcb`. Same original S2 `REFERENCE_TITLE × phytoregion × decade × area` blocks, at least two RAIL and two ROAD source plots per original block.
- Primary EIVE support >=3 independent source systems, original plot richness>=5 and EIVE-known at least 3/50% plot species. >=70% secondary. Original broad and strict blocks formed from source text alone **before thermal indicators**.
- Define `strict_papers` as source publication titles in post-EIVE strict HABITAT-only matched blocks (v344 5). Never add/drop papers based on temperature effect.
- **Contrast A**: strict HABITAT matched block outcomes in exactly the strict_papers.
- **Contrast B**: broad HABITAT+LOCALITY matched blocks, restricted to exactly the same strict_papers. Compare strict and broad mean only on `papers_available_in_both`; explicitly report source study counts and number of source blocks in each.
- **Contrast C**: same source S2 block keys where strict and broad source plot support both pass trait gate (intersection). Compare strict/broad outcomes in identical original publication×region×decade×area blocks, **but source plots still differ**, and report how many survive.
- **Contrast D (identity check)**: use only strict source plot IDs that ALSO carry the same broad RAIL/ROAD label, in the same shared source S2 block; if each original habitat arm has >=2, compute identical-plot CWM and check strict/broad exactly equal. Only used to verify that any observed difference between methods comes from a changed source plot cohort, not an estimator implementation bug.
- Original publication equal weighting, no paper×block pseudoreplication; any N <5 papers or <10 shared blocks = **descriptive audit only**. Report all support and signs, no field-randomized p. Never interpret joint paper effects as a causal result of labelling.

### Stopping decision
- If broad within the strict paper cohort is already negative (or not positive), broad 32-paper +0.049 and strict 5-paper −0.053 likely reflect sample selection/physical study differences; this does not show mislabeling creates spurious warming.
- If broad within the exact same papers/blocks stays clearly positive but strict becomes negative, exact source labelling/sampling definition is influential, still NOT proven causal: field plot selections differ.
- If overlapping study blocks too few, declare identifiable effect **unresolved**. In no case label CzechVeg historical RAIL text as true current ballast, train frequency, microclimate or railway-vs-natural matched plot data.
- Supplementary, post hoc response after observing v343/v344; not confirmatory global-effect p.

Only report aggregate study n and contrasts; do not re-publish raw author plot×species data, IDs, source titles or GPS. Stack in existing v342–v345 Draft PR #54, main remains untouched.
