# tetsu v28.5–28.7: Railway annual-filter global life-cycle source audit

Date: 2026-11-18 JST. Research branch only; main unchanged.

## Authenticated field evidence, distinct statistical targets
- Slovenia: Šabić et al. (2026) DOI https://doi.org/10.1007/s44473-026-00155-z, original 89 railway field sites/173 side-pairs and 15 author railway sections. Direct BS−BT strict annual species-fraction = +0.054703. Original-section cluster bootstrap 2.5–97.5%: [+0.017619,+0.086547]; sign invariant across 15 leave-one-section-out runs (see docs/Slovenia_rail_section_cluster_and_missing_trait_sensitivity_v105.md). Section bootstrap is not all-world randomized field confidence.
- France: Rakotoson et al. (2025) DOI https://doi.org/10.1111/wre.70025. Original 128 track/266 embankment vegetation quadrats; published annual/ruderal increase on track, source plot matrix and variance not freely published. Directional original-paper replication only; do not calculate faux standardized effect or SE.
- Kashmir India: Rashid et al. (2021) DOI https://doi.org/10.1111/1365-2664.13961, genuine Dryad version 99735 https://doi.org/10.5061/dryad.0p2ngf20h, original 31 sites, 4 distinct 100m2 zones, two years, 123 taxon columns; source SHA256 529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474. Full 4-zone 2014/2017 common cohort = 27 sites (15 station, 12 between). Original author's life-form Table3 excludes Track. Raw direct annual-fraction delta still requires taxon-level life span.

## Two independent potential trait sources
- Local alpine Kashmiri checklist: Rasray & Lone (2025) https://doi.org/10.5061/dryad.vmcvdnd4c, published source description lists 159 taxa (20 annual, 139 perennial) and a Life Span column. v28.6 workflow run https://github.com/zuizui0223/tetsu/actions/runs/44293386827 FAILED when the second original archive download returned HTTP 401; local XLSX not authenticated; matching taxa and effects unknown. Do not call this an observed replication.
- Global life-cycle source: Poppenwimer, Mayrose & DeMalach (2023) Nature DOI https://doi.org/10.1038/s41586-023-06644-x. 235,979 consensus species from 11 original databases, available via Figshare collection DOI https://doi.org/10.6084/m9.figshare.c.6176239; individual LifeForm Trait Database item 20764777. This is external trait information, NOT an independent railway vegetation study.
- v28.7 original Figshare CSV was successfully authenticated in Actions run https://github.com/zuizui0223/tetsu/actions/runs/44293988772: original filename LifeForm_Trait_Database.csv, source SHA256 83b23e7874d73387a0c24aeae3e2931c50ec7c78cad52ff26cdfc99d1b6623b1; columns Family, Genus, Accepted Name, LifeCycle, Woodiness; however initial crosswalk yielded 0/123 known traits (0/27 eligible sites). This is currently treated as a **source-name interpretation/normalization audit**, NOT real absence of matching species. A second diagnostic run is in progress.

## Frozen scientific estimands
- Primary source site effect: equal-original-site mean of per-site (2014/2017 mean of Track strict annual species fraction minus Interior Plot1 strict annual fraction), only sites where both source years and all four zones have genuine data. Trait-known incidence >=70% in both contrast zones in both years; 50% exploratory threshold.
- Secondary zones: Track−Interior2 and Track−Embankment, Embankment−Interior1/2, station and between-station spatial strata. Station is not observed frequency of freight trains.
- Life span mapping must be exact accepted species or a conservative, unambiguous same-species binomial; ambiguous, missing and conflicting values withheld. Genus-level life form inference prohibited.
- No gene, plasticity, freight delivery or railway-caused causal adaptation claims can be inferred from these compositional field results.
- The Slovenia BS/BT gradient and India Track/Interior gradient are physically distinct exposures and cannot be pooled as if they were one homogeneous effect. Meta-analysis predictive intervals and leave-continent-out external skill remain unestimated until true effect and SE from >=2 independently comparable studies. Seven or eight studies across three continents remain a recruitment target rather than a completed result.

## Reproducibility
- v28.5 six synthetic tests and real workflow SUCCESS: https://github.com/zuizui0223/tetsu/actions/runs/44292408854
- v28.6 source access failure but twelve synthetic checks pass: https://github.com/zuizui0223/tetsu/actions/runs/44293386827
- v28.7 strict source verification success but 0-match debug gate: https://github.com/zuizui0223/tetsu/actions/runs/44293988772
- Main script: scripts/analyze_Kashmir_global_lifecycle_source_v287.py; tests/test_Kashmir_global_lifecycle_source_v287.py; source workflow .github/workflows/tetsu-v287-global-lifecycle-Kashmir.yml. Only derivative crosswalk and aggregate JSON are written; no original site x species matrix redistributed.
