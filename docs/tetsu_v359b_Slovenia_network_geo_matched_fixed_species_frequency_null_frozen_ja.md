# v359b — actual rail+geographical distance matched Slovenian flora vs fixed source-species-frequency null

Frozen 2026-11-24 AFTER real v359 original botanical data Jaccard +0.0270333 was observed. **Post hoc source null sensitivity**, not prospective independent prediction.

## Fixed data and geometric triads
- Source Slovenian 2026 genuine publisher original 401,279-byte XLSX, SHA256 `26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f`, 89 original sites/15 genuine author source line sections, 2022–2023 physical ballast BS+BT flora.
- Pinned Overture 2026-11-18.0 true rail connectivity; v358 source site matching 86/89 + 3 physically ambiguous; v359 exactly 27 source anchor triads / 7 publisher sections where actual source geodesic distance difference <=2.5 km, original true network path difference <=5 km, each physical source distance <=40 km. Reference sites selected strictly from source geographical and physical rail network distances, before reading plant similarity; **same triads** as v359 real verified result, no retuning.
- Original Jaccard source anchor SAME author rail section minus DIFFERENT original author rail section = **+0.027033311999075067** equal-weighted over seven source sections. Reject if source values or cohort differ. No use of actual train OD or same living individual.

## Exact null models
Use original published 89 physical-site union BS+BT species sets, don't subset to only the 27 anchors to define source species prevalence (avoids conditioning original pool on positives):
- H0-A: 2×2 checkerboard swaps of species detections between site rows, conserving **each site's actual taxon richness** and **each species' exact total number of observed original source sites** across all 89.
- H0-B, stronger source-survey control: same checkerboard swaps restricted within original source year (2022 or 2023) × original `SINGLE/DOUBLE` layout (same groups as prior v28.1). Preserve each site's total richness and species prevalence *within each year×track-layout stratum*.
- Run **2 seeded correlated MCMC chains×99 retained draws**, fixed 3,000 burn-in proposals and 350 proposals between draws per source group; preserve source guild-free species identity taxonomic Jaccard on exact same original network+geo matched triads. Report original observed gap, null means, observed minus null, per original author section observed sign and null source section means, technical rank (NOT randomized field `p`) and MCMC chain diagnostics. Never treat 198 null values as independent field studies.
- Scientific stopping: if the true source +0.027 effect is not exceptional against within-year×layout occupancy-preserving null, **network-specific seed transfer is not supported**. Even an exceptional effect is only observational and could reflect environmental or regional species pools and original line characteristics, not randomized train movement.
- No extra tuned null grouping or number of draws chosen based on outcome. Null was frozen before this new run but after seeing v359 observed +0.027, so claim exploratory robustness only.

Output only original paper SHA and **aggregate** line-equal null statistics; no original taxa, GPS, paper station/line IDs or raw Overture geometry exported. Keep source Draft PR #58, main unchanged. Script `scripts/analyze_Slovenia_networkmatched_species_frequency_null_v359b.py` with original source CI.
