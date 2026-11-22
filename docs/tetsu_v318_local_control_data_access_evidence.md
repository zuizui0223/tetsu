# v318 source-access audit: local controls for railway climate-offset hypothesis (2026-11-22)

## Primary estimand
Within-site railway minus **genuine non-railway** vegetation community thermal affinity, and its interaction with macroclimate across independent climatic regions. Within-rail track vs shoulder/toe/embankment contrasts are mechanistic *secondary* analyses only; all ecological indicator values are ordinal preferences, not directly observed thermal microclimate.

## Source gates verified from publisher pages
| Source | Local comparison | Grain / climate replication | Original species-by-plot availability | Status for climate-offset test |
|---|---|---|---|---|
| Winkler et al. (2026), *Land* doi:10.3390/land15040523 | Railway yard, embankment, land vs field margin and adjacent arable land | 4 Czech localities, five habitats, 3 relevés of 10 m² per habitat, 3 surveys in 2021 **within one climatic region** | Publisher explicitly says full data **not publicly available** (supplement S1 taxon names/abbreviations only) | Strong local habitat design; cannot fit without unpublished row-level data; no independent region interaction |
| Dylewski et al. (2024 online; 2025 volume), *Land Degrad Dev* doi:10.1002/ldr.5383 | 25 abandoned railway sites paired with 25 adjacent semi-natural grasslands | Poland, distinct management context (abandoned) | Data available **from corresponding author upon reasonable request** | Best genuine matched natural reference, but not active-rail mechanism or multi-region alone |
| Rakotoson et al. (2025), *Weed Research* doi:10.1111/wre.70025 | 128 track vs 266 adjacent embankment surveys; 2×1 m, 2021, 2 lines around Paris | Rail infrastructure on BOTH sides; one regional climate | Data **on request**; 146/255 trait-characterized species in published analysis | Within-rail source, **not** a nonrail control |
| Šabić et al. (2026), *Plant Biosystems* doi:10.1007/s44473-026-00155-z | 173 same-side shoulder vs toe pairs at 89 field locations in 15 sections; both railway | Slovenia, one-country corridor survey | Original publisher XLSX public, exact SHA256 26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f, archived v102 source audit | v318 actual thermal indicator and missingness reanalysis; cannot demonstrate macroclimatic release |
| Haq et al. (2024), *Environmental Development* doi:10.1016/j.envdev.2024.101064 | Track, embankment, two interior plots at 31 sites, 2014/2017 | Kashmir, 119 km one corridor | Open Dryad doi:10.5061/dryad.0p2ngf20h, source hashes archived | Spatial railway edge effects, **not assured independent nonrail natural reference** |

## Decisions
1. **No valid multi-region railway×nonrail thermal-affinity interaction is yet identified.** Do not substitute national MIREN roadside sites hundreds of kilometres away for a matched local reference.
2. Run the fully traceable original Slovenia *within-rail* thermal indicator paired reanalysis, including >=0, >=0.5, >=0.7, >=0.9 two-habitat species EIV coverage; compute equal-location averages and cluster by original line section. Published inner-vs-outer-zone differences already exist, so this is a methodological robustness result, not a newly discovered general climatic phenomenon.
3. Before requesting restricted original data, specify exact desired columns: anonymized locality, paired plot ID, rail habitat vs nonrail habitat, plot area, original survey year/date, original species identities/presence or cover, plot latitude/longitude or locally associated regional climate, and redistribution restrictions.
4. If paired datasets remain restricted, pivot to a *separate* within-corridor thermal filtering paper rather than naming it a climate-offset discovery. Do not report published species totals or author-only datasets as newly collected tetsu data.

## Source links
- https://www.mdpi.com/2073-445X/15/4/523
- https://onlinelibrary.wiley.com/doi/10.1002/ldr.5383
- https://onlinelibrary.wiley.com/doi/10.1111/wre.70025
- https://link.springer.com/article/10.1007/s44473-026-00155-z
- https://doi.org/10.1016/j.envdev.2024.101064
