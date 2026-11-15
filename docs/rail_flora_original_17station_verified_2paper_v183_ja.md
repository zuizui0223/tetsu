# tetsu v18.1–v18.3: 2 original railway-station flora matrices with source-locked text/glyph QA

Date: 2026-11-15. Branch: `research/railway-seed-dispersal-v01`; PR #1 DRAFT. `main` unchanged. This is a **documented SOURCE-verified result**, **NOT** a regional-pool-conditioned multi-country functional-homogenization result.

## Verified sources and outputs

| Original independent paper | Spatial units | Taxa and timing | Source validation | Output |
| --- | ---: | --- | --- | --- |
| Májeková et al., Thaiszia 2016, East Slovakia (publisher 16-page PDF) | 2 stations with full Tab 2; 4 station×time strata | Košice contemporary 192 species; Spišská Nová Ves 191 species (published author counts); historical sources incomplete | Source PDF SHA256 `cd87a16377ef351b6aa1f51a3c9302f060f0f52431c5f8cac14ecfbbb94b19eb`; original raw `Setaria viridis` subspecies retained separately and species collapse reconciles source totals | [v18.1 SOURCE SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/43722187035), original station×taxon×period with source abundance, modern presence, audit JSON |
| Jasprica et al., Natura Croatica 2017, Croatia 13 stations + Bosnia/Herzegovina 2 stations (ONE paper) | 15 stations | 359 source original taxa (336 species and 23 subspecies); field survey 2013–14 | Publisher original 33-page PDF SHA256 `90d489ca8cbc29464d3e6a2b1d8b5d9568c9940e953b1ce720964991366ab956`, original rotated native PDF glyph records and LF source resolved for 359 | [v18.3 SOURCE SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/43723422179), 359×15 = 5,385 presence cells, each source taxon's original station numbers, LF and family; audit JSON |

Both source PDFs are independently published and attributed; the unpublished raw PDFs are **not committed** in GitHub. Original public source-data rows are regenerated in GitHub Actions artifacts, not mistaken for `tetsu` field sampling. Check artifact retention before citing downloadable derivative CSV; permanent rerunnable code and SHA lock remain in the PR.

## Source verification beyond number of taxa

**Croatia/Bosnia original paper, Tab 2 and narrative**:

- 359 unique taxon names; 359 taxon station-number lists present; 359 life forms and families present; all 15 stations represented; site richness from source counts ranges 47–152. The 15 site counts (code 1 through 15): 129, 128, 87, 93, 152, 136, 127, 85, 151, 129, 114, 71, 62, 58, 47.
- 72 taxa present in at least 8/15 stations, consistent with publisher original results; the sole taxon present at all 15 stations was `Chondrilla juncea L.` also consistent.
- **Two conflicts remain ON THE PUBLISHER'S PRINTED PAGE** and are retained transparently:
  1. Table 2 `Canna indica L.`: explicit station `12` (1/15, frequency should round to 7%), but printed frequency `53%`; 358/359 frequency entries agree with exact station lists. DO NOT silently change station numbers to match publisher percentage.
  2. The full printed 359 named-taxon×station lists have 99 taxa occurring at exactly one station and 260 other taxa; prose states 100/259. The original 359 list is **NOT** filtered or altered to force agreement. The published 259-taxon analysis is distinguished from our all-359 taxon descriptive matrix.
- Comparison among 15 stations: 105 pairwise taxonomic Jaccard similarities from the FULL 359 list have min=0.0797, median=0.2500, max=0.5238. Functional life-form profile overlap `sum_g min(p_{ig},p_{jg})` across exactly the same pairs has min=0.6317, median=0.8971, max=0.9735.
- Full source life form category totals: T (therophyte) 183, H (hemicryptophyte) 113, P (phanerophyte) 31, G (geophyte) 17, Ch (chamaephyte) 15. A high five-class overlap is **not** cross-country functional homogenization or necessarily above chance; observed taxonomic Jaccard and functional proportional-overlap scales must not be directly compared or differenced into an 'effect size'.

**Scientific sample cautions**: The 15 stations were sampled in 2013–14 in one original paper (Croatia/Bosnia not two independent meta analytic studies); sites span macroclimates and urban-rural gradients, plots 200–800 m² and generally operating tracks/platforms. Most visited spring–autumn; 4 southern stations surveyed only summer/autumn. This uneven detection effort biases comparisons. Site/area effects and local natural habitat availability (Tab 1 distance to nearby natural/semi-natural vegetation) are available but are **NOT** yet a standardized habitat-analogue or a matched regional species pool. Source floristic species, microhabitat, openness and rare-sampling must remain separate.

## Hypothesis and precommitment

**Can replicated anthropogenic habitats override regional species pools?** A railway site array can be climate/context-stratified but is not experimental standard substrate. Outcome axes:
- Species identity similarity or phylogenetic/generic overlap (species Jaccard distant biomes saturates near zero);
- Trait/functional composition: compare group proportions or ecological distance with adequate trait coverage, not only five lifeform classes;
- Independent regional pools (GIFT regional lists potentially too coarse; vegetation plots / habitat-specific local flora if available) sampled to match survey time, size, plot area and source effort;
- A station pair is not an independent unit; use original paper/climate region/habitat blocks and leave-original-study-out validation.

A **within-railway fixed-row/column randomization** is an important sanity check, but only asks whether site-level functional similarity exceeds what the observed RAILWAY taxon pool and occupancy margins predict. It is **NOT** a natural-regional-pool null and cannot show railway-generated functional convergence.

### Next execution priorities
1. Create a validated original Tab 1 15-station macroclimate/landscape context registry; test within-paper fixed margins while retaining climate-group permutation as a conservative sensitivity.
2. Extract source-native site×taxon from Italy Parma (28 stations), Austria East Tyrol (6), CZ and SI original study site data into a harmonized data dictionary; keep exact name, species rank, native/cultivated and sampling unit semantics.
3. Build matched natural habitat **independent** regional pools, estimate observed vs expected taxonomic/functional differences, first Europe and India, before any global claim. Do not interpret a 2-paper source dataset as global independent replication.

Primary paper links: [Croatia/Bosnia original PDF](https://hrcak.srce.hr/file/282429); [Slovakia original PDF](https://www.upjs.sk/app/uploads/sites/15/2022/10/173-188_Majekova_et_al-upr.pdf); [regional sources register](../data/literature/global_rail_flora_source_registry_v180.csv).
