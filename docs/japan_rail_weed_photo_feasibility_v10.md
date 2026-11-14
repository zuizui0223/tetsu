# v1.0 modern Japanese railway-weed photo source audit (2018–2025)

**This is an independent modern citizen-observation SOURCE audit, NOT a vegetation plot census, cutting/wind treatment or evidence of plant adaptation.**

## Reproducible primary source

[GitHub Actions #42437324005, successfully executed 2026-11-08 JST](https://github.com/zuizui0223/tetsu/actions/runs/42437324005), artifact `japan-rail-historic-weed-modern-photo-coverage`. Script [`scripts/audit_japan_rail_weed_photo_coverage.py`](../scripts/audit_japan_rail_weed_photo_coverage.py), official MLIT 2025 N02 rail geometry, iNaturalist 2018-01-01 through 2025-12-31, research-grade photographed observations. The program queries ALL results for each focal taxon and locked metropolitan rectangle, not only photos with visible flowers.

The 3 predeclared rectangles are Tokyo, Osaka and Sendai (the **same bounding boxes frozen in v0.4**, not municipal boundaries). The taxa were selected from the 2011/2012 Japanese railway weeds source *before looking at the modern photo near-rail distribution*:
- *Equisetum arvense* (horsetail, spore plant)
- *Solidago altissima* (tall goldenrod)
- *Setaria* (foxtail genus, **not one species**)
- *Miscanthus sinensis*, *Imperata cylindrica*, *Digitaria ciliaris*
- *Pueraria lobata* (historic botanical name for accepted *Pueraria montana var. lobata*)
- *Artemisia princeps*, *Erigeron canadensis* (taxonomic update of *Conyza canadensis*)

Full source retrieval: **27 of 27** area×taxon-group queries complete (no more than 2000 observations in any slice, no truncation). After rejecting obscured/private, explicitly cultivated/captive, unphotographed, missing date, >100 m coordinate uncertainty or wrong taxon, **597 records** passed; of these, **229** lay ≤250 m from mapped passenger railway tracks. Counts are from opportunistic botanical *photographs*, not independent site censuses or railway plant counts.

## Species and group source opportunity

| Taxon / group | Three-city photo records passing strict spatial QC | Of those within 250 m of MLIT passenger rail | % photographic records near rail |
|---|---:|---:|---:|
| *Equisetum arvense* | **78** | **40** | **51.3%** |
| *Solidago altissima* | **177** | **56** | **31.6%** |
| *Setaria* spp. (genus, several species) | **60** | **38** | **63.3%** |
| *Miscanthus sinensis* | **25** | **16** | 64.0% |
| *Imperata cylindrica* | **60** | **26** | 43.3% |
| *Digitaria ciliaris* | **1** | **1** | 100% (unusable as a comparison) |
| *Pueraria montana var. lobata* (source query P. lobata) | **110** | **18** | **16.4%** |
| *Artemisia princeps* | **76** | **30** | **39.5%** |
| *Erigeron canadensis* | **10** | **4** | 40% (sparse) |
| **TOTAL (not individuals or study sites)** | **597** | **229** | **38.4%** |

**CRITICAL:** Do NOT compare these percentages directly with the **2011 trackbed vs slope** report-share percentages in the published Abe 2012 survey. Here, ≤250m is only a point-to-passenger-rail **proximity**, not a verified trackbed or embankment plant location, and iNaturalist observer/site access may differ dramatically by taxon, city and flowering attractiveness. The snapshot is dominated by Tokyo; Osaka and Sendai have far fewer qualifying observations. No p-values or rail-enrichment claims follow from these numbers.

## Important taxonomic QA repair

Initial version excluded **all 143** *Pueraria lobata* API-returned observations from the three areas because iNaturalist's accepted canonical name was `Pueraria montana lobata`, without the historic name/var. punctuation. These are not absence records.

We inspected original API-accepted taxon fields (recorded as aggregate `unmatched_api_accepted_taxa` in the intermediate artifact) and aligned a **narrow accepted synonym** verified by Kew Plants of the World Online: `Pueraria lobata` = `Pueraria montana var. lobata`. The revised source code accepts the specific iNaturalist canonical `Pueraria montana lobata`, but **never merges all `Pueraria montana` varieties**. Result: 110/143 kudzu photo observations pass geographic QC. Repeated source audit **597 qualified, 229 rail-near** as above. Kew source: https://powo.science.kew.org/taxon/214449-2.

## Fundamental reproductive-biology QC

**Horsetail (*Equisetum arvense*) does not flower.** The plants produce fertile cone-bearing stems (tsukushi) in spring and regenerate from rhizomes; the Flowers/Fruits annotation is biologically irrelevant to its reproduction. All 78 qualified *Equisetum* photo occurrences must be considered **reproductive phase unknown** until fertile strobili/sterile shoot structures are independently annotated. This source audit now explicitly reports `not_applicable_spore_plant` for its iNaturalist flower phase.

Use:
- `scripts/extract_photo_roi.py` for flowering angiosperms and separately reviewed flower/inflorescence ROI; the script now **rejects Equisetum flowering labels/masks**.
- `scripts/extract_horsetail_spore_roi.py` for confirmed `Equisetum` plant foreground plus **fertile spore-cone ROI** from human-reviewed aligned masks, with different output column names, no false flowers.
- Existing fixed-image plant/flower ROI and documented cutting chronology remains ready, **but no actual repeated contemporary railway management photos or event logs have been supplied**.

Independent reproductive source: https://www.qsr.mlit.go.jp/chikugo/siryo/02-kawa/01_chikugo/syokubutu/sugina.html and https://research.fs.usda.gov/feis/species-reviews/equarv.

## Study implication and next gate

The 2011 published microhabitat aggregates have a strong **observational** trackbed affiliation for *Equisetum* (24/70 trackbed species report entries, versus 22/642 slope entries). The 2018–2025 citizen-photo source also contains numerous horsetail occurrences near mapped rails, but **the sources measure different phenomena and are not an independent quantitative replication**. No railway adaptation is established.

The next independent ecological test needs **verified real rail trackbed/slope/roadside species observations and treatment records**, with matched taxon (including *Equisetum* spore stages), site-specific vegetation management (`cutting` vs `herbicide` vs `growth_retardant` vs no management), observed repeated regrowth, standardized same-view ROI and actual safety/authorization. Previous 2007 field growth-retardant work is useful biological background, but not our new data.
