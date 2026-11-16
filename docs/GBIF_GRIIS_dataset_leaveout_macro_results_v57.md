# tetsu v5.7 — Full introduced-plant source vs leaving out the single dominant GBIF dataset

**All analyses completed using real GBIF indexed source facets and original frozen Overture 50km physical RAIL/ROAD graphs, November 10, 2026.** This is a source-dependence diagnostic. It is NOT evidence of actual railway-carried seeds, local railway verge colonization, historical freight operations, or plant evolution.

## New study question

The six-country GRIIS-introduced plants and physical RAIL/ROAD observations previously suggested a hypothesis-compatible sign in **5/6 country blocks at K80**, but not K20, in a posthoc source-record rarefaction sensitivity. The v5.0 GBIF DATASET_KEY source audit showed strong source concentration, especially US/Canada. Could the apparent plant-similarity signal survive a much stricter *source-removal* control, not just adjustment for between-cell dataset overlap?

We used the ORIGINAL v3.6 same-country introduced speciesKey list and 4×4 EPSG:6933 50km equal-area cells across US Midwest, France Burgundy, Canada Ontario, Britain Yorkshire, Germany Saxony and Australia Southeast; country GRIIS checked source-accepted plant keys, fixed two periods 2010–17 and 2018–25, and real independent Overture RAIL/ROAD physical source shortest routes.

## 1. Real country GRIIS alien SPECIES_KEY × GBIF datasetKey source partition (v5.6)

[Six successful actual indexed GBIF source country jobs #42831989148](https://github.com/zuizui0223/tetsu/actions/runs/42831989148), each passing 6 synthetic source/provenance tests. For each of the **83** previously GRIIS-introduced source-eligible 50km cells in BOTH periods, GBIF source speciesKey facets were completely fetched for (a) all source datasetKeys and (b) ONLY the dominant source DATASET_KEY identified *independently from all-plant recording frequency* (not selected based on the rail/plant outcome). The complement source speciesKey counts (c) all OTHER datasetKeys were derived by exactly subtracting (b) from (a), with strict count and current DETECTED/NOT_DETECTED source marginal consistency and a maximum 10% old-vs-current record count drift gate.

The new data are **current source species-count marginal arrays, not DOI-downloaded original occurrence rows**. The current source GBIF API does not permit reliable legacy occurrenceStatus=PRESENT filters; status facet checks were included but a small amount of NOT_DETECTED cannot be assigned to individual speciesKey without source event-level records.

| Fixed country macro block | Both-period original GRIIS source-eligible 50km cells | Cell-median fraction of 2018–25 recorded country-GRIIS introduced speciesKey RECORDS derived from the largest all-plant DATASET_KEY | Median fraction of 2018–25 introduced speciesKey categories with records in BOTH largest and remaining dataset groups | Cells retaining >=5 introduced keys and >=20 source records in BOTH periods after top dataset removed |
|---|---:|---:|---:|---:|
| US Midwest | 15 | **0.844785** | 0.277487 | 15 |
| France Burgundy | 16 | **0.354626** | 0.349768 | 16 |
| Canada Ontario | 10 | **0.893993** | 0.280275 | 7 |
| UK Yorkshire | 12 | **0.517571** | 0.572355 | 12 |
| Germany Saxony | 14 | **0.477360** | 0.404142 | 14 |
| Australia Southeast | 16 | **0.411828** | 0.335828 | 15 |
| **Total** | **83** | **Country/source composition varies** | **Dataset groups do NOT mean independent observers** | **79/83** |

**Positive observational data-quality finding:** Country-listed introduced flora RECORDS, not just the entire all-plant record corpus, are highly concentrated in one published dataset in US/Canada. A majority of indexed country-GRIIS introduced species keys lack supporting records in both source-dataset groups in **5/6 country regional medians**. That does not establish that those taxa are falsely identified or absent; it means multi-source corroboration for observed regional taxon categories is limited, and species source lists are not equivalent to repeated plot detections.

[Source code](../scripts/audit_gbif_griis_dominant_dataset_leaveout_v56.py), [strict facet tests](../tests/test_GBIF_GRIIS_dominant_dataset_leaveout_v56.py), [real source workflow](../.github/workflows/GBIF-GRIIS-top-dataset-leaveout-v56.yml).

## 2. Real same-50km-cell-pair comparison of full vs non-dominant source taxa (v5.7)

[Actual six-country source comparisons #42832865840](https://github.com/zuizui0223/tetsu/actions/runs/42832865840), SUCCESS, five synthetic source quality tests. Within each K=20,40,80 record-depth comparison, use *exactly identical botanical source cells, RAIL/ROAD physical shortest paths, and geographic+ROAD+full-GBIF indexed record scale+dataset overlap covariates*, then change only the RESPONSE:
- Full country GRIIS-listed introduced-plant records from all source datasets.
- GRIIS introduced source-records after subtracting the ALL-PLANT dominant dataset in that cell/period.

Cells must have >=K indexed introduced records and >=5 unique GRIIS listed introduced speciesKey categories in BOTH time periods, AFTER dominant dataset removal. For every period, route anchor strategy and rarefaction definition are held unchanged across the two responses. Original source country block is the aggregate independent unit, not shared-cell pairwise observations.

### Late 2018–2025 regional cross-source reproducibility — primary source quality result

| Country | K20 original vs without dominant dataset, same-pair rho | K80 original vs without dominant dataset, same-pair rho | K80 50km cells / paired routes after source exclusion |
|---|---|---|---|
| US Midwest | **−0.263712 → −0.012227** | **−0.415063 → −0.212381** | 9 cells / 36 dependent pair routes |
| France Burgundy | **+0.056691 → +0.078093** | **−0.059143 → +0.101852** | 16 / 120 |
| Canada Ontario | NOT ESTIMABLE (7 cells, 10 pair routes) | NOT ESTIMABLE (2 cells, 0 routes) | 2 / 0 |
| UK Yorkshire | **+0.056952 → −0.047043** | **−0.020952 → −0.295825** | 12 / 36 |
| Germany Saxony | **−0.069611 → −0.096261** | NOT ESTIMABLE (7 cells, 21 pair routes; minimum >=8 cells) | 7 / 21 |
| Australia Southeast | NOT ESTIMABLE (15 cells, only 15 pair routes) | NOT ESTIMABLE (12 cells, only 10 pair routes) | 12 / 10 |

A negative rho means relatively shorter physical RAIL vs ROAD source centroid routes are associated with greater recorded introduced-flora source similarity, conditional on the same source covariates. **No p-value**, because the 20–120 geographical comparisons within one macro region share source cells.

### Equal-country descriptive results on SAME K-specific source pairs

| K indexed RECORDS | Eligible country blocks after removing dominant source | Hypothesis-compatible negative directions among those blocks AFTER removal | Country median partial rank with all datasets, on identical pairs | Median after dominant dataset excluded, same pair set |
|---:|---:|---:|---:|---:|
| 20 | 4/6 | 3/4 | **−0.006460** | **−0.029635** |
| 40 | 4/6 | 2/4 | **−0.073495** | **−0.055990** |
| **80** | **3/6** | **2/3** | **−0.059143** | **−0.212381** |

**Crucial comparison:** the old K80 five-of-six apparent railway-compatible signs do NOT remain broadly testable when the dominant source dataset is physically removed. Only 3 country blocks meet the original source coverage and matched railway/road route thresholds, and **France flips from negative to positive** on the exact same K80 geographical pairs. Germany and Canada lose GRIIS source record depth; Australia's mapped rail/road pair overlap is too small. A median of 3 country regions does NOT imply replicated independent macroecological railway dispersal.

Also at K20, the original US partial rank almost vanishes when the dominant all-plant source dataset is removed (**−0.264→−0.012** on identical source pairs). This directly demonstrates that an apparently strong single-country observational signal is sensitive to the record-providing dataset, not just to the record-count control.

Reproduce [source matched model](../scripts/analyze_GRIIS_topdataset_leaveout_rail_macro_v57.py), [strict shared-pair tests](../tests/test_GRIIS_topdataset_leaveout_rail_macro_v57.py), and [original six-country combined workflow](../.github/workflows/GBIF-GRIIS-dominant-dataset-rail-leaveout-v57.yml).

## 3. What these findings do and do not establish

**Supported:**
- Regional country-GRIIS-listed introduced plant *records* often depend heavily on one source dataset; across country blocks, source-category corroboration and dominant dataset footprint differ strongly.
- Source dataset removal materially affects available comparison geography and signed descriptive railway–floristic record correlations.
- Previously proposed K80 biological pattern was posthoc, and the current independent region/data-provider sensitivity is insufficient to establish a robust generalized sign.

**Not established:**
- That railway lines actually transport species propagules, accelerate colonization, select wind/mowing-resistant genotypes, or shift flowering phenology.
- That source taxa observed in only one GBIF dataset are necessarily erroneous or biologically absent elsewhere.
- That a dataset UUID is an independent human observer, standardized ecological plot, or a full model of detection effort.
- That the physical 2026 Overture graph accurately represents active, directed or frequency-weighted 2010–2025 freight movements.
- That the 3 K80 surviving blocks have any statistically significant country-independent plant–rail mechanism.

### Next real methodological requirement

True independent-source biological confirmation requires **DOI-citable complete original GBIF DWCA event rows**, verified row-level DETECTED status, event/record duplicate/observer/dataset and geographic precision checks, comparable true terrestrial biome/ecoregion strata, and independently verified historic active freight service or actual station/rail-vehicle propagule measurements. Source GBIF DWCA request manifests for the original 10 blocks were prepared in v5.5, **but not submitted**, and no DOI/download was claimed. See [official GBIF request-preparation instructions](rail_macro_source_readiness_v55.md).

**Manuscript stance:** The clear positive result is a network and recording-source DATA-QUALITY finding, not railway-driven alien-plant dispersal. Any ecological source relation remains exploratory, selected-K-sensitive and noncausal.
