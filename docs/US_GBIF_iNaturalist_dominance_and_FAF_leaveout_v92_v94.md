# tetsu v9.2–v9.4: A single citizen-science GBIF dataset dominates US recorded plant source assemblages, but removing it does not recover a general rail signal

## Actual independently verified source workflows

- [v9.2 eight fixed U.S. states, SAME GBIF dataset UUID, two years #43198013903](https://github.com/zuizui0223/tetsu/actions/runs/43198013903) · [real original source artifact #13934370511](https://github.com/zuizui0223/tetsu/actions/runs/43198013903/artifacts/13934370511) · [code](../scripts/audit_US_GBIF_matched_provider_temporal_v92.py), 8 integrity tests PASS.
- [v9.3 ALL 48 quality-qualified U.S. states iNaturalist source concentration #43199018455](https://github.com/zuizui0223/tetsu/actions/runs/43199018455) · [source artifact #13933029917](https://github.com/zuizui0223/tetsu/actions/runs/43199018455/artifacts/13933029917) · [code](../scripts/audit_US48_iNaturalist_GBIF_source_dominance_v93.py), 6 integrity tests PASS.
- [v9.4 full versus only iNaturalist versus all NON-iNaturalist source-specific 2022 crop RAIL/TRUCK association #43200225106](https://github.com/zuizui0223/tetsu/actions/runs/43200225106) · [real model artifact #13934931904](https://github.com/zuizui0223/tetsu/actions/runs/43200225106/artifacts/13934931904) · [code](../scripts/analyze_US_FAF6_firstrecord_inat_leaveout_v94.py), 6 targeted tests PASS. Full repository descriptive smoke tests PASSED at the relevant commit.

**Official GBIF iNaturalist source**: dataset UUID **50c9509d-22c7-4a22-a47d-8c48425ef4a7** = **iNaturalist Research-grade Observations**, DOI [10.15468/ab3s5x](https://doi.org/10.15468/ab3s5x), dataset https://www.gbif.org/dataset/50c9509d-22c7-4a22-a47d-8c48425ef4a7. The same UUID is officially identified in GBIF's technical code constants and in iNaturalist's data-use documentation. SAME dataset UUID does NOT imply same observers, field sites, detection process, recording protocol, historical uploads or species validity.

**Data scope:** independently sourced GBIF 2026 interpreted positive indexed speciesKey record facet counts, U.S.-country-GRIIS-listed plant taxa, Census FIPS state names/area/geography, occurrence event years **2018–2021 before** and **2023–2025 after** (2022 excluded), original 2022 U.S. BTS FAF6 SCTG02 cereals and SCTG03 other agriculture economic origin→destination estimated crop RAIL and TRUCK freight. No actual seed contamination or local botanical establishment measurements.

## Strong REAL macro recording-process result

Eight v8.1 FIPS-preselected states (CA FL IL IA NY NC TX WA): select top 3 common dataset UUIDs purely by min(pre,post) ALL-Plantae source count, require ≥50 Plantae source records in EACH window, and ≥20 U.S.-country-GRIIS indexed source records in both to count a provider as usable. These same provider UUIDs retained median **99.7956%** of later ALL-Plantae GBIF indexed records across states and contained median **97.3078%** of later-only U.S.-country-GRIIS speciesKey categories recorded by all providers. iNaturalist was the **largest dataset in 8/8** states.

**For ALL independently prequalified 48 states:**

| State-level median share of ALL-Plantae GBIF source records from exactly one iNaturalist research-grade GBIF dataset UUID | 2018–2021 | 2023–2025 |
|---|---:|---:|
| Same 48 Census-state recording sources | **95.8247%** | **99.1155%** |

- iNaturalist's share increased in **44/48 states**; late-state shares range **91.9157%–99.9773%**.
- State median **94.5745%** of later-only U.S.-country-GRIIS source speciesKey categories appeared in the LATE iNaturalist dataset source too.
- **40/48** state sources retain ≥20 NON-iNaturalist GRIIS indexed RECORDS in both windows. RECORDS do not automatically imply ≥5 unique taxa or independent observers.
- **44/48** state sources have at least one later-only **within-iNaturalist** speciesKey category that already appeared earlier in ANOTHER provider's GBIF source. A provider-specific first record therefore can be a previously reported taxon.
- The earlier independently audited v9.1 national recording shift is related but distinct: per-year indexed plant source volumes increased in **48/48 states**, median **2.1426× all Plantae / 2.0640× GRIIS**. The contribution of new observers, backfills, true biological change and taxonomic changes cannot be separated here.

**Interpretation:** iNaturalist dominates the recorded U.S. plant source patterns. Matching the dataset UUID across periods does not create equivalent surveys. The large numbers of later-only source taxa are NOT estimates of true recent plant immigration.

## v9.4: Exact non-iNaturalist per-taxon source reconstruction

The original all-GBIF state-by-year source includes an indexed RECORD COUNT FOR EACH speciesKey, and the independently queried single UUID source contains the corresponding speciesKey record counts under the SAME GBIF filters. Thus for each GBIF speciesKey, state and year period, **NON-iNaturalist positive indexed source count = all-provider count MINUS iNaturalist-only count**. A priori numeric conservation requires all taxon counts to remain nonnegative and to sum to all original GRIIS source counts. This yields exact non-iNaturalist recorded taxon count *facets*, not a verified original independent observer dataset.

**State eligibility BEFORE inspecting ANY freight/plant coefficient:** from original 48 quality-qualified U.S. states, require in BOTH 2018–21 and 2023–25:
- ≥100 NON-iNaturalist all-Plantae indexed source records;
- ≥5 recorded U.S.-country-GRIIS speciesKey categories;
- ≥20 U.S.-country-GRIIS indexed records.

**40/48** states qualify. The eight source-ineligible states are Arkansas, Connecticut, District of Columbia, Kentucky, Maine, Oklahoma, Rhode Island and Vermont. They are neither ecological zeroes nor countries without native/introduced plants. Same selected 40 states are used for FULL-GBIF, iNaturalist-only and all-OTHER-GBIF source variants. Intersect valid first-recorded source candidate origin→destination pairs across the THREE source variants BEFORE calculating any source coefficient, yielding exactly **1,528** common directed state comparisons (which reuse the 40 source state datasets, and are NOT 1,528 independent ecological experiments).

**Outcome**: fraction of source-origin prior-recorded U.S.-country-GRIIS taxa (recorded in origin 2018–21 but not recorded at destination in 2018–21) source-recorded at the destination in 2023–25 under a 20-indexed-record inclusion/rarefaction model. Not true immigration, state-introduction or an independent plant absence survey.

**Model**: Partial rank correlation of 2022 crop RAIL economic OD volume to this recorded taxon source response; same-direction 2022 crop TRUCK OD, Census distance/region/area, GBIF recording-effort source, candidate species pool size, and source/destination crop freight node strength controls. For the PRIMARY comparison ALL source botanical covariates are held fixed, so the recorded response is the changed source dimension; secondary fits use each source-specific dataset's own recording effort covariates.

### Real results on EXACT same source pairs, with unchanged original all-GBIF covariates

| Source of counted plant record taxonomy (2018–21 origin / 2023–25 destination) | 2022 crop RAIL partial rho | 2022 crop TRUCK partial rho |
|---|---:|---:|
| ALL GBIF datasets | **−0.027949** | +0.047185 |
| iNaturalist-only GBIF original dataset | **−0.021543** | +0.044553 |
| ALL other GBIF datasets, excluding iNaturalist | **−0.055712** | +0.003699 |

The entire 40-state/1,528-directed source pair universe has **no positive RAIL correlation** in any of the three recorded-source variants; this does not establish a true negative causal transport effect or disprove railway seed dispersal. With each provider-specific source's OWN record controls, primary RAIL rho remains negative: all −0.027949, iNaturalist −0.021099, other −0.042011.

**MANDATORY conditional cohort/recording-source fragility:**

| SAME original 40-state geography and controls | All-GBIF RAIL rho | iNaturalist-only RAIL rho | All-OTHER RAIL rho | Directional pairs |
|---|---:|---:|---:|---:|
| ALL eligible source directions | **−0.027949** | **−0.021543** | **−0.055712** | **1,528** |
| Crop RAIL economic OD >0 | **−0.041313** | **+0.025003** | **+0.031110** | 128 |
| Both crop RAIL and TRUCK economic OD >0 | **−0.138917** | **−0.119155** | **+0.041179** | 106 |

Thus apparently favorable positive railway source signals on small, selected freight-connected pairs are **highly recording-source dependent**; choose neither a favorable source nor freight-cohort selection after inspecting the signed coefficient. The raw unrarefied source outcome differs yet further; original model outputs include it.

## Scientific conclusions

**Supported DATA results:** 2022 U.S. agricultural RAIL estimated freight has a highly concentrated, interstate-directed economic OD structure (independent FAF6 source); 2026 GBIF historical plant occurrence records are dominated by iNaturalist; a provider-excluded 40-state recorded source model is technically possible and gives no general positive RAIL freight/late first recorded plant correlation.

**NOT supported:** that actual viable alien plant seeds rode trains, were released at a particular destination, established in a true local alien state population, spread following a specific operated physical freight rail corridor, responded to rail mowing/ballast thermal/wind gust stress, changed plastic morphologically or adapted genetically.

**Crucial nonidentified sources:** 2022 economic OD shipment estimates are not observed track flows, early absence from source records is not biological absence, GRIIS country alien status is not state-level status, GBIF provider UUID does not identify observer protocol or event-independent plots, 4-year before vs 3-year after source windows differ, and state pairs share data. No valid pairwise p-values or cross-national replication of biological rail effects.

Research recommendation: retain both favorable small-subset and unfavorable nationwide estimates in Draft PR, and do not promote railway dispersal/adaptation to a confirmed ecological finding until independent original event provider/observer-fixed rail-yard and road/rail verge vegetation, seeded cargo/seed traps, dated active freight track routes and measured disturbance/fitness or common-garden evidence exist.
