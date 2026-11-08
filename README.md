# tetsu — Railway-mediated plant redistribution

**Research question:** Do railway networks connect distant plant communities through propagule transport, beyond railway-side habitat, roads, urbanization, ports and recorder bias?

**Status: 2026-11-07. Research design + exploratory code only. No empirical ecological results.**

## Mechanisms to distinguish
- **Transport:** seeds on freight/rolling stock, cargo contaminants/spillage, people and maintenance equipment.
- **Habitat:** disturbed railway verges allow ruderal/alien species to establish without any train-carried seed.
- **Confounding:** stations co-occur with ports, roads, cities, trade, gardens and heavy botanical observation.

## Predicted signatures
1. **Rail habitat:** railway sites have elevated alien plant occurrence/richness compared with matched road and background habitats.
2. **Rail transport network (primary):** connected, active railway regions share more alien plant species than otherwise comparable unconnected regions, after controlling for road travel, geographic distance, urbanization, climate, and recording effort; contrast native plant similarity.
3. **Cargo/ports:** signatures should be stronger for freight/port-linked regions and accidentally introduced taxa than for passenger-only rail and ornamental pathways.
4. **Historical test (conditional):** newly connected regions gain species records after line openings compared with credible temporal controls; first record is NOT arrival.

An association does not prove that a seed travelled by train.

## Real-data audit (2026-11-07)

[Inspect the measured data feasibility and caveats](docs/empirical_feasibility_2026-11-07.md).

The first GitHub Actions runs successfully retrieved the Japanese railway GIS, the **full 1,753-taxon** national alien plant checklist, official JR Freight container-handling site candidates, and GBIF occurrence counts. These are **data availability results, not positive rail-mediated dispersal results**.

- N02: 21,933 rail-section features and 10,234 station line features; freight-only lines/terminals omitted.
- 139 Japanese alien taxa have accidental-only introduction pathway flags; 351 have ornamental-escape flags without accidental ones.
- JR Freight web listing: 120 candidate sites, including 19 ORS (off-rail service), which could serve as road-logistics contrasts after validation.
- GBIF recently georeferenced Japanese plants: 249,745 records across **Plantae**, including nonvascular plants; name/effort screening is required before biological inference.
- Small pathway-stratified GBIF audit: only 4/12 accidental-group and 3/12 ornamental-group sampled taxa have >=20 qualifying recent records, showing that complete floristic incidence cannot be assumed.

Run the live source audit in GitHub Actions workflow `.github/workflows/public-data-audit.yml` or manually:
```bash
pip install requests openpyxl pyshp beautifulsoup4
python scripts/audit_public_data.py --strict
python scripts/summarize_alien_checklist.py
python scripts/audit_freight_stations.py
python scripts/audit_gbif_taxa.py --names-per-group 12
```

## Additional railway environmental filtering and plant phenology workstream (v0.3)

**New question:** Do railway mowing, train slipstreams, and hot/dry ballast conditions shift plant architecture, regrowth and flowering schedules, separately from trains moving seeds?

- [Mowing × slipstream × life-history hypotheses and interpretation limits](docs/rail_environment_traits_protocol.md)
- [Photo ROI and fixed-plot phenology input contracts](docs/photo_and_phenology_inputs.md)
- [Silhouette and floral ROI extraction from verified PNG masks](scripts/extract_photo_roi.py): width/height aspect ratio, pixel-corner silhouette solidity, visible flower fraction, position of flowers *within* plant silhouette, and explicit photo-framing QC. **Absolute plant size or biomass cannot be recovered from uncontrolled snapshots.**
- [Weekly observed flowering from complete repeated visits](scripts/summarize_photo_phenology.py): blooming vs nonblooming structured observations, observer/site effort, confirmed last cutting dates, and warning gates. **First flowering photo is not true flowering onset.**

```bash
pip install Pillow
python -m unittest discover -s tests -v
python scripts/extract_photo_roi.py --manifest your_manifest.csv --root path/to/PNG_masks --output results/photo_roi.csv
python scripts/summarize_photo_phenology.py --visits your_complete_visits.csv --events verified_cut_events.csv --output-dir results/phenology
```

These modules are input-driven methods **without real railway photo-treatment observations or any ecological effect estimates**. Within-individual repeated measures give more compelling evidence of a plastic *response* to environmental change; heritable adaptation requires additional common-garden/genetic evidence, and photo records alone cannot measure wind load.

## First real rail-proximity phenology and ROI photo audit (v0.4, 2026-11-07)

**[See the real-data report and exclusions](docs/real_photo_feasibility_2026-11-07.md)**, [fixed API/GIS protocol](docs/inat_rail_phenophase_pilot.md), and [licensed photograph review](docs/real_photo_roi_review_2026-11-07.md).

- [iNaturalist annotation × 2025 passenger railway GIS audit](https://github.com/zuizui0223/tetsu/actions/runs/42373021202) completed. Tokyo/Osaka/Sendai, 8 taxa, 2018–2025: **185 photo observations annotated Flowers, only 4 annotated No Flowers or Fruits** (possible overlap), 132 geo-quality passing, 37 within 250m of mapped rail. **Do not infer a railway flowering shift**, since nonflowering labels are far too sparse and are visibility statements, not repeated phenology visits.
- [Independently sampled CC0/CC-BY photograph review](https://github.com/zuizui0223/tetsu/actions/runs/42373519458) downloaded **13** licensed images of *Trifolium repens* and *Solidago altissima* from Tokyo/Osaka, across rail-near/rail-far strata. Photos and credits are retained as CI artifacts, not committed to the repo. Initial 280-px visual review: **0 clearly isolatable entire individuals; 1 ambiguous candidate; 7 multi-individual patches**. This is method feasibility, not trait selection.
- [Provisional photo-by-photo review](docs/photo_2026-11-07_provisional_review.csv) explicitly separates photographed patches/flower closeups from whole individuals.
- [ROI manifest extractor](scripts/extract_photo_roi.py) now requires `target_scope` and `mask_reviewed`: **patches and unreviewed masks cannot be reported as valid individual size/shape traits**, even if pixel area is computed.
- Rail GIS proximity **does not establish a railway-side microhabitat**, mowing history or slipstream load. None of the 13 reviewed image frames visibly confirmed track/trackbed.

**Next warranted dataset:** geographically replicated, same-species, independently site-verified photos at fixed plant/patch scale; repeated known flowering and nonflowering visits around documented mowing events. Actual railway effects, plasticity and adaptation **remain untested**.

## v0.5 verified Tokyo official-road control and mowing patch ROI (2026-11-07)

[Tokyo road-vs-rail empirical audit and corrected analysis](docs/real_rail_vs_road_tokyo_2026-11-07.md) · [Fixed analysis design](docs/rail_road_control_pilot.md)

- **Official N13 road GIS added** (2024 data): 1,943,251 features in the first Tokyo map mesh, 711,724 local road line parts indexed against N02 passenger rail.
- Real iNaturalist research-photo occurrences 2023–2025: 411 high-quality focal species records, 3,315 sampled other-plant reference records; 246 and 882 distinct observers respectively.
- Road proximity ≤100 m; rail-near ≤250 m vs rail-far ≥1 km: **124/134** focal photos near/far. Exploratory source-weighted composition ratio **0.1333**, month-stratified 0.1586; this is photographic recording composition and **not** a biological negative railway effect.
- **Spatial matching gate failed:** of 44 2-km cells with focal and reference photos, only 2 have both focus/reference across BOTH rail-near and rail-far. **Do not claim railway dispersal, adaptation or flowering shift based on this source alone.** Original incorrect short-page weight reported 0.1159 and was repaired; all figures above use corrected probability weight.
- [Patch ROI new valid analysis mode](docs/mowing_roi_event_study.md): fixed camera, independently reviewed **conspecific plant patch mask + flower mask**, repeated plot ID and documented mowing dates. Can calculate projected vegetation occupancy and flowering fraction changes without inventing physical height; `patch_pass` cannot be mistaken for valid whole-individual shape.
- [Paired-photo mowing CLI](scripts/summarize_mowing_photo_roi.py): requires ≥2 pre-cut fixed photos in days −28:−1 and ≥2 recovery photos in +7:+42 days; overlapping cuts and unknown photo quality are rejected. **Code has synthetic tests, not real documented cut measurements.**

## v0.6 — actual repeated flowering records and rail-site eligibility (2026-11-07)

[Real-data results, methods and caveats](docs/usa_npn_real_phenology_rail_v06.md)

- **USA-NPN** (actual monitored individual flowers present/absent): 2024, MD tuliptree 3,229 multi-phenophase records from 11 plants/3 sites, AZ saguaro 12,957 records from 51 plants/23 sites. These include many observed **no** entries, unlike the iNaturalist photo-annotation pilot.
- [Individual year status audit](https://github.com/zuizui0223/tetsu/actions/runs/42379297172): MD tuliptree Open flowers 36 yes/457 no, 10 plants with observed first-yes bracketed by previous no ≤21 d; AZ saguaro 174 yes/2,451 no, 26 bracketed first-yes and **9 potential yes → 2 no visits → yes observed sequences**. No data establishes that railway mowing or wind caused these natural flower patterns.
- [US FRA railway overlay](https://github.com/zuizui0223/tetsu/actions/runs/42379467494): **all 26** monitored sites in these two 2024 taxon samples were **>1 km from mapped rails**. They cannot identify a rail effect; do not mistake a large phenology dataset for a spatially balanced railway experiment.
- [Independent registered station proximity census](https://github.com/zuizui0223/tetsu/actions/runs/42379733758): of **1,041** MD and **2,918** AZ NPN registered sites, **105 MD + 65 AZ = 170** lie within 250 m of FRA rail. **Registered is NOT verified active with flowering data in 2024**. This provides a pre-exposure sampling frame for finding actual repeated-observation sites by species before effect testing.
- Reproducible Python source clients, strict yes/no/conflict rules, interval-censoring, future railway proximity eligibility and tests are under `scripts/audit_usa_npn_*.py`, `scripts/run_npn_individual_phase.py` and `src/tetsu/npn_individuals.py`. No private site/individual coordinates released in CI artifacts.
- The Japan freight/alien-network hypothesis and independent fixed-photo mowing/flower-ROI workflow remain separate and **unresolved**, not rejected by the USA-NPN source-feasibility diagnostics.

## v0.7–0.8 — Actual near-rail 2024 NPN census and multi-year matched flowering records

**This is evidence of DATA FEASIBILITY, not evidence of railway-related flowering, seed dispersal, mowing plasticity or wind adaptation.**

- [Full 2024 registered-site census](docs/npn_active_rail_2024_results_v07.md): all **170** NPN sites registered within 250 m of FRA rail in MD/AZ queried for actual 2024 status. Only **7** had any NPN records (3 MD + 4 AZ); only **6** any Open flowers record; only **2 AZ** near sites had at least one same-taxon individual with >=4 observed dates and both YES and NO. All 170 queries successful, zero network/API errors. Prior registered-site counts are NOT verified active monitor counts.
- [Multi-year site-phenometrics source audit](docs/npn_site_pheno_2023_25_v08.md): official USA-NPN national Open flowers first-YES site summary, 2023–2025. **Corrected a silent API state-filter failure**: one nationwide source per year, then independent MD/AZ registered-site-ID join (do not describe 3,770/3,670/3,561 nationwide rows as state totals).
- [Same-taxon/source + geographic controls](docs/npn_candidate_status_verification_v08.md): matching first by phenology data availability and <=75km geographical proximity, and THEN refetching independent NPN individual YES/NO time series, identified *Cercis canadensis* (MD) and *Prosopis velutina* (AZ). Each has only **2 rail-near independent active sites** in its best year(s). These are **woody trees**, not herbaceous plants for railway grass-cutting adaptation.
- [Official individual source run](https://github.com/zuizui0223/tetsu/actions/runs/42390743089) verified actual same-species monitored YES/NO at some matched far controls in five state/year cases; no statistical railway-effect estimators are reported. The later tentative interval midpoint summaries, if any, remain **site-level exploratory descriptors**, not biological onset or selection effects.
- [Observed onset-interval midpoint screening](docs/npn_observed_onset_midpoint_v08.md) on these outcome-independent selected woody NPN cohorts gave unadjusted near-minus-far days: *Cercis canadensis* **-4.75 (2023), -7.0 (2025)**; *Prosopis velutina* **+4.0 (2023), +0.75 (2024), +7.75 (2025)**. **Only 2 rail-near independent stations per year** and 1–4 valid far sites; date-interval width often 5–15 days, and climate/urban/observer variables and trackside plant exposure unverified. These are not demonstrated biological onset shifts or railway effects. No p-values.
- The **Japan freight/port/alien dispersal** main study and **fixed photo patch flowering ROI / actual mowing event** implementation remain valid independent research tracks. Do not substitute NPN trees or railway station distances for independently verified railway verge plants.

## v0.9 — independent Daymet GDD and Japanese maintenance contrasts (2026-11-08)

**[Real site-level thermal-time/source audit](docs/npn_daymet_thermal_control_v09.md)** · [GitHub Actions source result](https://github.com/zuizui0223/tetsu/actions/runs/42430207432)

- Independently sourced ~1km Daymet daily air temperature and precipitation joined to the **same preselected, outcome-blind NPN rail-near/far species×year cohorts**. All **16/16 unique site climate API requests succeeded**, with source errors = 0. Source variable schema, 365-day completeness, GDD calculations and fixed cutoff-day tests in `src/tetsu/daymet.py`.
- **Prebloom GDD5 contrasts (rail-near minus geographically matched rail-far sites):** *Cercis canadensis* Maryland DOY60: **−10.165 °C days** in 2023 and **+15.400** in 2025; *Prosopis velutina* Arizona DOY90: **−31.310** in 2023, **−7.885** in 2024 and **+20.310** in 2025. Each has only **2 independent near-rail sites** per year and 2–4 far controls. A 2026 FRA rail line near an NPN station **does not establish railway habitat or train-wind exposure**. GDD to day 120 can occur after first observed flowering and cannot be treated as a first-onset causal covariate.
- The previously observed first YES interval-midpoint timing differences and this source-backed temperature variation do **not** form a consistent railway-specific response or adaptation signal. Sites are too few, far comparators differ across some years, phenological onset is interval-censored, and climate/road/land cover are not balanced.
- **[Japan railway weed management historical baseline](docs/japan_rail_mowing_historical_baseline.md):** Ito et al. (1982) observed different distributions under repeated herbicide versus cutting, plus urbanization and embankment context; this is published observational evidence **not our own new cutting experiment**.
- Original freight/port seed transport, mowing-induced plasticity, train wind, flowering ROI and site-level adaptation remain separate scientific hypotheses, not confirmed or disproved by these environmental source audits.

## v1.0 — source-backed Japanese railway plant microhabitat filtering (2026-11-08)

**[New quantitative secondary-data reanalysis](docs/japan_rail_microhabitat_reanalysis_v10.md)** · [published-source table CSV](data/literature/abe_2012_railway_microhabitat_occurrences.csv) · [reanalysis code](scripts/reanalyze_abe_2012_rail_microhabitat.py) · [successful 6-test original-data reconciliation](https://github.com/zuizui0223/tetsu/actions/runs/42436745865)

- Reconstructed Abe et al. (2012) JR East 2011 **109-section** weed questionnaire Table 4, by warm/cold geographic grouping and **six rail microhabitats**: **1,334 taxon×habitat reported occurrences** (not independent individuals or survey plots), with **70 trackbed** and **642 embankment slope** occurrence mentions. All 22 taxon/region row totals, microhabitat columns and regional subtotals reconcile to original published table.
- **Original descriptive finding from published aggregates:** *Equisetum arvense* (スギナ) accounts for **24/70 = 34.3% of all trackbed plant mentions**, compared with **22/642 = 3.4% of slope mentions**: a **10.0-fold ratio of relative report shares**. The direction is present separately in warm (**13.5×**) and cold (**8.9×**) railway-region subgroups. Across all 1,334 occurrence mentions, horsetail accounts for just **98/1,334 = 7.35%**, giving a **4.67× relative trackbed overrepresentation**. **This is a newly computed OBSERVATIONAL composition ratio, not a causal adaptation estimate or new field experiment.**
- *Solidago altissima* (セイタカアワダチソウ) has 0 reported trackbed and 46 slope occurrences in the source Table 4, despite 90 reports across all habitats; **not species absence from all railway property**. Source taxa missing from one geographic top10 list are NOT coded biological absences.
- Independent 2002–2004 Okayama railway embankment published **growth-retardant** trial reported selective loss of large goldenrod/mugwort, replacement by smaller foxtail/crabgrass and ~30% less aboveground biomass by year 3, but this is NOT a direct test of wind or a national mowing-vs-herbicide comparison. 1982 Japanese railway survey described management-specific Setaria/Solidago/Imperata differences. [Evidence provenance](docs/japan_rail_microhabitat_reanalysis_v10.md).
- **[New real modern citizen-photo source audit](docs/japan_rail_weed_photo_feasibility_v10.md):** Tokyo/Osaka/Sendai in fixed city rectangles, 9 rail-weed taxa/groups, 2018–2025 iNaturalist photographic occurrences. All 27 queries complete; **597 geographic-quality-passing records, 229 within 250m of MLIT passenger railway lines**. **Modern photographic proximity is NOT rail vegetation, plant abundance, flowering onset, mowing exposure or independent confirmation of the 2011 survey.** Original accepted name `Pueraria montana lobata` had falsely excluded kudzu from initial counts; now Kew-verified narrowly matched.
- **[New horsetail spore ROI](docs/horsetail_spore_roi_inputs_v10.md):** *Equisetum* does not make flowers or seeds. Angiosperm flowering-mask pipeline now rejects Equisetum flower labels; [separate spore-cone/sterile shoot ROI](scripts/extract_horsetail_spore_roi.py) analyzes human-verified PNG masks without ever relabeling cones as flowers. Output includes projected spore-cone fraction, and discriminates individual vs fixed-plot patches. All new code has unit tests, no modern spore-cone field measurements claimed.

**Important mechanistic separation:** *Equisetum* cannot be explained by a **seed-only** railway-dispersal process, although its spores or rhizome fragments could potentially move. Its railway microhabitat association therefore sharpens the independent **establishment/management filtering** question, rather than confirming train-carried seeds. The original **angiosperm seed dispersal / freight / ports** analysis remains unresolved.

## Documentation
- [Study plan and falsifiers](docs/protocol.md)
- [Open data and important limitations](docs/sources.md)
- [Pilot CSV input schema](docs/pairwise_inputs.md)

## Exploratory pilot (synthetic example; not a biological finding)
```bash
python scripts/run_pairwise.py --incidence examples/synthetic_incidence.csv --pairs examples/synthetic_pairs.csv --output results/synthetic
python -m unittest discover -s tests -v
```

**Japan-first data strategy:** MLIT N02 + Egawa and Koyama's Japan alien plant list (1,753 taxa) + GBIF occurrences + independently validated rail freight/yard/port layers. Critically, MLIT N02 excludes freight-only lines and yards; do not equate its rail network with actual freight transport.

All real analyses must preserve the GBIF download DOI, raw-source license/version/date, coordinate/taxonomic cleaning log, and spatial-blocked validation.