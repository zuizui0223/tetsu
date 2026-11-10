# tetsu v5.2 — Recording-source concentration, dataset overlap and fixed-cohort rarefaction

**Actual data provenance: November 10, 2026 JST.** This is a methodological investigation of why the apparent relation between physical rail accessibility and country-GRIIS-listed introduced plant speciesKey similarity changes with record standardization. There is **NO established train-carried seed-dispersal mechanism**.

## Research question

In pre-existing six 200km geographically independent country source blocks (US Midwest, France Burgundy, Canada Ontario, UK Yorkshire, Germany Saxony, Australia Southeast), do source dataset overlap and rarefaction depth explain the sign instability of the observed introduced-plant floristic similarity vs physical RAIL/ROAD route rank correlation?

The geographic frame was fixed before inspecting these outcomes: 50km equal-area EPSG:6933 cells, 4×4 per source region. The original v3.6 country GRIIS introduced flora source data were used unchanged. The 2026 physical rail/road path source and two anchor strategies were unchanged. One country block supplies one descriptive coefficient: pairwise cell comparisons share sites, so no naïve pair-level p-values.

## 1. Actual new GBIF DATASET_KEY source distributions, 83/83 eligible cells

[Actual six-country full-source DATASET_KEY and OCCURRENCE_STATUS facet run #42796286437](https://github.com/zuizui0223/tetsu/actions/runs/42796286437) **SUCCESS**: 6/6 country workflow jobs succeeded, each passed six synthetic correctness tests; every source-eligible GRIIS-introduced flora cell had its historical periods' datasetKey facets fully paginated, and all passed source marginal data quality gates (previous index count drift <=10%, known dataset fractions >=98%, non-DETECTED/unknown status <=2%). **83 geographically distinct 50km cells × 2 periods = 166 time-indexed full dataset-source distributions**.

The source DATSET_KEY is a dataset UUID. An important technical point is that this concentration is for ALL KINGDOM PLANTAE SPECIES-rank indexed records, **not** only GRIIS-listed introduced taxa, and **not** deduplicated survey sites, observers or biological monitoring protocols. Each 50km cell can contain many records from the same observation campaign.

| Fixed country block | GRIIS source-eligible two-period 50km cells | Cell-median *late-period* share of the SINGLE largest dataset in all indexed plants |
|---|---:|---:|
| Canada Ontario | 10 | **0.943671** |
| US Midwest | 15 | **0.919509** |
| Australia Southeast | 16 | **0.543654** |
| UK Yorkshire | 12 | **0.534556** |
| Germany Saxony | 14 | **0.451257** |
| France Burgundy | 16 | **0.281344** |
| **Total** | **83** | **Do not pool raw records across countries** |

There is a major **source concentration** and potential source-coverage imbalance: US and Canadian 50km cells often have over nine tenths of all indexed plant records in a single dataset. The ecological recording process cannot be assumed comparable to France or Germany. A dataset name does not equal an independent recorder or a single ecological survey.

## 2. Actual v5.1 same-pair source-dataset-overlap sensitivity (K20)

[Real source dataset overlap and GRIIS plant association #42797029122](https://github.com/zuizui0223/tetsu/actions/runs/42797029122) **SUCCESS**, six strict synthetic source validations. For every source-qualified cell pair, compute an ALL-PLANT dataset composition similarity as the sum across common dataset UUIDs of the minimum of their two relative record shares: \(\sum_d\min(p_{id},p_{jd})\). Unknown/unattributed dataset source records are NOT treated as common across cells. This is a **lower bound on known dataset overlap**, not a cleaned introduced-only dataset similarity or an original event-level detection probability.

For 2018–2025, the *exact same eligible plant/road/rail pairs* were used before and after adding the dataset source overlap proxy to the rank residualization. Covariates previously included geographic distance, ROAD physical path, and the GBIF indexed recording count pair scale/balance.

| Same-source observational diagnostic | Six-country equal-weight median partial rank (negative is rail shortcut-hypothesis direction) | Country blocks in negative direction |
|---|---:|---:|
| Old source geographic / ROAD / GBIF record-count covariates | **+0.021259** | 3/6 |
| Plus all-Plantae DATASET_KEY compositional overlap | **+0.056600** | 2/6 |
| Plus overlap and maximum-dataset share imbalance | **+0.039974** | 2/6 |

US Midwest source K20 rho −0.264201 -> −0.263712 with dataset-overlap control; France +0.053809 -> +0.056248; Canada +0.201410 -> +0.191166; UK −0.011291 -> +0.056952; Germany −0.069825 -> −0.066717; Australia +0.411986 -> +0.345366. **No consistent hypothesis-compatible railway-specific floristic source effect emerges at K20**. Adding a source composition proxy is not a causal confound adjustment.

Reproduce from [v5.0 full DATASET source audit](../scripts/audit_gbif_dataset_composition_6country_v50.py), [v5.1 same-pair source comparative code](../scripts/analyze_dataset_overlap_adjusted_rail_plants_v51.py), [source tests](../tests/test_dataset_overlap_adjusted_rail_plants_v51.py).

## 3. Actual v5.2 JOINT K80-frozen cell cohort × dataset source comparison

[Real six-country fixed K80 + dataset-overlap analysis #42797571041](https://github.com/zuizui0223/tetsu/actions/runs/42797571041) **SUCCESS** with four synthetic source tests, all six pre-existing country blocks. The *same GRIIS-listed indexed source cells* must have at least 80 recorded introduced-plant source records in both historic periods; cells were selected by recording count alone, then rarefied at K20, K40 and K80. All-Plantae dataset recording overlap source covariate was computed on exactly those same pairs. Source 50km physical rail/road routes and both source-anchor definitions remain fixed.

| SAME K80-record-eligible original cells, recent 2018–2025 period | Previous source controls only: six-country median partial rank | Plus independent dataset overlap: six-country median partial rank | Direction-compatible country blocks **after dataset control** |
|---|---:|---:|---:|
| K20 | +0.021259 | **+0.056500** | 2/6 |
| K40 | −0.021333 | **+0.011045** | 3/6 |
| K80 | −0.064249 | **−0.069033** | **5/6** |

**Interpretation:** A hypothesis-compatible **exploratory sign pattern appears specifically at K80 with dataset overlap controlled**. It is not present at K20 or K40, and K80 was explicitly added **after observing earlier raw-vs-K20 sign discrepancies**. Thus this is a posthoc sensitivity-derived candidate, **not a preregistered replicated railway dispersal effect** or independent causal hypothesis test.

At K80, five country regions have negative descriptive partial ranks: US Midwest −0.193542, France Burgundy −0.058652, UK Yorkshire −0.020952, Germany Saxony −0.230200, Australia Southeast −0.079415; Canada Ontario is positive +0.140425. The six country coefficients are **not** six independent randomized train transport experiments, and their signs should not be given inferential significance without an independent study design.

The older period 2010–2017 under K80 fixed cohort and dataset source correction had country median +0.049069 with only 3/6 negative, further warning against generalized interdecadal transfer.

Reproduction: [joint fixed-record-cohort source script](../scripts/analyze_joint_rarefaction_dataset_overlaps_v52.py), [unit tests](../tests/test_joint_rarefaction_dataset_overlaps_v52.py), [real workflow](../.github/workflows/joint-recording-depth-datasets-v52.yml).

## 4. Out-of-analysis-region replication criterion frozen before new geography network outcomes

The [US Southeast independent out-of-analysis-region source protocol](holdout_us_southeast_v53_protocol.md) was committed on 2026-11-10 JST before its new physical RAIL/ROAD 50km source network geometry or dataset-source correlation had been evaluated. This block is pre-existing in the v2.9 map and source GRIIS v3.6, but did not contribute to original six-country ecological partial-rank coefficients. Its K80 main diagnostic is the *posthoc-chosen* source score; it is NOT a new independent country or a true preregistered transport mechanism test.

A geographically independent evaluation only happens if >=8 existing K80-eligible introduced-plant source cells and >=20 road/rail source-valid cell pairs survive. If not, the result is NOT ESTIMABLE and cannot be reselected or diluted into a positive p-value.

### Scientific status

- **Positive structural railway measurement remains real**: Overture rail resistor alternative-path electrical conductances and network shortest distances vary geographically and agree on synthetic networks with Julia Circuitscape.
- **Plant source status remains unresolved**: the K80 effect is chosen from a posthoc depth sensitivity, and the source dataset marginal is all plants rather than country-introduced-only data. GBIF legacy occurrenceStatus=PRESENT GET API gives zero or DETECTED HTTP400; current matched source DETECTED status facet counts and a bounded 7,200-row firstpage original-record quality audit are recorded in [v4.9 reproducibility report](GBIF_rarefaction_and_status_contract_v49.md).
- **Causal railway seed transport not demonstrated**: modern physical rail/road maps postdate historical source GBIF records and contain unknown active-service, station, freight, ports/urban/habitat and observer effort confounds.
- **Source analysis is still exploratory**: no pairwise p-values, no unplanned K selection to rescue an association, no generalization of selected-country median to global biogeography. For a true test, obtain a DOI-citable full original GBIF occurrence event sample, active train freight status, true ecoregion/landcover labels, then replicate in geographical regions chosen independently of this result.
