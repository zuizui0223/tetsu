# tetsu v5.3 — Frozen US Southeast out-of-analysis-region validation: NOT ESTIMABLE

**Source and analysis finalized November 10, 2026.** This note reports an unsuccessful *data-support* gate, not an ecological null hypothesis test. No railway-mediated seed dispersal is demonstrated.

## Locked before any new 50km network result

The [frozen US Southeast source protocol](holdout_us_southeast_v53_protocol.md) was committed before the new region's real Overture railway/road route or botanical similarity effect was computed. It selected the ORIGINAL v2.9 4×4 50km EPSG:6933 US Southeast geographical block at about longitude −78.60, latitude 35.80, outside the existing six-country source correlation region sample (the earlier US ecological block was Midwest). The main response and source controls were frozen from posthoc-selected K80:
- period 2018–2025, 80 indexed GRIIS-introduced GBIF source records per 50km cell in BOTH periods;
- country US GRIIS listed introduced speciesKey observed-only source composition;
- outcome = expected species intersection / expected species union under 80-record rarefaction;
- predictor log(physical RAIL source path length / ROAD source path length), using nearest source connector within each cell;
- rank source adjustment for geodesic separation, ROAD physical route, GBIF indexed source record balance/scale, and full all-Plantae GBIF datasetKey source composition overlap;
- minimum of **8 qualifying K80 botanical cells and 20 physically matched RAIL+ROAD cell pairs**, in addition to the original source record and 80% country-area screen;
- unqualified result must be **NOT ESTIMABLE**, not a zero effect and not a request to retune K for a desired sign.

## Real independent source acquisition: all three succeeded

[US Southeast source-only workflows #42798094566](https://github.com/zuizui0223/tetsu/actions/runs/42798094566), three parallel verified jobs SUCCESS:

| Data component | Actual source result | Evidence type |
|---|---:|---|
| Preexisting country-GRIIS K20 both-period plant source eligible cells | **12/16** | Frozen v3.6 GBIF GRIIS source opportunity, not actual plant population |
| October 2026 Overture physical RAIL eligible mapped segments | **4,713** | Source physical railway line geometry, not active train freight |
| RAIL physical route-valid 50km source cell pairs | **120/120** | Source physical topology, median augmented rail/geodesic path ratio **1.619842** |
| ROAD physical mapped source segments | **647,309** | Source physical roads, no legal traffic/vehicle restrictions |
| ROAD physical route-valid 50km source cell pairs | **120/120** | Source physical topology, median augmented road/geodesic path ratio **1.171015** |
| Full DATASET_KEY+OCCURRENCE_STATUS indexed source profiles | **12/12** source cells in BOTH periods | GRIIS-eligible botanical cells, no error |
| Top-one dataset share, late 2018–2025 cell median | **0.918335** | All-plants indexed record source concentration, not introduced-only |

**The K80 main support gate fails.**

[Actual original-source follow-up outcome-neutral gate #42798827089](https://github.com/zuizui0223/tetsu/actions/runs/42798827089), five unit tests passed and actual evaluation SUCCESS at completing the planned analysis:

| Precommitted K80 main diagnostic support | Observed |
|---|---:|
| Source-eligible GRIIS introduced flora 50km cells with BOTH periods >=80 RECORDS | **4 cells** |
| Physically RAIL+ROAD valid comparisons among these 4 cells | **6 dependent pairs** |
| Frozen minimum required | >=8 cells AND >=20 pairs |
| Primary K80 rank coefficient | **NOT ESTIMABLE (null)** |
| Positive/negative ecological association to report | **None** |

The same code can mathematically return correlations for less-restrictive K20 cohorts, but we explicitly did **not** relax the prerecord source gate after seeing K80 support. The 6 actual network-pair source observations do not become independent biological samples, and a trivial small-sample coefficient is not used to change the paper's sign.

### Scientific conclusion

The prior six-country K80 source sensitivity yielded 5/6 hypothesis-compatible partial rank directions only AFTER selecting K80 in a posthoc recording-depth analysis. It is a potentially interesting exploratory lead, but an original geographically out-of-analysis US Southeast source cannot provide K80 independent testing under the frozen source conditions. The five-of-six result has **not been externally validated**, and the source opportunity deficit is a material macro-replication blocker.

This negative DATA SUPPORT outcome is useful: a large total GBIF indexed occurrence source and a well-connected physical railway graph do not ensure that a country-GRIIS introduced-plant speciesKey survey has adequate effort at a strict equal-record depth. A valid next step is a source-independent, original DOI-backed botanical occurrence and recording-event dataset or new event-based sampling, followed by a separately frozen macro geography/biome coverage plan—not lowering K only for countries with desired source coefficients.

Reproducible:
- Original [locked protocol](holdout_us_southeast_v53_protocol.md)
- [Three actual source pipelines](../.github/workflows/us-southeast-holdout-sources-v53.yml)
- [Holdout field and support audit](../scripts/evaluate_us_southeast_geographic_holdout_v53.py)
- [Source support tests](../tests/test_us_southeast_geographic_holdout_v53.py)
- [Real support-gated holdout run](https://github.com/zuizui0223/tetsu/actions/runs/42798827089)

No current data source establishes the causal role of trains, freight, rail maintenance, trackside wind, flowering phenology, or heritable selection on plant traits.
