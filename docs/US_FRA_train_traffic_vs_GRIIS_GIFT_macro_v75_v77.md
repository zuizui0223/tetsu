# tetsu v7.5–v7.7 — railway TRAFFIC observation versus functional plant source patterns

**Verified genuine original connected-source runs, November 11, 2026.**
This is an exploratory macro-source comparison, NOT observed propagule transport,
railway-side mowing, mechanical gust effects or heritable local adaptation.

## Real traffic data: official FRA-crossing-form fields via publicly hosted ArcGIS mirror

- Predefined geographical source: 3 US blocks (US Midwest, US Southeast,
  US Pacific) from v2.9's original 24-block, 8-country 50km EPSG6933
  geography. Each contains 16×2500 km² cells.
- [Corrected actual FRA train census #43048727413](https://github.com/zuizui0223/tetsu/actions/runs/43048727413)
  and [artifact](https://github.com/zuizui0223/tetsu/actions/runs/43048727413/artifacts/13867289763):
  **9,564** real crossing inventory entries, validated paginated public
  service and source geometry. **3,121 closed** crossing rows explicitly
  excluded from operating-intensity summaries. Of **6,443 confirmed
  currently open** crossings, **5,211** have both day+night through-train
  counts and a reported train-count source year 2018–2026. Across
  48 frozen geographical cells, **42** contain some such source counts.
- This source provides *through trains at road crossings*, **not freight
  trains by number, annual train-kilometers, cargo tonnage, historic rail
  flow, actual seed/carriage dispersal, or train-induced wind speed**.
  Road crossings are strongly conditioned on local roads/urbanization.
  An OPEN crossing with reported through trains=0 is not an observed
  active railway segment. The mirror is not independently certified
  fully synchronized with all official FRA updates.

| Source block | Source crossing rows | Open crossing rows | Open with full year and 2 through-train fields |
|---|---:|---:|---:|
| US Midwest | 3,809 | 2,557 | 2,314 |
| US Southeast | 3,410 | 2,418 | 1,859 |
| US Pacific | 2,345 | 1,468 | 1,038 |
| Total | 9,564 | 6,443 | 5,211 |

## Real macro data integration: two original US GRIIS/GIFT source blocks (v7.6)

[Full real run #43050469923](https://github.com/zuizui0223/tetsu/actions/runs/43050469923)
and stricter [jackknife + nonzero-train-only sensitivity #43051783705](https://github.com/zuizui0223/tetsu/actions/runs/43051783705),
[artifact](https://github.com/zuizui0223/tetsu/actions/runs/43051783705/artifacts/13870415618).

Freeze US Midwest (14 eligible cells) and US Southeast (13 eligible
cells). 2018–25 GRIIS country-listed introduced GBIF species-key
source, pre-frozen 50km cell IDs, same independently obtained GIFT
1.2.1 growth forms and FRA *reported median total daily through trains
at currently open crossings*. Primary descriptive residualized rank
correlation controls log indexed all-Plantae GBIF records and log count
of open FRA crossing inventory rows per cell. Source endpoints are
not field botany plots; no causal model or pairwise pseudo-p-value.

| Trait source outcome | US Midwest rho (14 cells) | US Southeast rho (13 cells) |
|---|---:|---:|
| Raw recorded introduced species richness, all forms | −0.654119 | +0.118026 |
| K20 expected species richness, all introduced forms | +0.062601 | +0.304973 |
| Raw recorded herbaceous richness | −0.654119 | +0.156358 |
| K20 expected herbaceous richness | −0.082780 | +0.229904 |
| Raw recorded woody richness | −0.559832 | +0.012532 |
| **K20 expected woody richness** | **−0.377168** | **−0.381404** |

The raw/rarefied source responses can differ substantially. K20
expected species richness is derived by hypergeometric *indexed record*
rarefaction (not true field or detected species occupancy); source GBIF
effort control cannot remove nonrandom botanical provider, habitat
and spatial sampling bias.

### Leave-one-50km-cell-out sensitivity

On the original 14/13-cell source cohorts including reported zero-train
medians, K20 woody sign remains negative across **14/14** Midwest and
**13/13** Southeast single-cell exclusions:
- Midwest min/max jackknife rho **−0.817074 to −0.225759**.
- Southeast **−0.549301 to −0.237412**.
This is a *within-block observation sensitivity*, NOT an independent
country replication, confidence interval or hypothesis p-value.

### Crucial positive-train-only boundary

The original eligible source cohort includes **one** Midwest and **five**
Southeast cells whose OPEN-crossing median reported through trains = **0**.
Predeclared source sensitivity retains only cells with *positive*
reported through-train medians, while retaining the original
**>=10 unique 50km cells per geographic block** support rule.

| Country block | Cells with positive reported through-train median | Woody K20 partial rho, GBIF+open-crossing adjusted |
|---|---:|---:|
| US Midwest | **13** | **−0.225759** |
| US Southeast | **8** | **NOT ESTIMABLE** (<10 source cells) |

Thus **a repeated woody richness–train-intensity association in
TWO independently supported active-train macro blocks is NOT established.**
At most the exploratory negative sign remains in one selected Midwest
source block after excluding all zero-through-train-median cells.
No positive train intensity gradient can be asserted for Southeast
from its 8 positive-through-train source cells.

Current train-count report years span 2018–2026 whereas plant records
span 2018–2025; the report year and snapshot cannot establish
contemporaneous historical exposure. Freight-specific train movement
and actual stress remain unmeasured. Adding rough spatial x/y controls
does not solve urban/climatic/landcover confounding with n=13/14.

## Independent predefined Pacific source opportunity: v7.7

[REAL SUCCESS #43051476174](https://github.com/zuizui0223/tetsu/actions/runs/43051476174)
and [source artifact #13869297497](https://github.com/zuizui0223/tetsu/actions/runs/43051476174/artifacts/13869297497).

US Pacific was defined as one of the original v2.9 24 fixed geographic
200km source blocks, **before reading v7.6 signs**. GRIIS source
opportunity was assessed using the current GBIF indexed DETECTED and
NOT_DETECTED status facets (NOT unsupported GET/PRESENT legacy query),
all speciesKey facet pages, independent Natural Earth source-country
polygon fraction ≥80%, both 2010–17 and 2018–25 positive reporting
minimums, and published US country GRIIS list.

- **32/32** source cell×period facet requests in the 16 pre-frozen
  US Pacific cells completed, zero source errors.
- **8/16** original 50km Pacific cells pass both-year GBIF source
  record, accepted country-GRIIS taxon/record, and quality gates.
- Those eligible cells collectively contain **566** distinct indexed
  country-GRIIS speciesKey categories observed in the 2018–25 source.
  This is not 566 directly surveyed species, and some indexed records
  may include source statuses other than true DETECTED at species level.
- Intersecting those 8 GBIF-eligible Pacific cells with independent
  current FRA open-crossing cells with ≥5 year/day/night source
  observations leaves **7** unique 50km cells. This fails the same
  ≥10-cell descriptive macro eligibility rule **BEFORE looking at any
  Pacific train–plant coefficient**. No Pacific association is fitted.
- Pacific source used a CURRENT GBIF indexed facet method whereas
  the original 2-region legacy flora source used the earlier
  historical indexed methods; mixed-vintage results are not pooled.

## Ecological interpretation

**Supported source-level progress:** the physical railway map can now
be supplemented with genuine crossing-level through-train reports,
and 27 50km cells across two prespecified US regions support the
initial GIFT functional composition comparison. A negative
woody K20–train-intensity source correlation appears in those two
initial samples and is robust to any ONE cell's omission.

**Not supported:** two-region active train-traffic effects survive
strict nonzero-train support; a third Pacific independent confirmatory
effect is estimable; train-borne seed dispersal; freight-specific
transport; causal selection by mowing, gusts or ballast heat; adaptive
evolution, trait plasticity, or an overall ecological railway effect.

**Next actual identification step:** use independently dated freight
operations (service and measured freight movement or tonnage, not
crossing counts) and true rail-verge botanical occurrence/management
profiles, then expand prespecified country×ecoregion replication.
Record actual mowing/weed-control events, airflow measurements and
species regeneration/seed traits; route exposure + vegetation trait
composition alone does not identify selection versus arrival.
