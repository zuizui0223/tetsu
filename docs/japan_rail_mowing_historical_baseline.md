# Railway maintenance is a separate selection gate — documented Japanese baseline

**Publication-backed historical evidence, NOT a new `tetsu` experimental effect estimate.**

## Japan-specific vegetation management contrast

Ito, Ueki & Sakamoto (1982), *Studies on the Total Vegetation Control in Railroad I. Major Weeds and Factors Affecting Their Distributions*, Journal of Weed Science and Technology 27:41–48, DOI [10.3719/weed.27.41](https://doi.org/10.3719/weed.27.41). Original study: national questionnaire to **37 Japanese National Railways maintenance depots**, plus plant counts along approximately **370 km** of railway near Osaka observed from trains.

The abstract reports several distinct *qualitative associations*, not per-species randomized treatment effects:

| Functional comparison | Plant reported | Author-reported direction | Caveat |
|---|---|---|---|
| Herbicide-managed versus mown areas | *Setaria* spp. | More common in herbicide-managed areas; markedly less common under mowing | Herbicide/mowing and site type not randomized |
| Herbicide-managed areas | *Solidago altissima* | Tended to be less common | No experimental causal estimate |
| Herbicide-managed areas | *Imperata cylindrica* | Tended to be less common | No experimental causal estimate |
| Greater surrounding urban development | *Miscanthus sinensis*, *Imperata cylindrica* | Tended to decline | Urban conditions/observer and management confound |
| Greater surrounding urban development | *Solidago altissima*, *Pueraria lobata* | Tended to increase | Not evidence of railway-mediated arrival |
| Tall railway embankments | *Pueraria lobata* | More abundant | Habitat morphology / accessibility is a likely confound |

**Reinterpretation for tetsu:** herbicide management, cutting, exposed dry track ballast, embankment geometry and train-generated gusts are different exposures. A railway distance proxy measures none of these. Therefore split the local plant selection test by measured **maintenance regime** before attributing low cover or late flowering to train wind or mowing.

## Fundamental problem in a single rail/photo coefficient

The photographed/open-flower pattern depends on several sequential filters:

```
transported propagule
    -> lands in habitat
    -> survives/establishes despite management
    -> flowers at the observation time
    -> photographer selects the plant/flower
    -> image is georeferenced and licensed
```

As a formal factorization, `P(observed flower)` is a product of conditional probabilities along the chain. **A rail network can simultaneously raise arrival probability while cutting/herbicide suppresses establishment or the fraction photographed flowering.** Low near-rail iNaturalist photo proportions can therefore not be interpreted as a failed railway seed-dispersal mechanism. This is a logical mechanism distinction, not evidence that this opposite-direction mechanism happened in Japan.

## Field evidence required for an actual modern rail-mowing result

- independently documented site-level `cutting` vs `herbicide` vs `nonmanaged` (incl. number/date/cut height and vegetation-clearing history), not visual guesses from plant photographs;
- 10+ independently replicated corridors and matched roadside/background locations, within taxon and broadly comparable climate/urban strata (planning only, formal power analysis required);
- fixed-angle/quadrat image plant vegetation cover, photo flower pixel fraction and flowering/nonflowering **repeat visits** around each verified management event, retaining unflowered visits as observations;
- physical wind exposure: maximum instantaneous gust/ambient wind, track distance, obstacles, rail traffic schedule, train type, direction and safe authorized measurement; not inferable from N02 centrelines;
- repeated seeds/genotypes or common gardens to separate **environmentally plastic regrowth** from heritable evolution and between-species filtering.

## How to connect to the railway movement analysis

The eventual core model should have independently measured, separately identifiable components:

1. **Transport**: network connection and freight/port activity -> propagule arrival / geographic sharing;
2. **Establishment and management**: documented cut/herbicide, substrate/microclimate -> whether deposited plants survive;
3. **Reproduction**: observed flowering probability over DOY/GDD and fixed-photo display after maintenance;
4. **Observation**: plot visits vs opportunistic plant photos; photo-licensing and accessibility are selection processes.

These are multiple falsifiable claims, **not a single 'rail effect' regression**. A positive moving-seed result and a suppressed rail-side flower-display result could both be true without logical conflict. No new field outcome for these mechanisms is available in `tetsu` yet.

References:
- Historic Japanese rail vegetation survey (1982): https://www.jstage.jst.go.jp/article/weed1962/27/1/27_1_41/_article/-char/en
- MLIT operational guidance on plant removal along and outside rail boundaries: https://www.mlit.go.jp/tetudo/tetudo_tk7_000025.html
- Existing tetsu design: [protocol](protocol.md), [photo/phenophase plan](rail_environment_traits_protocol.md), [fixed plot cutting response](mowing_roi_event_study.md).
