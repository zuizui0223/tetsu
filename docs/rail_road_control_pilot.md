# v0.5 fixed rail-vs-road observation contrast (2026-11-07)

**Exploratory observational test; not ecological railway adaptation, seed transport, train slipstream or mowing effects.** This document freezes threshold decisions before inspecting outcomes from the [run](https://github.com/zuizui0223/tetsu/actions/workflows/tokyo-rail-road-pilot.yml).

## What changed from the one-axis v0.4 pilot

Added official **N13, 2024 road-centreline** spatial data, published by MLIT in April 2026; Tokyo first-order 5339 mesh contains 1,943,251 centreline features (verified archive; [N13 source audit](https://github.com/zuizui0223/tetsu/actions/workflows/n13-road-audit.yml)). Roads are NOT routes with actual measured road traffic. The attributes include type/width and map registration date (not opening date). Source is CC BY 4.0 and has reproduction/attribution conditions: https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N13-2024.html.

Rail from official N02 2025 is **passenger railway geometry**, omitting sidings/freight-only lines, and does NOT quantify rail traffic. Source: https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N02-2025.html.

## Source query frame and expectations

- **First location:** Tokyo v0.4 locked bounding box: lat 35.52–35.83, lon 139.56–139.92. No moving the polygon to maximize a result.
- **Years:** 2023–2025 only; map vintage approximates observation vintage.
- **Focal herb taxa (exploratory):** *Trifolium repens*, *Solidago altissima*, *Plantago lanceolata*, *Erigeron annuus*. No inferred introduction-pathway status or species identity outside these exact names.
- **Photographic target-group reference:** other Plantae iNaturalist research-grade observations with photographs (taxon ID 47126), in the same bbox and time interval. Target species are excluded from this reference to avoid denominator duplication.
- **Sampling:** FULL search for each focal taxon, but reference is ONE deterministic hash-selected iNaturalist observation ID page (up to 200) per each of the 36 year-month strata. The reference sampling is a one-page-per-month ID-cluster design, not independent record-level random sampling. Each record on the selected page has a *page inclusion probability* 1/number_of_pages; the correct inverse probability weight is number_of_pages, including when the chosen final page is short. Users/photos within one page may be correlated and residual selection bias remains. Source totals and page IDs (not individual observation IDs) are stored in aggregate reports. The first erroneous weight total/page_length (run #42376261914) was repaired in #42376881325.
- **Coordinate gate:** un-obscured, not marked cultivated/captive, dated, photographed, numerical spatial accuracy 1–100m; no imputation of unobserved plants. Original photo image pixels are not used or redistributed.

## Predefined contrasts

**Primary exploratory overlap:** among recorded plants **≤100 m from a mapped road**, compare observations **≤250 m from mapped rail** with **≥1,000 m** from rail. Exclude intermediate rail distances. Focus-taxon observations relative to other plant photographs are summarized as:

```
relative photo-record association
  = (target near-rail / reference near-rail)
    / (target far-rail  / reference far-rail)
```

For sparse tables, a 0.5 continuity correction is explicit. Do not interpret this ratio as biological occupancy or a local phenology response. Also report a calendar-month Mantel–Haenszel contrast (descriptive only).

**Fixed sensitivities:** rail ≤100 m and road ≤100 m, and rail ≤250 m / road ≤150 m. Use the same ≥1,000 m rail background. No threshold is retuned to change the result direction.

**Coverage gates:** (1) number of quality-passing target and other-plant observations; (2) distinct observers in each; (3) 2-km UTM cells containing focal/reference records with both rail-near and rail-far settings; (4) fraction excluded by high spatial uncertainty; (5) month-wise pages sampled and counts. If replication or source balance is poor, no ecological effect statement.

## Explicit alternative explanations and limitations

- Railways and roads co-occur with cities, bridges, gardens and parks. Road-centreline distance alone does NOT control urban land cover, public-access restrictions or ecology.
- iNaturalist plants are observer-selected and often popular flowering subjects; neither a photo's absence nor a missing annotation means nonflowering or a plant absence.
- Geographic coordinates may reflect park entrances and user GPS errors. 100m accuracy is far too coarse to prove plants occupy railway land.
- A positive association would be **rail-related photographic occurrence availability**, not train-borne dispersal. A null result with these proxies cannot rule it out.
- Observations and geographic cells are correlated: **no conventional independent-record p-values** are emitted. Spatial block/observer validation and a better reference sampling design are prerequisites for publication.
- The 2024 N13 road layer has ~25m reported position uncertainty for the layer. Rail 100m thresholds should be treated as sensitivity, not exact track microhabitats.
- Rail network topology, cargo volume, mowing timing, air-flow exposure, and common-garden genetic effects remain separate future workstreams.

## Scope of replicability

Scripts:
- `scripts/audit_n13_road.py` + `scripts/audit_rail_vs_road_inat.py`
- `src/tetsu/transport_overlap.py`
- Workflows: `n13-road-audit.yml`, `tokyo-rail-road-pilot.yml`

All downloads are API version/time dependent. Reports include time, source SHA-256, per-month counts, filters and failure diagnostics. Raw iNaturalist coordinates and IDs are never included in CI artifacts.

## Next improvement (before any hypothesis claim)

Sample another independently selected geographic region, replicate focal taxa or a properly justified trait class, include a target-group observer-effort offset, independently validated vegetation/land-use data, physically measured wind/grass-cutting events and standardized repeated plant photographs. Still distinguish environmental plasticity from selection of species and heritable adaptation.
