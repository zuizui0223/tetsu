# v1.9 — International rail graph × plant occurrence spatial-source audit

Current scope: two-country real Overture railway/road segments, international expansion scripts and physical-connector algorithm. **NO plant railway-dispersal biological effect is established.**

## Source provenance and first actual result

Overture Maps transport segments from **2026-10-21.0**, Amazon S3 global geospatial parquet, OpenStreetMap and Overture Maps contributors (ODbL). Classes counted in two preselected source windows: rail standard/narrow gauge/light rail and road motorway/trunk/primary/secondary/tertiary/residential/unclassified. These are **bbox-intersecting segment numbers**, NOT route kilometers, freight flows, active track or measured botanical exposure.

[Actual first Japan/Germany source run #42611543610](https://github.com/zuizui0223/tetsu/actions/runs/42611543610):
- JP Sendai: **429** eligible mapped rail segments, **64,100** mapped major-road segments.
- DE Leipzig: **1,176** eligible mapped rail segments, **15,174** mapped major-road segments.
- Initial GBIF species facet parsing failed, causing the overall joint-source run to be marked FAILURE. No joint ecological relationship was computed.
- GBIF facet source reports upper-case SPECIES_KEY, rather than always camel-case speciesKey; corrected source parser normalizes underscores/case. Source now requests observed PRESENT taxa explicitly. Corrected ten-country real-source run succeeded, with source details below.

**Correction and complete 10-country source result:** [actual successful Actions run #42612721887](https://github.com/zuizui0223/tetsu/actions/runs/42612721887), 2026-11-09. Six targeted unit tests passed, **10/10 source windows succeeded, 0 GBIF / Overture errors**, and geography-matched indexed source RECORD counts + source segment classes were returned:

| Source urban window | GBIF georeferenced 2018–25 Plantae species-rank PRESENT occurrence-record search count | Overture eligible-class rail segments | Overture selected major-road segments |
|---|---:|---:|---:|
| JP Sendai | 1,475 | 429 | 64,100 |
| DE Leipzig | 39,406 | 1,176 | 15,174 |
| KR Daejeon | 5,982 | 350 | 24,843 |
| TW Hsinchu | 40,780 | 333 | 16,646 |
| US Madison | 82,299 | 111 | 23,408 |
| CA Ottawa | 82,746 | 291 | 29,432 |
| GB York | 39,030 | 220 | 10,042 |
| FR Dijon | 37,233 | 745 | 15,509 |
| AU Canberra | 169,109 | 31 | 38,097 |
| NZ Hamilton | 12,757 | 111 | 8,629 |
| **Ten source windows total** | **510,817** | **3,797** | **245,880** |

The 510,817 items are raw matched **GBIF record counts**, not downloaded rows or distinct species/observed individuals. Only top 50 GBIF speciesKey facet categories were requested per region; **this top50 subset cannot serve as complete plant assemblage or a verified alien plant status dataset**. Different urban windows have different land areas, efforts and segment splitting. Do **not** interpret record or segment sums as species richness or infrastructure density, and do not claim an ecological association from a joined source availability table.

Original fixed geographic study windows before seeing taxon outcomes: Sendai JP, Leipzig DE, Daejeon KR, Hsinchu TW, Madison US, Ottawa CA, York GB, Dijon FR, Canberra AU, Hamilton NZ. These 10 urban samples are feasibility frames, not national ecological replication samples.

Source scripts: [v1.9 Overture × GBIF geographic audit](../scripts/audit_overture_gbif_overlap_v19.py), [six tests](../tests/test_overture_gbif_overlap_v19.py), [2-country workflow](../.github/workflows/overture-gbif-overlap-v19.yml), [10-country workflow](../.github/workflows/overture-gbif-tenregion-v19.yml).

## Actual topology code is not just distance-to-track GIS

New [rail graph utilities](../src/tetsu/rail_graph.py) and [seven synthetic tests](../tests/test_rail_graph_macro.py):
- Build physical railway connector graph from Overture segment ID and connector_id / at linear references, not from crossings in an image.
- Distinguish graph connected/disconnected, compute shortest connector path in meters, count components, prevent duplicate/invalid source IDs, wrong class inclusion and zero-length edges.
- Modeled edge lengths currently assume normalized linear reference scales with source segment length; source geometry calibration is needed before numerical ecological models.
- Window-edge graph truncation, tunnel vs at-grade, active/disused status and freight traffic remain UNKNOWN. Source data only gives mapped geometry, not propagule flux.

## Full macro hypothesis and strict evidence gate

**Across independently sampled countries/ecoregions, do railway-connected pairs of plant assemblages share more non-native taxa than road-connected pairs matched in Euclidean distance, climate, ports, urbanization, vegetation cover and record sampling intensity?**

For genuine scientific inference, first obtain DOI-citable source GBIF occurrence rows (not country count totals or top-50 species facets). Validate accepted vascular plant taxa, regional country-specific alien/native status, coordinate uncertainty, dataset/observer and effort, duplicates, and cell-level 25/50km equal-area geographic placement. Next compute **rail and road network connectivity** from explicit connectors and validate enough road-matched pairs in independent geographic blocks. Use block holdouts/stratified network nulls, not independent-point pseudoreplication. Only then estimate a railway coefficient and its uncertainty.

A top-50 GBIF facet is not a complete plant assemblage. An indexed record is not a living individual or biological absence when unobserved. A mapped railway is not an active freight train. **No train-mediated dispersal, plant adaptation or flowering effect was confirmed in v1.9.**

Official references:
- https://docs.overturemaps.org/guides/transportation/segments-and-connectors/
- https://techdocs.gbif.org/en/openapi/v1/occurrence
- https://techdocs.gbif.org/en/data-use/api-downloads
