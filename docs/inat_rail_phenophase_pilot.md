# iNaturalist railway phenophase feasibility pilot (v0.4)

## Motivation

Evaluate whether there is **enough spatially precise, taxon-validated, photo-annotated evidence** to justify a within-species rail-versus-background flowering analysis in Japan. This is **not** the flowering shift, plant plasticity, or evolutionary adaptation test itself.

Important iNaturalist annotation IDs (documented by iNaturalist):
- `term_id=12`, Flowers and Fruits.
- `term_value_id=13`: ≥1 open flower visible; **not** a count of flower abundance or flowering onset.
- `term_value_id=21`: no flowers, fruits, or seeds **visible in the images**; this may be an incomplete view or a nonflowering plant.
- `15`: flower buds; `14`: fruit/seed (not requested in this initial focused pilot).

Authoritative definitions: https://help.inaturalist.org/en/support/solutions/articles/151000191830

## Immutable first-pilot geographic definitions

Bounding boxes (SW latitude, SW longitude, NE latitude, NE longitude):
- Tokyo: `35.52,139.56,35.83,139.92`.
- Osaka: `34.53,135.35,34.83,135.76`.
- Sendai: `38.16,140.72,38.40,141.04`.

These are deliberately **not administrative city boundaries**. They were fixed before viewing API result counts. Future expansion to multiple random geographic blocks needs a different protocol and version.

Focus species:
- Japanese nationally documented accidental-only alien-introduction labels: `Briza minor`, `Lepidium virginicum`, `Rumex crispus`, `Sonchus asper` (based on Egawa & Koyama v1).
- Other common urban plant *exploratory species only* without introducing unverified pathway labels: `Trifolium repens`, `Plantago lanceolata`, `Solidago altissima`, `Erigeron annuus`.

Some taxa may be scarcely recorded or differ between urban areas; **lack of records is a feasibility finding, never an ecological absence**.

## Script and real raw-source queries

Run:
```bash
pip install requests pyshp 'shapely>=2,<3' pyproj
python scripts/audit_inat_rail_phenophases.py --output results/inat/rail_phenophase_audit.json --limit 100 --delay 1.10
```

The source calls use https://api.inaturalist.org/v1/observations , bounded by exact scientific taxon name, research grade, photo present, 2018-01-01 to 2025-12-31 and each locked geographic bounding box. A count query records all observations for the taxon (many **lack** relevant annotations). Each annotated category is fetched once, with the first 100 IDs in ascending order, while recording the API-reported total. This is deliberately **capped and NOT a random sample**: where `sample_truncated=true`, do not compare proportional abundance in returned records to the unqueried records.

iNaturalist advises approximately one API request per second and observation exports/GBIF datasets for bulk analysis:
https://www.inaturalist.org/pages/api+recommended+practices

GIS: official MLIT 2025 N02 passenger rail shapefile https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N02-2025.html, fetched with a checksum per run. Separate local UTM transforms used for each city. Use nearest polyline distance in metres as a **descriptive rail proximity marker**, not rail passenger activity, rail freight traffic, train-induced wind stress or herbicide/cutting intensity. N02 excludes freight-only line geometry and sidings.

## Input quality controls

A photograph observation contributes to a spatial proximity count ONLY if:
- un-obscured and not geoprivacy-private;
- not explicitly captive/cultivated;
- numerical coordinate accuracy >0 and <=100 m;
- valid lat/lon and observation date;
- photo present.

API totals are recorded separately from **coordinate-quality filtered records**. The output only exports taxon × area × annotation × railway proximity bins (0–250 m, 250–1,000 m, >1,000 m), group sizes, distinct observer counts, API truncation indicators, rejection diagnostics, counts with at least one CC0 / CC-BY / CC-BY-SA / CC-BY-NC or CC-BY-NC-SA licensed photo. No photo is downloaded, and no raw coordinates or observation IDs are uploaded as output.

Note that the classification 0–250m and >1km is *about georeferenced observation proximity to a mapped rail line*, not proof that the plant was growing on rail land. An accuracy threshold of 100m is still wide relative to rail verges, so this is not sufficient for microhabitat hypotheses.

## Feasibility gate

- **PASS for further measurement design** requires actual quality-filtered data in rail-near and well-separated background strata, multiple independent observers and sites, non-obscured locations and an adequate number of full-plant CC-licensed photographs.
- **FAIL for railway flowering shift inference** if samples are sparse, truncated, dominated by single observers or if apparent near-rail plants may just be station gardens/parks. An API's `Flowers` annotation is observational evidence of a flower in a photo but does NOT estimate flowering *onset*.
- **True cutting effect:** fixed individual or quadrat images before and after independently verified mowing. Opportunistic iNaturalist cannot supply reliable unobserved cutting history.
- **Wind and adaptation:** measured peak gust/airflow, shielding, site light/temperature/ballast, and ideally common-garden progeny/reciprocal-transplant needed. Rail distance should not be used as a proxy for heritable wind adaptation.
- **Generalization:** use independently replicated corridors rather than treating photos or taxon-by-city samples as replicates of rail systems.

## Status

Source-code / data-access gate: evidence to be added after GitHub Actions reports real values. Do not silently reuse prior artificial photo masks as field observations.
