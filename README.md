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