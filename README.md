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