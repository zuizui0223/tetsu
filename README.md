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