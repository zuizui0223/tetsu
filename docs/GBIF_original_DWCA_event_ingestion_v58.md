# tetsu v5.8 — Original GBIF DWCA event ingestion is implemented, but authentic occurrence downloads have NOT been obtained

**Verified November 10, 2026 JST.** This provides a safe, reproducible importer for the next source-quality stage of the 50km multi-country railway and country-GRIIS introduced flora macro project.

## Status of real ecological evidence

Real GBIF indexed source facets are available for 83 50km country GRIIS introduced-plant cells in six country regions, plus 2026 source physical railway/road shortest paths and Circuitscape-equivalent electrical network redundancy. The latest real [dominant-dataset leaveout test #42832865840](https://github.com/zuizui0223/tetsu/actions/runs/42832865840) shows a posthoc K80 5/6 railway-compatible source sign does NOT survive removing the dominant GBIF dataset: only three country regions retain K80 sufficient original source opportunity; on exactly matched cell pairs France changes sign (−0.059143 vs +0.101852), and at K20 the US Midwest source partial-rank almost vanishes (−0.263712 vs −0.012227).

**This leaves the biological railway seed-dispersal hypothesis unconfirmed.** It does not mean actual train seed delivery is absent; the source indexed record process is insufficiently independent.

## GBIF official original ZIP structure

GBIF official [DWCA download format documentation](https://techdocs.gbif.org/en/data-use/download-formats) specifies a compressed archive containing **interpreted occurrence.txt**, **verbatim.txt**, **meta.xml**, **metadata.xml**, and relevant citation/rights files. The interpreted occurrence.txt fields include GBIF record identity, datasetKey, eventDate, speciesKey, coordinate and coordinateUncertaintyInMeters, and taxon/country/status metadata.

The interpreted DWCA record-level occurrenceStatus generally uses PRESENT / ABSENT even though current GBIF REST indexed OCCURRENCE_STATUS facets report DETECTED / NOT_DETECTED. The [documented v4.9 API incompatibility](GBIF_rarefaction_and_status_contract_v49.md) must not be silently concealed or solved with an invalid REST filter. The importer explicitly checks actual row-level status rather than assuming that downloaded event rows are positive biological observations.

## Implementation: streaming, privacy-safe original archive QA

- [New GBIF original DWCA streaming importer](../scripts/ingest_GBIF_DWCA_original_events_v58.py): requires user-provided original GBIF DWCA ZIP, corresponding SUCCEEDED download metadata with download key/DOI and official HTTPS host, the preexisting country GRIIS accepted speciesKey registry, and one of ten original geographic source blocks. It does **not** request downloads or perform background API work.
- [Source behavioral tests](../tests/test_ingest_GBIF_DWCA_original_events_v58.py) and [offline-only GitHub workflow](../.github/workflows/GBIF-DWCA-original-ingest-QA-v58.yml): actual CI SUCCESS on FOUR synthetic test checks. The tests use a temporary constructed ZIP and a synthetic DOI-shaped string explicitly flagged as test-only. **No real GBIF DOI event file was downloaded or analyzed in this workflow.**
- Stream imported source interpreted TSV without reading raw recorder identities into retained output. Check original GBIF ID duplicates, duplicate dataset+occurrenceID, same dataset+event+species for repeated event observations, exact country, source year and matching original 50km WGS84 geometry.
- Annotate ONLY accepted country-GRIIS speciesKey categories (non-members remain status UNKNOWN; never automatically native).
- Output per 50km equal-area source cell and 2010–17 / 2018–25 period the number of original source events per GRIIS plant speciesKey under all located coordinates, strict <=5,000 m stated coordinate uncertainty and sensitivity <=25,000 m. Missing/zero uncertainty is NOT assumed precise.
- Discard non-PRESENT/DETECTED event statuses before deriving botanical source occurrence. Store only aggregate counts and taxon keys, no individual GBIF IDs, recorder names, exact observation coordinates or raw event files. Temporary SQLite event-key hashes are removed at process exit.
- Input source archive SHA256 and download metadata identity are recorded. Source metadata provided by the user is marked **not independently authenticated against a live GBIF API**; official DOI citation must be verified via GBIF download information rather than accepting arbitrary metadata.

These are stricter QA controls than indexed speciesKey facets or first-page GBIF search. They still do not make an eventID a standardized survey or infer physical train vector transmission.

## How a user-supplied real archive would be processed

This is a prepared workflow, NOT an instruction to submit credentials into a public code repository. After a GBIF account holder requests and completes an official authenticated source download, obtains the ZIP and official GBIF download metadata JSON and reports the source DOI:

    python scripts/ingest_GBIF_DWCA_original_events_v58.py \
      --block us_midwest \
      --dwca-zip /private/path/to/completed-gbif-original.zip \
      --gbif-download-metadata /private/path/to/official-download-metadata.json \
      --griis-accepted-keys /private/path/to/griis_9country_alien_backbone_keys_v24.json \
      --output results/audit/us_midwest_verified_event_source_v58.json

The source original download should use one of the ten already prepared [official download request manifests](../scripts/prepare_GBIF_original_DWCA_download_manifests_v55.py) in [real manifest generator run #42800153221](https://github.com/zuizui0223/tetsu/actions/runs/42800153221). Those manifests were PREPARED but NOT SUBMITTED; no actual DOI has been obtained by tetsu in this workflow. An authenticated GBIF account is necessary for a new official download request. Do not commit individual record archives, passwords, API credentials, observer identities or precise geolocation to GitHub.

## Limits for downstream macro botanical inference

A DOI-backed original GBIF archive and row-level metadata would allow species observation events and independent data sources to be audited much more carefully than indexed aggregation alone, but would STILL NOT prove rail-mediated seed dispersal without historical active freight lines and actual seed transfer evidence. Country/biome, climate, urbanization, rail verge and mowing intensity remain open ecological confounders. Regional summary is the independent macro unit, not 120 within-region overlapping plant cell-pairs.

**Scientific decision:** the demonstrated positive results so far concern mapped infrastructure geometry and recording-source data quality, NOT a causal positive plant introduction effect. Continue with source event quality and historically verified rail services before presenting a biological macro causality claim.
