# v0.7 — Actual 2024 NPN plant phenology activity at registered rail-near sites

**Updated 2026-11-08 JST.** Research-quality-gate and coverage study, **NOT evidence that trains affect seed dispersal, flowering, grass cutting, or evolution**.

## Frozen question and purpose

The previous v0.6 report identified 105 Maryland (MD) and 65 Arizona (AZ) *registered* USA-NPN sites within 250 m of a 2026 FRA railway line. **Registration is not evidence of any observed plant in 2024.** The scientific next question is whether real repeated Open flowers status observations exist for **the same taxon** in both rail-near and independent rail-far sites. Determine this **before** fitting any railway flowering effects.

Source-backed previous registrations: [#42379733758](https://github.com/zuizui0223/tetsu/actions/runs/42379733758). Source: USA-NPN `https://services.usanpn.org/npn_portal/stations/getAllStations.json` and FRA railway geometries `https://services.arcgis.com/xOi1kZaI0eWDREZv/ArcGIS/rest/services/NTAD_North_American_Rail_Network_Lines/FeatureServer/0`.

## Data sources and sampling frame

- States **MD and AZ**, initial year **2024**, fixed before any new outcomes.
- Near: a registered NPN site centre within **0–250 m** of mapped FRA line.
- Far: a registered site centre **≥1000 m** from FRA rail.
- Full registered near-site census up to **120 candidate sites per state** (covers the previously reported 105 MD and 65 AZ, if the upstream registry remains stable).
- Choose a unique far registered station closest to each near centre in the same state, using a projected UTM coordinate distance and requiring **≤75 km** between paired centres. Match geographically before reading ANY taxon status. Retain unmatched near sites in the **census** even when they cannot enter paired contrasts.
- Fixed SHA256 ordering of near site IDs; no selection according to taxonomic richness, flower phenology, positive/negative result or presence.
- Far site candidate pool is only the registered stations with valid coordinates. Independent urban land cover, roads, station access, wind or cutting history are **not yet controlled**. Geographic matching alone is NOT sufficient to infer a treatment effect.
- FRA network data are from **2026**, two years after the observation year; mapped lines need historical operation verification. Distances relate to **station centroids**, not the registered plants' true microhabitat.

The USA-NPN R client documents `station_ids`, `species_ids`, `phenophase_ids` and `states` filters for status data: https://usa-npn.github.io/rnpn/reference/npn_download_status_data.html. Requests use the NPN `observations/getObservations.ndjson` endpoint and explicitly check that every returned record has the requested `site_id` and observation year. Failures are reported and **never turned into no records**.

## Eligibility definitions set before the full census

- **Registered**: listed in NPN station registry.
- **Site actually active in 2024**: NPN returned ≥1 raw status observation for that site in 2024.
- **Site with Open flowers observation**: ≥1 2024 raw NPN row with `phenophase_description == "Open flowers"`.
- **Individual/year eligible for phenophase analysis**: same site, species, individual ID and year; ≥4 distinct dates with verified phase status in {0,1}; both ≥1 observed YES and ≥1 observed NO. Status -1 is unknown, not a NO. Duplicate/conflicting same-day states are deduplicated or rejected.
- **Site×species eligible**: has ≥1 individual with the above 4-date and observed YES/NO eligibility.
- **Common species across rail exposure classes**: at least one eligible near and at least one eligible far site **of the same NPN species_id**. Counts per taxon and independent site are reported, **not pooled photos as independent replicates**.
- **Prospective site replication gate for an eventual rail-near versus far biological comparison**: for each taxon ≥10 **independent verified active** rail-near stations and ≥10 far controls, at least 3 independent matched geographical/landcover blocks with overlapping treatments, observations in multiple years, and independently measured temperature/GDD/land cover. A lower sample size remains descriptive, not a powered treatment test. This 10+10 threshold is a planning heuristic, not power analysis.

**Observed phenology is not railway habitat:** a site 200 m from rail may be a school garden, parking lot, park or sheltered backyard. Without site inspection, mowing logs and gust exposure, there is no defensible inference about railway-specific management or evolutionary adaptation.

## Already completed small feasibility audit

[First 12 geographically frozen pairs per state (48 station requests), run #42389091804](https://github.com/zuizui0223/tetsu/actions/runs/42389091804):
- MD: 12 near + 12 matched far queried successfully, **0** near and **0** far with an Open flowers record in the 2024 query.
- AZ: 12 near + 12 matched far queried successfully, **0** near and **1** far with an Open flowers record.
- No shared eligible species; this is a **source-coverage finding**, not evidence of no biological flowering or rail effects.
- Code/tests subsequently tightened so missing/malformed species IDs cannot silently be counted as zero taxa.

## Full 170-near census (separately labeled new stage)

Run `.github/workflows/npn-active-rail-species.yml` with:
```bash
python scripts/audit_npn_active_rail_sites.py --pairs-per-state 120 --output results/audit/npn_rail_active_sites_2024.json
```

Reports by state:
- number of registered near-site positions sampled vs independently matched far sites;
- number of successfully queried vs failed near station requests;
- number with **any NPN statuses** (actual 2024 activity);
- number with real **Open flowers status records**;
- number containing at least one repeatedly monitored species with both flower present/absent;
- taxa common in eligible near and far sites, with site replication counts; no treatment effect sign used in their selection;
- explicit no-matched-far and geographic filtering counts.

Only aggregate status and taxonomy counts are retained; no NPN site coordinates, station names or individual identifiers are published as CI artifacts.

## What a failed coverage gate means

A zero count for a fully queried, empty site means **no suitable records in this downloaded year/filter**, not that the plants were missing or never bloomed. A failed service request means **unknown**, not zero. Zero common taxa means **cannot estimate a within-species railway contrast from this chosen sampling frame**; it does not refute train-mediated seed transport, environmental filtering or plasticity in field populations.

### Separate research threads retained

- **Japan rail network / seed redistribution**: freight and ORS, plants and ports, independent road and land-cover controls, species assortments and native/alien statuses.
- **Rail mowing and photos**: independently documented cutting events; verified same fixed-quadrat conspecific plant-mask and flower-mask measurement; pre-/post-cutting comparisons with uncut controls.
- **Wind response**: independent measurements of actual train-induced gusts, shielding and plant biomechanics.
- **Phenotypic plasticity vs heritable adaptation**: repeated individual responses vs common-garden or reciprocal transplant/genetic differentiation.

Never claim ecological adaptation from source-proximity or source-availability alone.
