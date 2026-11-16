# v0.8 — NPN national Site Phenometrics joined to rail geography (2023–2025)

**Validated source run:** [GitHub Actions #42390409359](https://github.com/zuizui0223/tetsu/actions/runs/42390409359). Source artifact `usa-npn-2023-2025-site-openflowers-rail-eligibility`.

## Critical API correction and provenance

The first [v0.8 experimental run #42390038083](https://github.com/zuizui0223/tetsu/actions/runs/42390038083) sent `state[1]` but the `getSiteLevelData.ndjson` source **ignored the state filter**. For each year, Maryland and Arizona responses were **byte-identical**, including SHA256 and row count. Those totals are **US nationwide**, not state totals. They must **never** be presented as Maryland or Arizona counts.

The corrected code [`scripts/audit_npn_site_phenometrics_rail.py`](../scripts/audit_npn_site_phenometrics_rail.py) intentionally downloads **one nationwide source per year**, filters to `phenophase_id=501` (*Open flowers*), joins each record's site ID against independently retrieved **MD or AZ registered station IDs** and FRA rail proximity, and records national totals separately from state-matched site×taxon entries.

- Nationwide derived Open flowers site-phenometric rows: **3,770 (2023)**, **3,670 (2024)**, **3,561 (2025)**.
- State-matched *unique site×taxon* source rows: Maryland **48, 39, 30** respectively; Arizona **310, 226, 184** respectively.
- Site-level phenometrics represent **observed positive flower-status events**, possibly averaged across several monitored individuals, **NOT direct 0/1 observed nonflowering denominators**, true first blooming or documented cutting effects. Reference: https://usa-npn.github.io/rnpn/reference/npn_download_site_phenometrics.html.

## Source-backed SAME-TAXON coverage (not an ecological effect)

These are recorded *positive site phenometrics* joined to independently assigned 0–250m vs ≥1000m distance from 2026 FRA railway line. Listed taxa are **NPN numeric IDs** until validated against NPN's taxonomy API.

| State | Year | Candidate NPN taxon ID | Near registered sites with positive site phenometric | Far registered sites with positive site phenometric |
|---|---:|---:|---:|---:|
| Maryland | 2023 | 7 | **2** | **8** |
| Maryland | 2024 | No qualifying shared taxon | 0 | — |
| Maryland | 2025 | 7 | **2** | **5** |
| Arizona | 2023 | 317 | **2** | **21** |
| Arizona | 2024 | 317 | **2** | **14** |
| Arizona | 2025 | 317 | **2** | **11** |

Other less replicated candidate IDs exist in AZ, but *taxon ID 317* consistently has the most complete near/far **data-availability** overlap among the predefined AZ states/years. This ordering was based on site replication, not flowering date or change direction.

## What still prevents a railroad phenology inference

1. **Only two near sites** per leading taxon in these state-year slices. That is not a powerful spatial contrast regardless of how many individual-day observations are made at the same sites.
2. No verification that the individual plants were growing along the rail right-of-way, affected by rail maintenance/weed cutting, train gusts, ballast heat, or active 2023–2025 rail traffic.
3. The geographic spans of the far comparators may be large. **Before computing any difference in blooming date**, match actual taxon-positive far sites to taxon-positive near sites using a maximum 75km spatial separation, then verify true 0/1 individual records at those sites. Site phenometrics are only a screening stage.
4. Any cut/flower ROI claims require independent site visits and repeatedly imaged conspecific plots under documented cutting events.
5. Train-carried seed dispersal is a separate network/community hypothesis and is **not addressed** by NPN flowering observations.

## Reusable modules and validation

- `scripts/audit_npn_active_rail_sites.py` with full 170 near registered-site 2024 actual-status census. [Results](npn_active_rail_2024_results_v07.md).
- `scripts/audit_npn_site_phenometrics_rail.py` for nationally sourced annual positive site-phenometric coverage, explicitly corrected for ignored source state filter.
- `scripts/audit_npn_candidate_flowering_status.py` for **independently selected taxon/site/75km-comparable-site pairs**, rechecking raw same-individual 0/1 records. **No railway phenology effect size** is permitted until there are enough independent site pairs and credible climate/road/habitat controls.
- CI retains source checksums, failure records, and only aggregated site/taxon counts, without precise station coordinates and individual identifiers.

This dataset represents a reproducible **partial positive data feasibility** finding (observed same species on both rail distance strata in some years), not a positive railway biological hypothesis test.
