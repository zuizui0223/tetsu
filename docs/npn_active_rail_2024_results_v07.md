# v0.7 verified 2024 USA-NPN rail-near site census — results and identification gate

**Full reproducible source run:** [GitHub Actions #42389422185](https://github.com/zuizui0223/tetsu/actions/runs/42389422185), artifact `usa-npn-2024-rail-near-matched-active-site-species`. Run successful; no API errors. Date of execution 2026-11-07 UTC (2026-11-08 JST). Script `scripts/audit_npn_active_rail_sites.py`.

**Scientific status:** actual source coverage and taxon-eligibility counts only. These are **not railway effects on flowering, mowing response, wind adaptation, train-mediated dispersal, or genetic changes.**

## Previously reported 170 registered near-rail sites: all checked for real 2024 records

| Exact eligibility gate | Maryland | Arizona | Combined |
|---|---:|---:|---:|
| Registered NPN sites within 250m of 2026 FRA rail | **105** | **65** | **170** |
| Near registered sites queried for real 2024 records | 105 | 65 | **170** |
| API read successes / failures | 105 / 0 | 65 / 0 | 170 / 0 |
| Near sites with **any actual 2024 NPN status record** | **3** | **4** | **7** |
| Near sites with ≥1 real 2024 "Open flowers" status record | **3** | **3** | **6** |
| Near sites with ≥1 same-taxon individual having 4 distinct observed flower-status dates including YES and NO | **0** | **2** | **2** |
| Raw NPN status rows in near-site census | 139 | 8,270 | **8,409** |
| Raw Open flowers rows in near-site census | 14 | 514 | **528** |

The status rows are **observations, not independent plants or sites**. A small number of AZ monitoring sites generate a disproportionate number of raw observations, so a large total record count must not be equated with many independent railway site replicates.

## Geographically matched far-site comparison

Using uniquely chosen nearest registered site >=1km from mapped rail, same state, <=75km separation:
- Maryland **101** near/far registered pairs; 4 near candidates had no far registration within 75km, but **all 105 near candidates remained in the independent near-site census**. Among 101 far sites, **8** had Open flowers observations, **7** had a species meeting repeated YES/NO eligibility, across **9 eligible species IDs**. However, **none** of the rail-near MD sites had an eligible taxon.
- Arizona **65** geographically matched pairs; **1** far site had any Open flowers observation and **0** had an eligible repeated YES/NO species. **2** near sites contained **3 eligible species IDs** but no eligible far counterpart.
- **0 shared eligible same-species contrasts** in either state for this single year. No pre-/post-rail or rail treatment analysis was fitted.

## Interpretation

The initial v0.6 registered-site count appeared promising but **registered sites were not valid denominators for observed repeated phenology**: 170 registered near points yielded only 7 actively observed sites in 2024. This is a **data-source/observability constraint**, not negative evidence for railway effects. The initial 12-pair sample (run #42389091804) saw 0 near active Open flowers; full near census avoided falsely generalizing from 24 selected near registrations.

Even if the same species had been observed in both bins, registered-site proximity to a railway centreline does not establish a railway-adjacent **plant**, rail maintenance exposure, train wind exposure, or **active train service in 2024** (FRA line GIS is from 2026). A valid railway trait study needs independently confirmed site conditions, matched road/urban land cover, and actual mowing/wind data.

**Do not select favorable flowering shift direction after viewing these availability results.**

## Next source-backend upgrade: years and actual activity instead of registrations

Use independent USA-NPN **Site Phenometrics** `getSiteLevelData.ndjson`, filtered to Open flowers (`phenophase_id=501`) over fixed **2023, 2024, 2025**, *before* looking at flowering dates or changes. This quickly identifies **actual site × species with a positive site-phenometric record** across the entire state. Such records are not a valid absent-flower denominator and **cannot directly estimate railway-associated flowering shifts**. Only if a species is present in both rail-near and far active sites do we retrieve raw individual 0/1 status observations and check repeated visits and measured environmental contexts.

The detailed fixed eligibility definitions are in [`npn_active_rail_eligibility_v07.md`](npn_active_rail_eligibility_v07.md). Never reinterpret NPN site phenometrics or real status coverage as proof or refutation of the independently maintained Japan rail/freight seed-distribution hypothesis.
