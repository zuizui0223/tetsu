# v0.8 — Real individual 0/1 flowering eligibility after species/data and geography screening

**Executed 2026-11-08 JST.** [Successful source run #42390743089](https://github.com/zuizui0223/tetsu/actions/runs/42390743089), artifact `npn-2023-2025-common-taxon-open-flowers-individual-status-check`.

**What was actually validated:**
1. Annual national USA-NPN derived Site Phenometrics for **Open flowers** (phase ID501) were downloaded once for each of 2023, 2024 and 2025; registered Maryland and Arizona station IDs were joined locally (the upstream site-phenometric state query filter does NOT work as expected).
2. Within each state×year, select a common taxon **solely on source availability**: at least 2 rail-near (≤250m) sites and ≥3 rail-far (≥1000m) sites with a positive Open flowers site-level phenometric. Rank by minimum number of sites per group, then total sites, then taxon ID. The sign or timing of blooming was **not used** to select the taxon.
3. Choose the nearest **unique** far stations to each near station within **75km**, up to 2 far sites per near site. Railway exposure is relative to a *registered station centre* and FRA 2026 line geometry, NOT a known rail-side plant location.
4. Re-query official NPN **individual Open flowers status** (`station_id`, `species_id`, `phenophase_id=501`, year), and check each returned site/species/year/phenophase filter for correctness. Require ≥4 distinct observed dates from the same individual and at least one observed YES and NO; missing/unknown statuses never become NO. All sites queried without errors.
5. Only aggregate counts are exported, not site/individual IDs or exact coordinates.

## Verified species identities

Source `https://services.usanpn.org/npn_portal/species/getSpeciesById.json` identified:
- Maryland selected taxon ID **7** = *Cercis canadensis* (eastern redbud; woody tree)
- Arizona selected taxon ID **317** = *Prosopis velutina* (velvet mesquite; woody tree)

**These are NOT herbaceous railway mowing adaptation candidate plants.** They are valuable for validating actual individual phenology status retrieval, but do not prove railway vegetation management or plant wind tolerance.

## Actual eligibility results

| Taxon × year | Sites with positive site-phenometric near rail | Positive far sites statewide | Far stations geographically matched ≤75km | Near sites passing 4-date yes/no | Matched far sites passing 4-date yes/no | Eligible near individual-years | Eligible matched far individual-years |
|---|---:|---:|---:|---:|---:|---:|---:|
| *Cercis canadensis*, MD 2023 | 2 | 8 | 4 | **2** | **4** | 2 | 13 |
| *Cercis canadensis*, MD 2025 | 2 | 5 | 2 | **2** | **2** | 2 | 2 |
| *Prosopis velutina*, AZ 2023 | 2 | 21 | 2 | **2** | **2** | 7 | 5 |
| *Prosopis velutina*, AZ 2024 | 2 | 14 | 2 | **2** | **1** | 7 | 3 |
| *Prosopis velutina*, AZ 2025 | 2 | 11 | 2 | **2** | **2** | 7 | 4 |

**Independent site sample size remains tiny** (just two rail-near sites for each species/year). Repeated individual-year counts do NOT raise independent railway treatment replicates. For MD 2025, closest match was about 73km away, which leaves serious climate and geography confounding. All three Arizona years had a closest pair ~1.4–2.2km apart, but the different plants may still be in unlike habitats and their railway access is unverified.

## Additional descriptive raw record structure

Actual individual status records allow the source-calibrated algorithm to identify:
- first observed YES with an observed prior NO within 21 days (an *interval bracket*, not exact biological first flowering);
- potential monitoring sequences of YES → at least 2 documented NO visits → YES. These may be phenological reflowering, but can also reflect observation disagreements, distinct flowers on a plant, or changes in observation schedule. No cut events or rail gust data are supplied.

| Taxon-year | Near first-yes brackets | Far first-yes brackets | Near potential interruptions/reappearances | Far potential interruptions/reappearances |
|---|---:|---:|---:|---:|
| *C. canadensis* MD 2023 | 2 | 6 | 0 | 0 |
| *C. canadensis* MD 2025 | 2 | 2 | 0 | 0 |
| *P. velutina* AZ 2023 | 6 | 4 | 1 | 1 |
| *P. velutina* AZ 2024 | 7 | 3 | 0 | 3 |
| *P. velutina* AZ 2025 | 7 | 3 | 0 | 3 |

These counts are **not a statistical hypothesis test of rail effects** and were not used to select candidates. Apparent near/far differences in proportions cannot be interpreted without repeated independent plots, within-year visit schedules, matching climate and land cover, and verified railway vegetation exposure.

## Hard ecological inference limits

- **Tree phenology, not grass cutting:** Neither candidate is a valid surrogate for annual roadside weeds under repeated mowing.
- **Spatial site unit:** NPN registered station coordinate may represent a garden, park, schoolyard, housing area or monitoring plot. A 250m nearest-rail distance is not evidence of the species being on railway land.
- **Historical exposure:** Source FRA network is dated 2026; no verification that rail operations, train frequency, or maintenance were active at these sites in 2023–2025.
- **Confounding:** 75km nearest-station matching by itself does not balance land cover, elevation, roads, heat island, precipitation or growing-degree days; the selection also conditions on having at least one observed flowering event, so estimating flowering *probability* from this matched subset would be biased.
- **Replication:** Two rail-near registered sites per taxon/year cannot distinguish site idiosyncrasy from an underlying railway-specific response. Even if effects look directional, reserve them as exploratory unpublished diagnostics.
- **Evidence sources:** USA-NPN status documentation https://usa-npn.github.io/rnpn/reference/npn_download_status_data.html, NPN site phenometrics https://usa-npn.github.io/rnpn/reference/npn_download_site_phenometrics.html and 2026 FRA railway centreline GIS `https://services.arcgis.com/xOi1kZaI0eWDREZv/ArcGIS/rest/services/NTAD_North_American_Rail_Network_Lines/FeatureServer/0`.

## Next justified paths

1. **Rail-specific mowing / plasticity:** prioritize independently verified railway, road and uncut plots with repeated fixed-camera *herbaceous* plant and flower masks plus actual dated cutting logs; use the existing `patch_pass` ROI and `summarize_mowing_photo_roi.py` pipeline.
2. **General flowering phenology:** replicate independently observed NPN active site×species comparisons in additional states and with more rail-near sites, then fit climate/GDD and land-cover corrections after verifying actual plot exposure. Avoid choosing sites because flowering appears earlier/later.
3. **Transport:** maintain the Japan freight/port/rail network and alien plant assemblage test separately from flowering. A tree's observed open flower is not evidence of train seed dispersal.

**Claim enabled:** Repeat monitored yes/no observations and some geography-matched same-species source cohorts **exist**. **No claim enabled:** Train-induced seed transport, grass cutting recovery, wind stress adaptation, phenology shift or heritable railway local adaptation.
