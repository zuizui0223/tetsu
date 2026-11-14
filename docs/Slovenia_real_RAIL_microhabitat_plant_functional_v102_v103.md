# tetsu v10.2–v10.3: Actual Slovenian railway plant microhabitats and life-history filtering proxies

## Source: genuine public railway field-survey data, not MIREN roadside plots

Original authors: **Šabić, Azra; Grašič, Mateja; Rovanšek, Anže (2026)**, *Vegetation patterns, functional traits and invasion dynamics across railway microhabitats in Slovenia*, *Plant Biosystems* 160, article 154, [DOI 10.1007/s44473-026-00155-z](https://doi.org/10.1007/s44473-026-00155-z), CC BY 4.0, [public supplementary file1 XLSX from publisher](https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs44473-026-00155-z/MediaObjects/44473_2026_155_MOESM1_ESM.xlsx).

**Genuine independent source success:**
- [v10.2 author-published 2026 original XLSX obtained and CRC/SHA256 verified #43538799032](https://github.com/zuizui0223/tetsu/actions/runs/43538799032), [unchanged original supplement + machine-readable hash provenance artifact #14090917979](https://github.com/zuizui0223/tetsu/actions/runs/43538799032/artifacts/14090917979), [source code](../scripts/audit_Slovenia_rail_vegetation_supplement_v102.py).
- **[v10.3 site/track-side matched ballast shoulder vs toe functional response SUCCESS #43539536809](https://github.com/zuizui0223/tetsu/actions/runs/43539536809)**, [original-run JSON result #14090649879](https://github.com/zuizui0223/tetsu/actions/runs/43539536809/artifacts/14090649879), [reproducible analysis](../scripts/analyze_Slovenia_rail_microhabitats_v103.py), [8 original-source method tests](../tests/test_Slovenia_rail_microhabitats_v103.py) PASSED. Original supplement is read only; no original workbook cell is edited or rewritten.
- The source supplementary workbook's SHA256 is **26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f**, exactly **401,279 original bytes**. MIREN road data were NOT used to fit the Slovenian model.

### Actual source measurement grain and rigorous quality checks

The original real-RAIL XLSX comprises six sheets: README, 92 survey locations, 482 railway microhabitat field-zone measurements, **6,331 species-by-zone original observations**, 312 species trait rows, and 312 species codebook records. There are 15 real railway line segments; observations made **2022 and 2023** after a 2021 public infrastructure herbicide policy ban. These are not before/after management treatment experiments. Sampling points were intentionally distributed along railway corridors based on visible established vegetation, NOT a random national plot sample.

**Reconciled source discrepancies** (do not hide): 6,331 raw species-by-zone records have **36 repeated EPPO codes within identical zone**. Actual zone richness matches the resulting **6,295 unique species-code×field-zone incidences** exactly for all 482 zone rows, after deduplicating repeated taxon entries. There is an original taxonomic codebook versus traits mismatch: trait/source observations label *Lotus herbaceus* code **LOTHE**, while published codebook gives **DOYPH** for same Latin name; 4 original occurrences use LOTHE. We retain the trait and occurrence source key (LOTHE) for measurement and flag the mismatch, rather than silently substituting codebook names or throwing away biological observations.

One original railway-zone **cover** source value equals **200**; do not call every original cover field a 0–100% percentage. Cover is compared as published **source cover units**; exclude the >100 source observation in a documented sensitivity.

**Microhabitat design:** railway Sleeper Zone SZ (within ties), Ballast Shoulder BS, Ballast Toe BT, Intertrack Area IA (double tracks). BS and BT are adjacent 1–2m transition zones along a 100m rail transect; physical zone width and sampled area are **NOT exactly documented as equal**. For the primary analysis, pair **BS and BT on exactly the SAME survey location AND SAME L/R side**. This retains **173 zone-side pairs from 89 unique railway field locations, across all 15 route sections**. Left and right sides are averaged WITHIN LOCATION before any geographic sign summary. SZ and IA are not treated as equivalent independent BS/BT samples because number of zones and plot areas differ for single/double track; lack of comparable area and route-level replication prevents naive zone-richness causal conclusions.

**Functional variables:** fraction of observed species with source species-level strict Therophyte life form label Th or strict annual lifespan label annual, using denominators of species with source-known relevant trait codes. Mandatory broad-coded alternative includes mixed Th/He and annual/biennial labels. Species with BLANK Nativeness.status remain **UNKNOWN**, not Native. Explicitly 'alien' marked Nativeness.status taxon fraction is an annotation only and not a complete locally validated alien/native composition.

## Actual 89-field-location mean within-location same-side BS − BT differences

| Original RAIL measured microhabitat variable | BS − BT equal-location mean | Same-site locations BS greater | Same-site locations BT greater |
|---|---:|---:|---:|
| Source cover (original arbitrary / source units, not proven %) | **−35.179775** | 7/89 | **81/89** |
| Species richness count per zone (sampled zone areas not exactly equal) | **−2.078652 taxa** | 29/89 | 58/89 |
| **Strict Therophyte ('Th') proportion among trait-classified recorded taxa** | **+0.058871 (+5.89 percentage points)** | **62/89** | 27/89 |
| **Strict annual lifespan fraction** | **+0.054703 (+5.47 percentage points)** | **60/89** | 29/89 |
| Broad coded Therophyte fraction including Th/He etc | +0.049214 | 61/89 | 28/89 |
| Broad-coded annual lifespan fraction including annual/biennial | +0.040012 | 55/89 | 34/89 |
| Explicitly labelled alien taxa fraction of all recorded taxa, unknowns not native | **−0.019014** | 34/89 | **53/89** |

A very clear **RAIL microhabitat structural pattern** is that ballast shoulder exhibits lower reported vegetation cover, while sites with adequate trait metadata have modestly higher relative frequencies of STRICT annual/Therophyte species. Invasive/nonnative-labelled recorded taxa do NOT show positive concentration in the shoulder; this counters any simplistic "harshest track zone has the most aliens" hypothesis. The absolute species richness difference is secondary because true zone sampled areas differ.

### Region-/track-/year-stratified source checks

There are 15 source rail track sections but only **11** with >=3 independent field locations for segment sign summaries; four are sparse and are NOT treated as negative region outcomes.

| v10.3 segment-level sign source among 11 sections with >=3 field sites | BS greater | BT greater |
|---|---:|---:|
| Raw reported cover | **0/11** | **11/11** |
| Strict annual lifespan proportion | **8/11** | 3/11 |
| Strict Therophyte proportion | **8/11** | 3/11 |
| Explicitly alien-coded species fraction | 5/11 | 6/11 |
| Raw recorded zone species richness | 5/11 | 6/11 |

**Survey-year sensitivity**, BS minus BT at field location mean (NOT herbicide-policy treatment effect):
- **2022, 53 locations**: cover **−35.943396** source units, strict annual fraction **+0.064755**, strict Therophyte **+0.071053**, explicit alien fraction −0.027170.
- **2023, 36 locations**: cover **−34.055556**, strict annual **+0.039903**, strict Therophyte **+0.040935**, explicit alien −0.007006.
- **SINGLE-track, 65 locations**: cover **−35.915385**, strict annual **+0.046123**, Therophyte **+0.052095**.
- **DOUBLE-track, 24 locations**: cover **−33.187500**, strict annual **+0.077939**, Therophyte **+0.077221**.
- **Exclude any location with cover source >100, retain 88**: cover **−34.642045**, strict annual +0.054716, Therophyte +0.058660, explicit alien −0.018598.

These subgroup summaries support the direction of relative annual life-form and cover differences, but with only 11 sufficiently sampled track sections, independence and sample selection do not justify global biological generalized p-values.

## What this adds to tetsu and what it cannot establish

**NEW direct RAIL field source**: Unlike prior GBIF state-aggregated records or MIREN ROADside field plots, actual Slovenia railroad ballast microhabitats show source-verified microhabitat structuring of the measured plant composition and cover. These facts are secondary independent reanalysis of 2026 published authors' data, not novel raw data collected by tetsu.

**Complementary independent ROAD plot source**: MIREN v9.9–v10.1 shows high Alien richness near mountain roads, plus recorded turnover vs nestedness among 100m² equal-area Plot1/Plot3 sites across many mountain region codes. That is biologically different from Slovenia's unequal-width ballast zone data. **Do not calculate a causal RAIL-vs-ROAD difference** by subtracting source numbers from these two incompatible sampling designs.

**Unidentified railway mechanisms**: SZ vs BT area/side differences; actual traffic intensity on each segment; wind gust/slipstream and mowing/herbicide dates; ballast exposure/temperature/soil moisture; freight contamination and plant propagule arrival; local seed mass/regrowth; fitness/common garden/genetics. The source article sampled the years after a 2021 herbicide ban, but without source pre-ban microhabitat repeated baseline, effects cannot be attributed to that policy. Higher annual/Therophyte proportion in BS is consistent with habitat filtering but does NOT show selection CAUSED by trains, mowing, slipstream or adaptation.

**Next rigorous macro identification gate**: independently sourced additional countries/railway original field quadrat datasets with the SAME track zone/side/time and trait taxonomy, paired matched road corridors in the SAME ecoregions with measured source survey effort; validate freight line use and direct seed contamination or deposition in seed traps, and independently record mowing, ballast soil moisture/light and gust wind. Then test a prespecified cargo-delivery × habitat-filtering × dispersal/regrowth trait interaction, with rail sections/regions as independent geographical units rather than pairwise pseudoreplicates.
