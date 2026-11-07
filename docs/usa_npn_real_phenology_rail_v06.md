# v0.6 — Real USA-NPN phenophase × federal rail GIS feasibility (2026-11-07)

**Status:** Source-backed repeated-plant phenophases and rail site proximity, not a confirmed railway flowering effect. Original Japan seed movement and mowing/ROI workflows remain open and **are not negated**.

## Why extend the Japan photograph study to the United States?

The Japan iNaturalist photo-annotation pilot had **185 annotated Flowers vs only 4 annotated No Flowers or Fruits**, and the Tokyo roadside-controlled photographic association was not spatially identifiable (2 complete grid blocks). Photos are observer-selected, and "no flower visible" is not a biological absence.

USA-NPN Nature's Notebook has actual *repeated individual observations*, with status **1=phenophase present, 0=absent, -1=uncertain**, not just attractive flower photos. Its published portal has >47 million phenophase status records across the USA as of November 2026 (dataset totals will update): https://www.usanpn.org/data. Documentation: https://usa-npn.github.io/rnpn/reference/npn_download_status_data.html ; cleaning https://usa-npn.github.io/rnpn/articles/VIII_data_cleaning.html.

## Three real source-backed audits completed

1. [Source/field audit — run #42378936396](https://github.com/zuizui0223/tetsu/actions/runs/42378936396):
   - Maryland 2024 *Liriodendron tulipifera* (species ID 82): **3,229** all-phase records at **11 individual plants / 3 registered sites**.
   - Arizona 2024 *Carnegiea gigantea* (ID 210): **12,957** all-phase records at **51 individual plants / 23 sites**.
   - Combined 16,186 real records with individual_id, site_id, observation_date and phenophase_status. Original API response can be array-format even when NDJSON URL is used; fail-closed parser now accepts correct array/NDJSON object schemas.

2. [Real individual flowering diagnostic — run #42379297172](https://github.com/zuizui0223/tetsu/actions/runs/42379297172):
   
   | 2024 source | Tuliptree, MD | Saguaro, AZ |
   |---|---:|---:|
   | 'Open flowers' raw records | 493 | 2,673 |
   | Valid observed 'yes' | 36 | 174 |
   | Valid observed 'no' | 457 | 2,451 |
   | First observed yes preceded within 21 days by observed no | 10 individual-years | 26 individual-years |
   | Monitored yes → ≥2 distinct no visits → yes, with each transition ≤21 days | 0 candidate individuals | 9 candidate individuals |

   These are **observed-state transitions**. The 9 saguaro sequences represent potential reappearance of open flowers separated by documented non-flower visits, NOT necessarily induced flowering by mowing, train wind, genetic change, or even truly separate continuous flowering episodes. Unknown/conflicting same-day records never become biological zero.

3. [Rail GIS audit — run #42379467494](https://github.com/zuizui0223/tetsu/actions/runs/42379467494):
   - All **3** Maryland monitored sites and all **23** Arizona monitored sites in these *two taxon-specific 2024 series* are ≥1 km from mapped FRA 2026 rail lines. They **cannot** be used for a rail-near vs rail-far test.
   - Rail source: U.S. DOT/FRA National Transportation Atlas Database, updated **2026-10-19**, with PASSNGR, NET, line ownership, track data: https://services.arcgis.com/xOi1kZaI0eWDREZv/ArcGIS/rest/services/NTAD_North_American_Rail_Network_Lines/FeatureServer/0. Later rail vintage is not an accurate historical 2024 exposure without verification.

4. [Registered-site census — run #42379733758](https://github.com/zuizui0223/tetsu/actions/runs/42379733758), **taxon- and flowering-outcome-independent**, using USA-NPN `stations/getAllStations.json?state_code=...`:
   
   | Registered monitoring sites (not active in 2024 necessarily) | MD | AZ |
   |---|---:|---:|
   | All listed | 1,041 | 2,918 |
   | ≤250 m of FRA 2026 railway | **105** | **65** |
   | 250–1,000 m | 272 | 256 |
   | ≥1,000 m | 663 | 2,594 |
   | Excluded geolocation or duplicate ID | 1 | 3 |

   Thus **170 railway-near registered-site candidates**, but neither active 2024 flowering data nor railway verge habitat verified. The initial two species are among the wrong spatial/monitoring cohorts for the rail contrast.

## Hypothesis structure: one research ecosystem, independent evidence tracks

- **Transport hypothesis (Japan-first):** validated freight vs passenger vs road infrastructure + region species assemblages + introductions/ports. NPN phenology does **not** test train-carried seeds.
- **Habitat/phenology hypothesis:** among plants **independently verified to occupy railway verge, roadside and background**, flowering event/interval states differ at matched GDD, species, site-year and observer. USA-NPN can provide repeat phenophases, but station rail distance cannot alone classify the habitat.
- **Mowing plasticity hypothesis:** same marked individual/fixed quadrat before/after **documented actual cutting**; repeat flower display/cover from human-reviewed masks. NPN status times without management dates cannot identify cutting.
- **Train-slipstream/selection hypothesis:** separate measured peak gust/shelter/thermal conditions, with common-garden or genotypes to claim heritable adaptation.

## Next strict eligibility gate (without choosing taxa by effect direction)

1. Enumerate **active repeat-visited NPN station×taxon×year records**, independent of flowering outcome; geography/visit-effort gate first. Registered sites with no observation in a year must be excluded, not treated as nonflowering.
2. Within a **common taxon**, keep both ≥250m separation groups with enough individuals and explicit yes+no records, not purely sites sharing gardens/rail corridor proximity.
3. Prefer **≥10 independent active sites per rail-near vs rail-far group, ≥3 geography/landcover blocks with both groups**, and independent observers, as feasibility thresholds; no significance from insufficient balance.
4. Source independently matching roads, land cover, urbanization, temperature/GDD and historical rail service. Contrast passenger/freight when verified, not FRA network flag guesses.
5. Estimate **observed onset bracket** (prior no → first yes) with interval-censoring, and possible interruption/reflower under repeated surveillance. Do not label first photo as blooming onset or claim causal adaptation.
6. A field/management log is mandatory before attributing phenophase interruptions to mowing.

## Reproducible code
- `scripts/audit_usa_npn_status.py`: bounded API real source record schema.
- `src/tetsu/npn_individuals.py` and `scripts/run_npn_individual_phase.py`: normalized Open flowers yes/no, interval-censored transitions, conservative potential interruption.
- `scripts/audit_npn_rail_proximity.py`: NPN measured-site locations vs FRA network.
- `scripts/audit_npn_registered_sites_rail.py`: independent full registered-site rail proximity screening. No point coordinates or individual IDs published in artifacts.
- Tests and CI workflows under `.github/workflows/usa-npn-*.yml`.

### Scientific literature guiding cautious wording

- Rakotoson et al. (2025), *Plant communities responses and effects in railway environment*, actual Paris RATP plant inventories at 394 locations, 255 observed plant species, demonstrating that track and embankment are different habitats: https://doi.org/10.1111/wre.70025
- *Ecological Processes and Nature-Based Solutions in Urban Railway Corridors* (2026), seasonal phenological observations at three Perth railway settings, confirming railway phenological comparisons are plausible but geographically limited: https://www.mdpi.com/2073-445X/15/5/714
- *Grassland Management in a Changing Climate* (2026), cut height/frequency/timing and plant community response on Dutch road verges: https://doi.org/10.1111/avsc.70085

**Do not report a railway flowering shift or stress adaptation from the above coverage audits.**
