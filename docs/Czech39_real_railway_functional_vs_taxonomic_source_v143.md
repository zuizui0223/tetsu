# v14.3 — Independent Czech 39-station field flora and Pladias four-guild functional–taxonomic identity diagnostic

## Verified original published biological source

Original author study: Kutlvašr et al. (2024), *Railways as a source of alien plants*, Preslia 96:247–266, DOI https://doi.org/10.23855/preslia.2024.247. Authors sampled 39 real railway stations along five source stretches in 2020–22. The paper reports 764 taxa and 11,139 original taxon observations across three station-habitat zones. This is an **original-authors' field dataset secondarily reanalyzed by tetsu**, not new railway field sampling by tetsu.

Original actual publisher XLSX data: [S1 species-by-station and history](https://www.preslia.cz/appendices/P24302KutlvasrSupplTableS1.xlsx) source 79,261 bytes, SHA256 d f 6 b c 0 b 9 4 9 e 7 f 6 a 5 c b 5 7 7 5 5 d c 8 f 9 c 8 6 1 8 7 2 2 3 f 6 6 d d f b e f 9 5 9 c f 0 e 0 3 e f 7 c 4 6 b 0 6 (without spaces); [corrected S2 erratum](https://preslia.cz/appendix/file?id=540), 13,533 bytes, SHA256 00c7cfdf5047854ba99d47dc8ec8c55b1b7df72731ce647c22b286ddef24a5a2. All source file hashes, OOXML checksums, exact 39 station codes and 6,941 contemporary species × station presence marks were verified in successful genuine GitHub Actions.

**Real original publisher-data run: [v14.3 SUCCESS #43669088184](https://github.com/zuizui0223/tetsu/actions/runs/43669088184)**. [Source analysis](../scripts/analyze_Czech_2024_39station_functional_identity_v143.py), [source math and provenance tests](../tests/test_Czech_2024_39station_functional_identity_v143.py), [workflow](../.github/workflows/Czech-real-39station-Pladias-functional-identity-v143.yml). Actual external XLSX downloaded in ephemeral GitHub Actions, not republished.

## Taxonomy/trait coverage audit

**Important integrity difference:** The paper reports 764 taxa across original zone-level study records. The released S1 39 contemporary REC station columns have **761 genuinely marked species**, 6,941 original X marks, and 1,053 total taxon rows including historical-only observations. **Three taxa in the paper's number are not in the published current-station matrix and must NOT be invented.** The original S2 total station richness also disagrees with source S1 in three stations (S1 minus S2: DV +4, SV +1, TI +1). No correction by filling imagined presences. Historical 33 station HIST columns have nonuniform sampling and cannot be treated as comparable current or past absences.

Exact full taxonomic name join with Pladias/Sádlo author-published dispersal archetype and Raunkiaer life form yields **544/761** actual contemporary taxa with BOTH unambiguous traits, **5,081/6,941 =73.20%** real contemporary species×station incidence coverage, minimum 82 trait-known taxa at every station. Source strategy archetypes Epilobium/Phragmites are guild labels, not botanical genus-based name matching. Raunkiaer therophyte is **not** the exact same life-history variable as the Slovenian author's strict annual lifespan code. Four coarse source guilds: therophyte×wind, therophyte×other, non-therophyte×wind, non-therophyte×other. Functional similarity measures proportional overlap of these bins rather than a Gower/continuous hypervolume. Missing traits NEVER assigned to the nonwind group.

## Authentic observed source similarities

39 original stations make 741 dependent pairs, **143 same author railway stretch** and **598 different stretch**; five source stretches belong to ONE country, not five independent experiments.

| Source comparison (source location = entire original station) | Same source stretch | Different source stretch | Same−different |
|---|---:|---:|---:|
| Taxonomic Jaccard of 544 known-trait species | **0.384862** | **0.360920** | **+0.0239415** |
| Four source functional guild proportional overlap | **0.930571** | **0.934462** | **−0.0038906** |
| All author-published S1 REC species (761 taxa) Jaccard | **0.392277** | **0.369262** | **+0.0230143** |

Exploratory 300 correlated Markov-chain checkerboard null draws preserve EACH station's source four-guild species richness and EACH original source species' total frequency across 39 stations. Observed known-species Jaccard difference **+0.0239415**, source null average **+0.0000971**, excess **+0.0238444**, source chain 2.5–97.5-percentile interval **[−0.0046542,+0.0047747]**; plus-one upper reference fraction **0.003322**. This is NOT a randomized train dispersal p-value: original railway stretch membership and habitat were not randomized, MCMC draws are autocorrelated, and the study contains no rail network freight or regional vegetation pool adjustment.

## Critical geographical confounding and limits of interpretation

Czech source S1/S2 lacks independently verified station GPS/coordinates. Entire rail STATION source areas differ dramatically, with three habitats (tracks/embankment, buildings, wider surroundings), not the Slovenian 100m matched BS-vs-BT field method. Thus even an apparent **species identity surplus beyond coarse guild profile** could be a simple product of geography, soil, local flora, habitat/sampling area, urbanization and original railway stretch location. No train-carriage seed samples, cargo OD, historic train activity or same-scale ROAD comparison in S1/S2.

Independently published Slovenian 2026 89 actual 100m railway source transects had raw Jaccard same-section vs other gap **+0.052644**, but author original GPS showed same-section median separation 17.20km vs other-section 75.24km. Predeclared equal-geodesic-and-guild source matching reduced the gap to **+0.009855**, and preserving actual author survey-year×track-type source species frequencies made observed−null excess **−0.008858**, not a robust railway dispersal result. [Real source Slovenian v14.1–v14.2 geographic falsification report](Slovenia_functional_vs_taxonomic_spatial_source_v141_v142.md), [real original source CI](https://github.com/zuizui0223/tetsu/actions/runs/43668179797).

**Conclusion:** The Czech positive raw source guild-adjusted species identity gap is exploratory, not an independent confirmation of actual freight seed dispersal. The original Czech spatial/site grain and Pladias therophyte proxy differ from original Slovenian site/year and lifespan categories. NO same-effect meta-analysis should pool raw +.02394 Czech with original raw +.05264 Slovenia. Credible macro inference requires original station GPS, local regional source plant species pool, geography/environment-matched connected/unconnected rail route pairs, road controls, actual active freight traffic and original site management measurements. Actual native-plant mowing/strong-wind genetic adaptation NOT measured by either study.

Research output in the [Draft PR #1](https://github.com/zuizui0223/tetsu/pull/1); default main not changed.
