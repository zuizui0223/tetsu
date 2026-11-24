# Railway typical taxa and climate-offset hypothesis — research handoff (2026-11-22)

## Status and provenance
Exploratory synthesis of earlier source audits and the 2026-11-22 conversation. **This commit does not rerun the original species-by-study analysis**: prior local scratch files and derived CSVs are not available in this execution environment. Numerical species counts below are **previously reported, pending independent reproduction from archived original inputs**. Do not label this commit a new source-verified calculation.

## Reported typical-taxon sensitivity
| Rule | Prior reported count |
| --- | ---: |
| Positive railway-minus-reference frequency in at least four studies, original name matching | 25 |
| Same, checked synonym reconciliation | 26 |
| Positive rank difference in at least four studies | 10 |
| In at least four studies simultaneously: rail occupancy >= 0.20, frequency excess >= 0.10, rank excess >= 0.10 | 4 |

Previously reported stringent candidates: *Erigeron canadensis* (6/6 positive frequency; 4 stringent), *Erodium cicutarium* (5/5; 4), *Lactuca serriola* (6/6; 5), *Lotus corniculatus* (4/5; 4). These denominators count **eligible** studies; missing/unmatched reference observations must not be coded as negative evidence. A previous sensitivity check reported 25 -> 16 candidates when changing the Czech reference from urban to roadside vegetation. These counts are exploratory and not corrected for selection/multiple testing. The six railway studies are not six matched nonrail experiments.

## Priority hypothesis A — climatic offset, not just warm-affinity taxa
Primary response for comparable rail and local nonrail vegetation: Delta T = community-weighted thermal affinity(rail) - community-weighted thermal affinity(nonrail), on a common trait scale. Primary test: whether Delta T increases as regional long-term growing-season temperature decreases (study/region hierarchy, spatial dependence, survey area, plot year, habitat openness, urbanization and substrate controlled where supported).

Competing predictions:
- Positive, climate-invariant Delta T: rail habitat filtering, **not** evidence of stronger climate decoupling in cold sites.
- Larger Delta T in colder sites: supports a *climate-dependent* rail contrast; still does not alone prove seed movement, colonization beyond physiological limits, or long-term population persistence.
- No contrast after matched controls: no detectable railway-specific thermal offset.

Data gate: Slovak 2016 19 units cover 484–673 m, with geographic distance and elevation difference highly collinear (r approximately 0.992); cannot separate both effects convincingly within that line. Slovenian source has numeric temperature indicator for 142/312 taxa; 135 have nonnumeric x, 35 blank; ~49% of shoulder and toe incidence is covered. Do not impute x as a thermal value. Report coverage, complete-case and missingness sensitivity. Slovenian shoulder vs toe is an *within-railway* contrast, not a rail-vs-nonrail counterfactual. Temperature indicators are ordinal affinities, not degrees Celsius. Prior papers already discuss warm-affinity plants on railway substrates; novelty requires a **matched cross-climate interaction**, not a simple mean difference.

## Alternative hypothesis B — habitat replacement vs incidental occurrence
Czech study reports 85 red-list taxa across sampled railway-related zones, not 85 confirmed ballast-dependent endangered species. Partition threat categories, native status, and actual plot zone. Test whether taxa associated with declining gravel-bar/dry-grassland/sandy habitat persist or reproduce in railway habitat beyond matched local source-pool expectations. Occurrence alone cannot establish a demographic refuge or habitat substitution. Preexisting railway-grassland conservation literature makes a simple richness claim non-novel.

## Next executable source-backed work
1. Retrieve the original archived source matrices and previously generated candidate CSV (if accessible), verify SHA/checksums and taxonomy; reproduce per-taxon eligible-study denominator and all four counts.
2. Freeze a single nonrail reference design before inspecting rail identities. Keep Czech urban vs roadside sensitivity explicit; separate Kashmir paired within-site controls from cross-country MIREN roads.
3. For A, construct a plot-by-species incidence and temperature-affinity table with observed trait coverage and an explicit missingness sensitivity; evaluate whether any *genuinely local* nonrail comparison exists at the same grain.
4. Run matched rail–nonrail climate interaction only when >=2 independent climate regions and adequate plot-level controls exist. Use leave-one-study-out evaluation and spatial/block resampling; do not treat pairs or repeated survey years as independent studies.
5. For B, audit plot-zone-specific red-list records, then join original habitat affinities and historical habitat-loss layers; distinguish richness, persistence and reproduction.

## Reporting restrictions
Do not claim train-mediated dispersal, climate-driven range extension, a railway-specific microclimate mechanism, or demographic refuge from observational taxon-frequency comparisons alone. This is a design and evidence-gate document; **no new ecological model was fitted in this commit**.
