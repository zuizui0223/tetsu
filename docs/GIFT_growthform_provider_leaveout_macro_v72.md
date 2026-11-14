# tetsu v7.2 — The woody-flora railway signal does not survive GBIF source-provider exclusion

**2026-11-11 confirmed with actual independently preexisting GBIF/GRIIS/GIFT/Overture source artifacts, not simulated ecological data.**

- Original six-country growth-form macro source and 50km RAIL/ROAD routes: [v7.1 run #42980383943](https://github.com/zuizui0223/tetsu/actions/runs/42980383943).
- Same-country/cell/trait FULL-vs-WITHOUT dominant GBIF dataset source reanalysis: **[v7.2 live SUCCESS #43009955904](https://github.com/zuizui0223/tetsu/actions/runs/43009955904)**, [source artifact](https://github.com/zuizui0223/tetsu/actions/runs/43009955904/artifacts/13850570616).
- Code: [analyze_growthform_source_leaveout_v72.py](../scripts/analyze_growthform_source_leaveout_v72.py), including fixed source eligibility / both-growth-group / both-provider same physical-pair checks.
- This is a **descriptive 6-country source-sensitivity macro study**, not proof of seed dispersal or plant adaptation.

## Earlier positive source-only macro finding: growth form (v7.1)

From 2,356 distinct 2023-Backbone GRIIS-listed GBIF keys detected in original six-country macro-cell source, 2,089 received an unambiguous species-level GIFT 1.2.1 herbaceous/woody classification (**88.67% union mapping**). Across six selected country blocks, five had sufficient same-physical-pair source data for both groups under 2018–2025 nearest-in-cell route anchoring.

- Herbaceous: railway-hypothesis compatible *negative* partial rho in 3/5 independent country blocks, country-median **−0.052343**.
- Woody: negative in **5/5**, country-median **−0.144935**.
- Herb−woody descriptive partial-rho contrast positive in 4/5 (country median **+0.053926**).
- The more restrictive **nearest-larger-component** physical RAIL anchor already changed the woody result to 4/6 negative (median **−0.053687**) and herbs to 2/6 negative (median **+0.086353**); the older 2010–17 period had only 3 qualifying woody blocks, 2/3 negative. No independent biological mechanism identified.

## New stricter real recorded-source leaveout — primary 2018–2025, nearest-in-cell, K20

For each original eligible 50km source cell, the *largest ALL-Plantae GBIF dataset* was selected without looking at rail/plant correlation. The country-GRIIS speciesKey counts were partitioned into FULL and EVERY OTHER DATASET from real indexed GBIF source records (original v5.6). The source v7.2 output holds **exactly the same 50km cell-pair cohort**, both groups, physical RAIL/ROAD route endpoints, and geographic/ROAD/all-plant recording-scale covariates when replacing FULL with WITHOUT dataset. The larger original v7.1 all-record results are not on the same restricted cohort.

| Trait group and independent country-block direction | FULL source on common cohort | After single largest dataset excluded |
|---|---:|---:|
| Herbaceous, railway-hypothesis compatible *negative* signed country rhos | **3/4** | **2/4** |
| Woody, negative country rhos | **4/4** | **2/4** |
| Herbaceous median country partial rho | −0.070427 | +0.038920 |
| Woody median country partial rho | −0.091771 | **−0.024878** |
| Herb−woody median same-pair rho contrast | +0.050535 | +0.066678 |

**Individual woody coefficients (FULL → WITHOUT, strictly identical pair cohorts):**
- **US Midwest:** −0.003709 → +0.154987 (sign reversed)
- **France Burgundy:** −0.086555 → +0.089708 (sign reversed)
- **UK Yorkshire:** −0.096987 → −0.445059 (more negative)
- **Germany Saxony:** −0.144935 → −0.139463 (largely unchanged)

**Coverage gate:** after requiring at least 8 independent 50km source cells and at least 20 within-region physical RAIL/ROAD paired routes (not independent biological replicates), only these **4/6** macro blocks remain estimable for both forms and both GBIF source groups. Canada Ontario had 16 same-cohort connected physical pairs, despite 8 distinct paired cells; AU Southeast only 7 distinct paired cells on valid physical pair routes. Nonestimability does not mean a null ecological effect.

**Anchor sensitivity:** For nearest-of-largest-component source rail network anchoring in 2018–25 and the identical-provider/group K20 cohort, FULL woody source negative signs 3/5, WITHOUT 2/5. Median woody partial rho shifts **−0.003709 → +0.089708**. Historical 2010–17 K20 primary nearest-in-cell has only 2 eligible country blocks, so no older-period macro replication.

## Scientific conclusion

1. **Supported**: growth-form-stratified indexed floral source similarity can exhibit an apparently multi-country woody-RAIL pattern, but this signed pattern is not stable to independent source-provider composition or physical network anchoring.
2. **Not supported**: a generalized causal woody-plant railway transport advantage, selective establishment caused by grass-cutting or train wind, morphological plasticity, or local adaptation. The loss of 5/5 source-level pattern under provider leaveout blocks any such ecological claim.
3. **Next testable theory**: connect a separately measured active/freight railway *delivery* process to species-specific propagule traits and independent documented rail-verge mowing/gust/heat *filtering*, with road-only and nonoperated-line controls, site×species repeat observations, and country×ecoregion replication; use recorded genetic/field phenotype contrasts only for evolutionary adaptation. Do not manufacture wind or mowing stress variables from rail density.
4. **Prior data limitation**: GBIF indexed facets may contain NOT_DETECTED categories that cannot be assigned to individual species without original event records. Source dataset UUIDs are not independent observer or field-survey identities; removing the largest provider is a robustness diagnostic, not perfect detection adjustment. Original study relies on physical 2026 mapped rail routes for plant records from 2010–2025.

No pairwise p-values were computed; 50km pairs share cells. No train/freight historical operating status or actual seed flux is measured. Do not treat GRIIS nonmembership as native status or unknown traits as negative records.
