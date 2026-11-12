# tetsu v8.0–v8.3 — Freight network opportunity and recorded introduced-flora geography are decoupled in state-scale U.S. source data

**Verified 2026-11-12, actual government+GBIF source data, no biological causal claim.**

- [Official FAF6 crop RAIL/TRUCK economic OD network v8.0 SUCCESS #43142647675](https://github.com/zuizui0223/tetsu/actions/runs/43142647675) · [OD source JSON](https://github.com/zuizui0223/tetsu/actions/runs/43142647675/artifacts/13910351407)
- [Census FIPS×GBIF/GRIIS 51-state plant source v8.2 SUCCESS #43144396513](https://github.com/zuizui0223/tetsu/actions/runs/43144396513) · [51-state species source JSON](https://github.com/zuizui0223/tetsu/actions/runs/43144396513/artifacts/13910004390)
- **[State-pair aligned exploratory freight×plant source model v8.3 SUCCESS #43145635279](https://github.com/zuizui0223/tetsu/actions/runs/43145635279)** · [raw results JSON artifact](https://github.com/zuizui0223/tetsu/actions/runs/43145635279/artifacts/13911985381) · [model code](../scripts/analyze_FAF6_GRIIS_US_state_OD_alignment_v83.py) · [6 targeted tests](../tests/test_FAF6_GRIIS_US_state_OD_alignment_v83.py). The complete repository smoke test also PASSED, including regression against independently sourced original data.

## Main question

Does **estimated agricultural freight tonnage moving between U.S. states by rail** correspond to a greater **recorded similarity of GRIIS-U.S.-country-listed plant taxa** between those states, once the same-OD TRUCK tonnage, state geography and indexed GBIF recording-scale covariates are held fixed?

This moves beyond prior 50km physical railway-map connectivity, but **does not observe viable alien seeds, dated introduction or ecological adaptation**. The route is an economic origin→destination state pair, not a particular railway line. The floristic source response is symmetric state-pair similarity: since it contains no dated arrival history, the source directional FAF commodity OD is symmetrized as tons(o→d)+tons(d→o). Directional actual dispersal cannot be inferred.

### Independently verified four-data-source join

1. **2022 BTS official FAF6 state origin–destination source:** 599,529 original records; state-coded 51 geography; 2022 *estimated* tonnage; commodity SCTG02 cereals/grains (including seed) and 03 other agricultural commodities; transport modes RAIL versus TRUCK. Actual v8.0 source ZIP checksum and national tonnage match source v7.9.
2. **2025 U.S. Census TIGER state geography:** 51 state FIPS+state names exactly matched to prior Census and FAF economic state codes; authoritative state representative centroid coordinates, Census REGION/DIVISION, land area. These centroids are geography controls, **not** physical freight route lengths or the location of botanical sites.
3. **2018–25 GBIF US stateProvince indexed source:** all 51 state/DC names queried; 49/51 pass preexisting >=100 all-plant records, >=3 publisher dataset categories, >=5 country GRIIS-listed keys, >=20 GRIIS indexed records and <=2% NOT_DETECTED/unknown source status. Fails are **Alaska (18.50%)** and **Colorado (5.63%)** of indexed statuses non-DETECTED/unknown, NOT biological absence. 2,660 unique GRIIS-U.S.-country-listed speciesKey categories appeared in any of the 51 state source queries. **State-level introduced status NOT verified**; country-listed GRIIS taxon may be native to some states.
4. **Macroscopic source response:** union/intersection of recorded speciesKey categories between each of the 49 eligible source states; predefined primary K20 ratio of expected incidence intersection/union after random drawing of 20 indexed GBIF records within each source state, with raw recorded Jaccard as a mandatory sensitivity. K20 is a source-record sampling model, not true site occupancy/detection adjustment.

A subset of state pairs was chosen only by independent 2022 source RAIL/TRUCK estimated nonzero flow, never based on the floristic response sign. Source pairs reuse states, so hundreds of pairs do **not** constitute hundreds of independent ecological replicates.

### Primary real results

Residualized rank-based (partial Spearman-like) correlation of **log1p(estimated 2022 state-pair agricultural RAIL freight tons)** with recorded floristic similarity, conditional on:
- log1p(estimated same-OD crop TRUCK freight tons);
- Census centroid geodesic separation and Census administrative region/division matching;
- state land-area geometric mean/ratio and GBIF all-Plantae record geometric mean/ratio;
- the two states' overall source RAIL and TRUCK economic crop freight connectivity strength.

**Positive** partial rho is in the railway-transport similarity-hypothesis direction. The regression adjusts measured proxies only; source publisher, observer, climatic, agricultural-area and colonization history remain unmeasured.

| Original 2018–25 GRIIS country-listed species source metric | ALL quality-eligible interstate state pairs (49 states; 1,176 shared-state pairs) | Pairs with positive RAIL crop OD only (141 pairs; 40 states) | Pairs with both RAIL and TRUCK crop OD positive (133 pairs; 39 states) |
|---|---:|---:|---:|
| **Predeclared K20 expected-inclusion similarity** | **−0.035831** | **+0.048024** | +0.062162 |
| Raw observed speciesKey Jaccard | −0.042230 | +0.020715 | +0.009820 |

**Mandatory jackknife sensitivity of all-state primary K20 result**: Leave out every eligible Census state once, then refit from its remaining 48-state pair universe. **49/49 refits negative**, coefficient range **−0.068123 to −0.023088**. This is **within-US node-deletion stability**, not 49 replicated nations nor a valid p-value. Source state pairs share endpoints and all observations depend on GBIF reporting effort and state management histories.

**Scientifically important pattern:** The **economically observed agricultural rail freight network is highly interstate, sparse and destination-hub concentrated (v8.0)**; however this does NOT translate into a general positive national link between reported freight flow and GRIIS-listed plant-category similarity. A small positive signed relationship occurs if we condition on a pair *already having* positive rail freight, but the magnitude is only rho 0.048/0.062 and it disappears/sign-reverses in the full source state pair universe. This suggests a transport opportunity versus recorded flora association **disconnect**, not causal biological filtering.

### What the data do NOT identify

- Genuine local alien plant colonization in a state after a particular dated cargo movement. GBIF plant records cover 2018–25 around a single FAF 2022 source base year, **without treated-before/treated-after exposure**.
- Species-specific viable diaspores actually released from trains, containers, grain shipments, vehicles, rail freight yards or port connections. **FAF tonnage is not viable-seed counts.**
- Specific *physical* railway connections carrying source freight: the FAF state OD estimate is not an Overture railway shortest-path traffic allocation.
- Observer-neutral, survey-verified state plant absence or stage-specific introductions. GBIF stateProvince is an interpreted source text field, not spatially validated Census state polygon records.
- Railway-edge **mowing, wind gust stress, herbicide, ballast heat**, regenerative/reproductive trait selection, phenotypic plasticity or genetic/evolutionary adaptation.
- Pseudoreplicated significant p-values: 1,176 state pairs are functions of 49 shared state plant lists. Country is one broad continental unit, and cross-sectional source effects cannot establish universality.

### Next decisive macro research gate

Stop adding more posthoc binary rail-mapping metrics. For **delivery**, acquire independently dated actual freight line/service flow or validated origin cargo contamination/seed assays and direction-specific *post-arrival* plant records. A commodity control such as resource SCTG10–12 can be assessed on the same state pairs to ask whether any association is **agriculture-specific**, but it remains a negative-control *candidate*, not proven sterile freight. For **filtering**, acquire repeated rail and road verge plant surveys with actual documented mowing/gust/thermal regimes and species-level regrowth/seed aerodynamics, then test a preregistered cargo-supply × disturbance × trait interaction. For **evolution**, require within-species repeated phenotypes and common-garden/genetic evidence.

The present supported result is a descriptive transport network structure and a source-controlled ecological **non-generalization**; do not conflate it with trains actually planting alien species.
