# v365 — 実線路最短経路は地理距離以上に植物相の類似度を予測するか

**2026-11-24に、v365の予測結果を見る前に固定。** 継承したv364ソースは、本物のSlovenia 2026原著89地点・15原著線路区間、Overture 2026-11-18.0実connector源の測地線長補正済み。元v364 [Actions #45358570723 SUCCESS](https://github.com/zuizui0223/tetsu/actions/runs/45358570723) は、旧UTMの距離が中央値39.4m最大100.9m変化しても、v359の27地点/7元source区間マッチとJaccard+0.0270/年×layout頻度null残差+0.00848に変化が無いことを確認した。**このv365は同一路線ラベルの肯定結果を再利用するより、量的なsource物理network shortest distanceが予測力を増すかを検証**する。

## 固定データとsource consistency
- True original Šabić et al. (2026), DOI 10.1007/s44473-026-00155-z, original XLSX 401279 bytes SHA256 26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f. Real 89 botanical physical field locations within 15 source author sections; pool BS shoulder+BT toe union (primary); BS-only and BT-only as separately labelled sensitivities.
- Overture physical GERS source `2026-11-18.0`, local source field original source 8886 rail-subtype, 6612 conventional active-candidate segments with bridge/tunnel maintained, 12346 connector nodes/14181 source edges. Physical source-site assignment 86 unique, 3 ambiguous; no ambiguous imputation.
- Real original GERS `connectors[].at` weighted by source WGS84 ellipsoid geodesic length; closest source original site GPS projected to nearest segment with original WGS84 source fraction reconstructed, not lonlat planar nor UTM line-length normalized. All 4334 original ordered pairs with network path ≤300km reaudited, no newly dropped/gained source paths vs v359. Remaining original pairs without physical source paths are MISSING, not 0km and not inferred disconnected plant populations.

## Estimand (before viewing v365)
For each **distinct pair of original field sites**, if both real source physical track assignments exist and GERS physical network path ≤300km, response **taxonomic Jaccard similarity** of two real BS+BT field plant lists. Predictors:
- geo: `ln(1+d_geodesic_WGS84_km)` and its square, same source survey year indicator, same original SINGLE/DOUBLE track layout indicator;
- effort/richness sensitivity: mean log(1+source observed site richness) and absolute difference log(1+observed richness);
- **new independent physical network distance**: `ln(1+d_original_GERS_shortest_path_km)`;
- original author line membership: same source section dummy, a THIRD SEPARATE comparator (not a named OSM route or actual train OD).

Models fixed:
M0 = geo log+square + same source survey year/layout + source richness mean and abs difference;
M1 = M0 + true source network log distance;
M2 = M0 + original author same-section indicator;
M3 = M1 + original author same-section indicator.
Sensitivity S1 without richness covariates for all four models, same source pairs.

Technical fitting: weighted ridge least squares on Jaccard (predictions clipped to [0,1]), source-pair weights averaged inverse within-training incidence of the two source sites, then normalized to sum1. Intercept NOT penalized, each predictor train-weight standardized, penalty `lambda=0.05` on each standardized slope; **same hyperparameters for all 4 models and no tuning after results**. This is a descriptive conditional prediction of site-pair Jaccard, NOT estimated individual fitness, vegetation abundance, niche, seed transport, causal regression, generalized additive spatial field or exhaustive variable selection.

### Entire-source-line heldout transfer (not arbitrary pair random CV)
For each of 15 original author rail sections, train on source reachable pairs with **neither** endpoint in the heldout source line; test on reachable source pairs with **at least one** endpoint in that originally heldout source line. Test pairs with same source sites only in one heldout source section; no field site observed taxon sets in model training if that source site is in heldout line. Minimum train pair records ≥200 and heldout test pair records ≥10; otherwise omit that original section and mark unsupported, not zero model gain. Pairwise comparisons are **not independent observations**; report equal-heldout-section average test mean-squared-error and MAE, number of heldout sections with positive MSE gain, LOO across sections and original source fold support. The same test pairs across models, and the source-year/layout info is known even at source-heldout sites. A pair between heldout new field site and previously observed site is not a completely novel two-site pair, so prediction validity is scoped to line-heldout and NOT global-country transfer.

### Predeclared contrasts
Primary incremental predictive value for physical network path over geography: `delta_MSE=MSE(M0)-MSE(M1)`; positive means physical network length helps; negative is a proper failed network hypothesis. Secondary after source same-section label: `delta_MSE=MSE(M2)-MSE(M3)`, positive means network path provides independent information when author line identity already known. Report M0–M3 test losses, source equal-section line improvements positive/negative, and no p-value based on dependent ~2000 source pairs. Technical "n" for inference = at most 15 source original line sections nested in ONE country, not 4334 ordered track pairs.

Before evaluating ecology, audit actual heldout-support and source geographic/network distance predictor ranges; if original identified source sections <5 or representative train/test support low, return NOT_IDENTIFIABLE rather than an effect. Named OSM `route=railway/tracks` relation is not equivalent to the author's `section`; separate v362 original named-line analysis has limited support (unique-route panel only 4 anchors/2 route relations) and strict source year/layout frequency null nonexceptional (any-overlap residual +.01461, technical rank .181).

### Interpretation
Even if true network path improves heldout ecological predictability, no causal demonstration of train seed carriage, railway environment adaptation or within-species evolution. Network shortest path correlates with shared regional flora, mountain valleys, urbanization, substrate, same railway operator and management. A valid macroecological claim would require independent plot×species rail-biotope surveys across other countries and matched nonrail habitats, real line/traffic OD and climate/land-cover controls. Czech source 39 stations have original 5 publisher stretches and species lists but no authoritative contemporary original GPS; do not invent physical graph coordinates. Unknown source train passenger/freight services remain unknown.

Code `scripts/evaluate_Slovenia_global_network_distance_heldout_predictive_v365.py`, tests, true original-source Actions workflow. Only aggregate study results committed, not raw species, source plot GPS or original rail geometry. Draft stacked on PR #58, main unchanged.
