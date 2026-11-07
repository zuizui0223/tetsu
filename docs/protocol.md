# v0.1 Study design — trains as a possible seed dispersal network
2026-11-07. **Pre-results design**, not externally preregistered.

## Three competing processes
- T: seeds physically transported by rolling stock, freight spillage, passengers or maintenance.
- H: tracksides are disturbed habitats favoring ruderal and alien plants, without transported seed.
- U: rail co-location with roads, ports, trade, cities, horticulture, climate and GBIF observer effort.

A nearest-track correlation cannot tell T apart from H or U.

## Analyses
**A — local habitat:** compare rail-side, road-side, and matched surrounding habitat for standardized alien presence/richness and native controls. Match urbanization, land use, climate and all-plant recording effort. This identifies habitat associations only.

**B — network signature (primary):** for pairs of independently sampled 5–10-km cells or local floras, calculate shared alien and native plant species (Jaccard and sensitivity metrics). Use actual network connectivity and rail travel distance, not geometric line proximity alone. Fit:

- M0: geographic distance + road accessibility + climate difference + urban / land-use / port context + recorder effort.
- M1: M0 + railway-side proximity.
- M2: M1 + verified active rail connection and travel impedance.
- M3: M2 + freight flows/yard function, port-to-inland access and introduction-pathway interactions, **only if independently measured**.

Report whether M2/M3 improve **spatial-blocked out-of-sample prediction**, and whether the alien-versus-native contrast is robust. Pairs share endpoints; ordinary independent-pair p-values are invalid. Fit species-level/region random effects or use endpoint/spatial-block uncertainty and geography/degree-constrained network permutation. Test cell grains 5/10/20 km, occurrence thinning and survey-source splits.

**C — cargo and port discrimination:** independently verify port-connected freight routes, passenger-only lines, roads and closed lines. Predict stronger signal for unintentional introductions (e.g. contaminants) than ornamental planting. If there is no reliable freight information, report a yard/port *proxy*, not freight volume.

**D — chronology (conditional):** where actual opening dates and reliable time-structured floras exist, use not-yet-opened/never-opened localities as controls with pretrend and effort diagnostics. N02 snapshots from 2005 are not a historical railway-buildout dataset; first observation is not arrival.

## Interpretation gates
1. Scientific names unified; known alien/native status validated in region. Unknown status excluded, never treated automatically as native.
2. Remove poor coordinates, duplicates, cultivated garden specimens and impossible dates; record source datasets and GBIF DOI.
3. Distinguish passenger railway, freight-only service, yards, abandoned routes; MLIT N02 **excludes freight-only** features.
4. Control road distances, port proximity, urban history, land cover, climate, record effort and spatial autocorrelation.
5. Contrast transport network to similarly routed roads; use spatially and topologically constrained null graphs.
6. Direct seed-transport attribution requires seed collection, genetics/provenance or strong quasi-experiment. Positive observational associations only support **consistency with** transport.

## Decision rule
A promising network signature needs spatial-blocked improvement of M2 over M0/M1 plus robust alien/native contrasts after road, city, climate, port and observation checks. Failure of M2 does not erase separate habitat associations. Never fabricate positive results.
