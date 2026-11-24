# v363: worldwide physical rail GERS connector position identity QA (source topology validation)

2026-11-24. Frozen AFTER the v361 global giant-component result and BEFORE v363 positional source QA, so exploratory audit, not biological preregistration.

## Why
v361 source 2,054,855 rail segments, 4,970,641 connector UUIDs; largest physical component ~1,548,770 km (all original rail) and ~1,369,687 km (conventional active-candidate), with most rail centerline length in top 10 components. A wrongly shared connector UUID far apart could artificially merge unrelated subnetworks.

## Verified source inputs and fixed sampling
Original global source Overture 2026-11-18.0 in 16 genuine rail-only GeoParquet artifacts from complete all-128-source original Actions run 45342573085; source original all-world rail segment count 2,054,855. Enforce artifact SHA checks against genuine complete original source census before any result. Source license Overture Maps / OpenStreetMap contributors (ODbL).

GERS connectors[].connector_id UUIDs are physical graph nodes; connectors[].at are expected normalized WGS84 LineString positions 0–1. For each original connector UUID select deterministically if SHA256(first 32 bits) modulo **256 == 0** (about 0.39% of all original connector UUIDs), independent of location, plant identity or observed mismatch. For every source occurrence of a selected connector, interpolate original source line geometry at the source normalized fraction; compute WGS84 geodesic separation from FIRST source occurrence of the same connector UUID. Report n sampled UUID, n multi-incidence connector UUID, maximum per-UUID source point discrepancy and bins >5, >25, >100, >1000 metres, including exact same source segment duplicate anomaly counts. Original source missing or invalid normalized positions audited and not silently zero.

Interpretation: a low >25m mismatch rate supports source graph node positional consistency but does not certify an actual active switching path or name of a railway route. Any >100m far-location source connector identities require extra graph component sensitivity excluding those UUIDs before a robust giant-network statement; no retrospective deletion to strengthen source claims. Overture source sample remains all original rail classes; conventional class analysis follows as separate class sensitivity.

Publish only aggregate hashed-source counts and distance bins, no original connector IDs, precise source railway geometries, source rare botanical GPS, plant names or raw source species rows. No botanical species Jaccard or train dispersal causality estimated.
