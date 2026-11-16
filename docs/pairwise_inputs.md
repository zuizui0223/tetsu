# Descriptive pilot CSV contract

**This pilot calculates similarity, not train-mediated dispersal or causal effects.** Example data are synthetic.

1. `incidence.csv` columns: `region_id,taxon,status`. `status` must be exactly `alien` or `native`, and be vetted for the study region. Unknown taxa are not native by default. At most one row for each region/taxon. For a Japan-only study, the same taxon has one national status; international studies need region-specific status.
2. `pairs.csv` columns: `region_a,region_b,euclid_km,road_km,rail_km,climate_distance,urban_difference,effort_a,effort_b`. Distances in km, `rail_km` blank if the places have no active network connection. Effort is independent total botanical sampling opportunity. Use each unordered pair once. A rail line passing near two regions is NOT sufficient to mark them connected.
3. `scripts/run_pairwise.py` creates `pairwise.csv` (Jaccard alien and native overlaps and their difference, connectedness and simple diagnostics) and `summary.json`. **No p-values**, since pair comparisons are dependent and network connectivity is observational.

The provisional `rail_directness` diagnostic is log(1 + road_km) - log(1 + rail_km) for connected regions, zero for disconnected. It is not measured train flow, does not imply seed transport, and is NOT used by the summary as proof.
