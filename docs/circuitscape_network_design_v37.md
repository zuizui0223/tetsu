# v3.7 — Circuitscape NETWORK mode for railway-driven plant redistribution

**Status 2026-11-10:** Circuitscape-compatible exact network circuit solver and source-only Overture pilot implemented in the tetsu research branch. This is **not** a claim that Julia Circuitscape has been executed, that the spatial electrical model predicts actual seeds, or that railway plant dispersal has been detected.

## Why circuit theory, and where not to use it

The original v3.0 50km EQUAL-AREA rail-present cells in Tohoku JP, Midwest US and Burgundy FR showed 120/120, 105/105 and 120/120 *mapped physical* rail-connected pairs. Thus a binary connected/unconnected predictor is saturated. In v3.2 source routes, rail path length was nonconstant. Circuitscape network mode, unlike shortest-path routing, admits **all alternative physically connected routes**, producing pairwise effective resistance and, optionally, branch currents. See official [Circuitscape network-mode guide](https://docs.circuitscape.org/Circuitscape.jl/stable/usage/) and [network input format](https://docs.circuitscape.org/Circuitscape.jl/stable/options/).

The source **must be the rail graph**, with source-defined Overture connector IDs representing real physical junctions, NOT a raster in which all rail-side pixels have low resistance regardless of actual tracks, tunnels or missing junctions. This avoids inventing free cross-rail transfers where two drawn railway lines merely overlap on the map.

A meaningful circuit sensitivity is:
- Set each mapped rail segment resistor (r_e propto ell_e), where 1 km of physical segment is assigned a hypothetical 1 resistance unit.
- Solve the graph Laplacian (L=B\,\mathrm{diag}(1/r_e)\,B^T) for 1-unit source and sink current at each mapped source-cell anchor pair. (R_{\mathrm{eff}}(i,j)=V_i-V_j).
- Compare **effective resistance (R_{\mathrm{eff}})** with the equivalent *shortest cumulative rail resistance* (D_{\min}) on the same node pair. For connected pairs, (0<R_{\mathrm{eff}}\le D_{\min}). The dimensionless redundancy (Q=1-R_{\mathrm{eff}}/D_{\min}) is zero on a tree-like single path and increases when alternative parallel routes exist.
- Keep source-cell centroid-to-network snap distances separate and never let a source unit mismatch between meters and kilometers introduce artificial redundancy. Sensitivity: add both endpoint access costs to numerator and denominator (fixed local access branch only).

For source-only transportation networks, a higher modeled current or lower effective resistance is **potential accessibility**, not a measured number of seeds or actual frequency of train/carriage movements. Every edge resistance would ultimately need physically and biologically justified parameters. In real rail networks an isolated branch with no cycle does **not** gain connectivity information from the all-paths circuit calculation beyond its shortest path; test incremental information before investing in a new biological model.

## Current executable implementation

- [Exact sparse grounded-Laplacian network solver](../src/tetsu/circuit_graph.py): preserves parallel-track conductances, prunes dead-end branches and contracts unprotected degree-2 nodes exactly, retaining all 50km focal terminals, then solves per connected component. This is *the same circuit mathematics* used by Circuitscape network mode; it is not a recorded Julia package run.
- [Source data audit](../scripts/audit_overture_50km_circuit_v37.py): directly fetches actual Overture 2026-10 mapped railway segments via source frozen 75km halo and 50km EPSG:6933 equal-area regional anchors (same as v3.2). Compares true effective resistance against Dijkstra's source shortest rail path, outputs counts of pairs with >=1% and >=5% non-shortest-path redundancy, and both fixed centroid anchor choices.
- [Unit tests](../tests/test_circuit_network_v37.py) include tree series (same as shortest path), triangle (effective resistance 2/3 for unit edges), diamond (two alternative 2-unit branches in parallel give 1), two same-endpoint tracks in parallel (0.5), disconnected pairs (null, not zero), and optional Circuitscape textual input format.
- [Actual source CI](../.github/workflows/overture-macro-circuitscape-network-v37.yml): installs SciPy/Overture geospatial stack, verifies the resistor laws, and runs three frozen rail-map blocks. Only aggregate pair metadata is archived; raw Overture line geometry and source connector IDs are not part of the default output artifact.
- CLI: `python scripts/audit_overture_50km_circuit_v37.py --blocks jp_tohoku us_midwest fr_burgundy --output results/audit/50km_rail_circuit_resistance_v37.json`. Optional **local** Circuitscape pairwise network text export: add `--local-circuitscape-export-dir /your/private/source/directory`; this export is deliberately not uploaded by CI. It produces integer-node resistor triples and focal node lists for Circuitscape network pairwise runs.

**Important:** Circuitscape NETWORK format consists of 3 columns `node_id_a, node_id_b, resistance` (one direction per undirected edge, to avoid *double* conductance from mirrored records), plus focal IDs. [Official text list format](https://docs.circuitscape.org/Circuitscape.jl/stable/options/). Physical rail **parallel tracks** with the same source connectors must be summed as conductances, whereas Dijkstra uses the shorter path; conflating them causes systematically false electric-network resistance. The source Overture map includes unknown freight operations and possibly nonoperating physical tracks.

## What would count as a MACRO ecological finding?

1. **Geometry gate:** For the 3 already frozen 50km blocks, quantify (Q) variability across source-cell pairs and test if the circuit result is measurably different from source shortest-path distance. A near-zero (Q) and little among-pair variation means Circuitscape adds almost no ecological predictor value. Record this null result even if disappointing.
2. **Same-scale ROAD negative control:** Repeat the same circuit/source model on real physical road segments and measure the variation beyond road Dijkstra paths. Compare rail vs road **on identical pairs**, not independent country medians, and account for massive city-road graph scale. Road legality/traffic unknown means the network mode is structural only.
3. **Introduced plant response:** Use the already downloaded country-GRIIS introduced-only **speciesKey × 50km cell** matrices. Dependent pairwise species similarities are not 100 independent samples, so estimate at independent country×ecoregion block and use held-out blocks for verification. Adjust at minimum direct geography, ROAD circuit/shortest-path terms, source recording, land and ecology.
4. **Biological parameterization:** Only after independent actual active freight routes, seed contamination/cargo or rolling-stock/seed-trap measurements should conductances be linked to actual propagule transfer, flow direction and introduction propagule pressure.
5. **Falsification:** Permute source rail segment resistances within route class/length strata, and consider "rail ∼ road only", native/unclassified country flora (never automatically infer GRIIS nonmembers are native), and source/target restrictions. No p-values from overlapping pairwise records.

## Model limitations and stopping rule

Electrical random walks are undirected and diffuse across *all feasible physical edges*. Freight trains instead follow scheduled/directed services; rail movement may traverse only a tiny subset of tracks and require station/cargo transfers. (R_{\mathrm{eff}}) is therefore an explorative connectivity index, not a mechanistic seed-flow model.

The 75km graph halo truncates real source long-distance alternatives, potentially overestimating resistance or reporting missing connections. A larger halo sensitivity, actual line availability and documented freight transport are essential before any globally generalized conclusion.

The 2026 transport map postdates both 2010–2017 and 2018–2025 indexed GBIF observation periods. It cannot identify when a plant arrived or when a rail corridor became operational.

**Predeclared decision:** adopt Circuitscape resistance as a secondary network predictor ONLY if the source graph shows meaningful, spatially varying redundancy over shortest-path exposure and the same biological effect improves replication in independent geographic blocks, not because a single city yields the desired positive coefficient.
