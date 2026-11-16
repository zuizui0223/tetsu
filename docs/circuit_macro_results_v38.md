# tetsu v3.8 — Circuit-theory railway network: source novelty versus biological test

**2026-11-10 verified real-source research update.**

## What is computed, and what is not

Circuitscape official documentation supports **network input** and **pairwise effective resistance** from a resistor list, in addition to raster mode. The current tetsu implementation uses a mathematically equivalent SciPy sparse grounded-Laplacian solver (NOT the Julia Circuitscape executable) and provides optional local Circuitscape-compatible network triples / focal node text lists. See [official network mode](https://docs.circuitscape.org/Circuitscape.jl/stable/usage/) and [3-column input format](https://docs.circuitscape.org/Circuitscape.jl/stable/options/).

Original source is 2026-10 Overture railway edges with explicit connector IDs. Each link's **hypothetical resistance** is equal to its length in km-equivalent units. Parallel track links sharing source connectors are combined by summing conductances; shortest paths continue to use the minimum edge resistance. This is required mathematics, but is NOT direct observation of train seed movement or freight/service frequency. Overlapping drawn rails without a common connector ID are never connected.

The source true physical graph includes current, potentially inactive rail infrastructure; 75km graph halo, centroid anchors and resistance choice remain uncalibrated to plant dispersal.

## 1. Source-only circuit geometry — positive technical result

[Actual Overture three macro block source run #42766118921](https://github.com/zuizui0223/tetsu/actions/runs/42766118921) **SUCCESS**, all nine synthetic analytic network checks, all three pre-frozen 50km source blocks, 0 actual source errors. Tests include tree, triangle, diamond, parallel paths, disconnected nodes, Circuitscape text output and physical source unit consistency. A discovered critical unit mismatch was corrected *before* reporting the real-source results: Dijkstra's path is returned in meters, whereas normalized circuit resistance is expressed in km equivalents. A physical (R_{\mathrm{eff}} \le D_{\min}) bound is now enforced for every source pair.

| Region (200km box, 50km cells) | Original Overture rail connector nodes | Physical mapped multigraph cycle rank | Source-cell pairs reachable in circuit | Pair count with circuit redundancy >=5% | Median (R_{\mathrm{eff}}/D_{\min}) | Median (Q=1-R_{\mathrm{eff}}/D_{\min}) |
|---|---:|---:|---:|---:|---:|---:|
| JP Tohoku | 8,581 | 1,211 | 105 | **105** | **0.414158** | **0.585842** |
| US Midwest | 41,380 | 5,817 | 105 | **104** | **0.474843** | **0.525157** |
| FR Burgundy | 40,075 | 7,671 | 120 | **120** | **0.411538** | **0.588462** |

This is a **real positive result for measuring network structural redundancy**: a circuit effective resistance often differs substantially from shortest path, not only from a binary connected flag. The multigraph cycle rank includes parallel links as cycles; Overture duplication and link representation must be audited before assigning a biological role to these alternatives. Trees yield no circuit-vs-shortest novelty (synthetic test passes). Current rail corridors do show non-tree source alternatives.

No ecological plant data was used to choose or calculate these numbers. The reported (Q) is an index of *hypothetical physical alternative path conductance*, not verified parallel train services or seed delivery.

## 2. Actual macro ecological comparison — not replicated positive

[Actual four-source US/FR plant-vs-circuit incremental run #42767060549](https://github.com/zuizui0223/tetsu/actions/runs/42767060549) **SUCCESS**, all five synthetic quality checks, and the same frozen 50km source cells, GRIIS-introduced plants, physical rail shortest paths, road shortest paths and physical rail effective resistance. Road network and plant data are on the *same source cell pairs*, rather than comparing unpaired medians.

**Outcome:** GBIF-indexed GRIIS-listed introduced species-key similarity, rarefied to **20 positive indexed records per source cell**, using the *ratio of expected species-key intersection to expected union* (not (E[Jaccard])).

**Adjustment:** partial ranks of geography, physical ROAD shortest path, source GBIF indexed record counts and count balance. For *redundancy at a fixed shortest path* additionally adjust the RAIL shortest path. Each cell occurs in multiple pairs; **no pairwise p-values or inflated confidence intervals**.

| Region | Period | n dependent source cell pairs | Partial rho of shortest RAIL/ROAD route vs source flora similarity | Partial rho of circuit RAIL/ROAD effective route vs similarity | Partial rho of extra path redundancy (Q), controlling shortest RAIL |
|---|---|---:|---:|---:|---:|
| US Midwest | 2010–17 | 91 | −0.044934 | +0.095766 | −0.230679 |
| US Midwest | 2018–25 | 91 | **−0.264201** | **−0.073106** | **−0.181667** |
| FR Burgundy | 2010–17 | 120 | +0.148512 | +0.295627 | −0.231848 |
| FR Burgundy | 2018–25 | 120 | **+0.053809** | **+0.016237** | **+0.066599** |

**Interpretation of sign:**
- For route-ratio correlations, **negative** means shorter RAIL relative to ROAD aligns with more similar recorded introduced flora.
- For redundancy (Q), **positive** means more alternative physical RAIL paths align with more similar recorded introduced flora, conditional on the shortest RAIL route.
- In the later period, the US circuit route association is **weaker** in absolute magnitude than the source shortest route (−0.073 vs −0.264). France is approximately zero under both circuit and shortest (+0.016 vs +0.054).
- The specific redundancy (Q) association has **opposite signs** in the later period: US −0.182, FR +0.067. The early period is negative in both blocks but lacks independent replication, and direction changes within France between periods.

This is not a generalized positive signal supporting railway-driven introduced plant homogenization. Circuitscape network equations are **useful structural measurement tools here**, but **do not yet improve biological explanatory evidence**. Differences between rarefied and raw observed-key Jaccard metrics are reported in the source output; neither is systematic occupancy data.

Both predeclared source-anchor strategies (nearest within a cell, nearest on largest physical component) were evaluated and returned matching summary correlations in these two regions. No selection for a desirable direction occurred.

## 3. Why no causal seed-transport interpretation

- Overture physical rail lines can be inactive, closed, nonfreight, or without actual train stops, and are only observed in a 2026 map, not over the 2010–25 flora sampling periods.
- A railway can be a habitat and/or propagule delivery route, and both correlate with ROAD density, freight ports, municipal disturbance, accessible plant observations and climate.
- Under an undirected resistor model, seed dispersal is analogized to a random walk across all physical tracks. Real cargo and rolling stock move directionally and by service schedule; equal conductance per km is only a structural assumption.
- The sample of **two independent 200km blocks** is not sufficient for a biogeographic macro significance claim, even with 211 pairwise cell comparisons (those pairs share source cells).
- GRIIS country status varies in completeness, and all source plant observations are GBIF faceted positive indexed occurrence keys, not DOI-downloaded event-level surveys, without observed absences, geographic coordinate uncertainty checks or observer effort.
- No comparison to actual *ROAD circuit effective resistance* yet, only physical road shortest path. A cross-network circuit index is not fully matched until both rails and roads use comparable resistance/availability parameterization.

## 4. Decision and next evidence gate

**Keep Circuitscape NETWORK as a secondary physical connectivity axis, not the current manuscript's mechanistic centerpiece.** It demonstrated measurable structural redundancy that shortest path misses, but has not added replicable botanical evidence across these two broad blocks. Next advance:
1. Determine whether source parallel conductance comes from genuine independent tracks/routes vs Overture mapping duplication; carry track-class and operational status when available. Report an alternative resistance model with unit links and/or validated train frequencies, never tune resistances to favor a plant coefficient.
2. Make source ROAD circuit effective resistances using comparable independently validated network equations; verify the extra rail circuit axis is not merely a universal infrastructure density/cycle measurement.
3. Expand the plant and continuous infrastructure join from 2 to geographically independent blocks, using the already validated [v3.6 8-country 105-cell source-eligible introduced flora](https://github.com/zuizui0223/tetsu/actions/runs/42761140912). Prioritize independent blocks/true country–ecoregion coverage rather than multiplying non-independent within-block cell pairs.
4. Obtain historical/active freight corridor evidence, source observer/georeference quality, and physical seed vector data before ever arguing propagule transport causation.

## Reproducibility

- [Sparse network circuit computation, parallel-track preservation, series contraction and optional official-format export](../src/tetsu/circuit_graph.py)
- [Actual source physical railway resistance audit](../scripts/audit_overture_50km_circuit_v37.py)
- [Nine analytic tests](../tests/test_circuit_network_v37.py)
- [Real 3-region source workflow](../.github/workflows/overture-macro-circuitscape-network-v37.yml)
- [Real 2-block biological incremental diagnostic](../scripts/analyze_macro_alien_circuit_increment_v38.py)
- [Five synthetic-source tests](../tests/test_macro_alien_circuit_increment_v38.py)
- [Real four-source posthoc workflow](../.github/workflows/circuit-increment-plants-v38.yml)
- [Circuitscape-network interpretation and source parameterization](circuitscape_network_design_v37.md)

**Note on software:** The exact grounded-Laplacian solver implements network electrical effective resistance in SciPy. An actual Julia Circuitscape executable comparison is still pending; the tests validate the mathematical identities, not a claim of full software numerical identity on every source track graph.
