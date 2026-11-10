# tetsu v4.0–v4.2 — independent circuit-network controls

**Real source and software verification performed 2026-11-10 (JST).** This is a method and transport-graph structural evidence record. No train-carried seeds, plants adapting to rail wind/mowing, or general positive railway-alien flora association is inferred.

## The ecological question

Does source RAILWAY network electrical redundancy uniquely help explain GRIIS-listed introduced plant community similarity relative to ROAD transport controls, or is high circuit redundancy merely a generic dense transport network property?

Physical Circuitscape NETWORK resistances are stipulated equal to road/rail map edge length in km. They have no biologically measured freight volume, train departures, seed deposition rate or station transfer probabilities. Even a mathematically valid electrical model can be biologically meaningless if its source edge conductances are not identifiable.

## 1. Skeptical map-parallel-link audit — verified source SUCCESS

[Actual Overture three-block rail sensitivity #42775249149](https://github.com/zuizui0223/tetsu/actions/runs/42775249149), three targeted independent synthetic math tests, 3/3 source regions complete and zero errors.

Source assumption A: two Overture segments sharing the SAME two connector IDs have two independent conductances. Source assumption B: conservatively collapse all such parallel source links into only the SHORTEST one, while preserving true alternate loops via different connector IDs. This does *not* prove those original parallel segments were mistaken source records.

| 200km ×200km macro block | Source mapped parallel connector links | Median \(Q=1-R_\mathrm{eff}/D_\min\), sum all parallel conductances | Median Q after one shortest parallel link | Median Q increment from keeping parallel mapped links | Source cell pairs with Q >= 0.05 after collapse |
|---|---:|---:|---:|---:|---:|
| JP Tohoku | 180 | 0.585842 | **0.584302** | +0.001957 | 105/105 |
| US Midwest | 480 | 0.525157 | **0.522844** | +0.001273 | 104/105 |
| FR Burgundy | 633 | 0.588462 | **0.587940** | +0.000598 | 120/120 |

**Interpretation:** The source circuit advantage over shortest paths is **not mainly caused by collapsing parallel links sharing the same mapped connector IDs**. Most mapped circuit redundancy remains even under this conservative source treatment, suggesting distinct alternative mapped routes through different connectors. This does NOT validate the independently operated rail services, historical operating status or actual seed movement along those routes.

Reproduction: [v3.9 physical source audit](../scripts/audit_rail_circuit_parallel_source_v39.py), [parallel-link geometry tests](../tests/test_rail_circuit_parallel_source_v39.py), [actual workflow](../.github/workflows/rail-circuit-parallel-source-v39.yml).

## 2. A genuine ROAD circuit structural comparator — US real source SUCCESS

[US Midwest Overture MAJOR ROAD backbone circuit #42775739386](https://github.com/zuizui0223/tetsu/actions/runs/42775739386), five targeted synthetic checks, completed true source graph:

- Geographic frame: identical outcome-frozen 50km equal-area 16-cell US Midwest 200km block and 75km map halo as the rail source.
- Road classes: **motorway, trunk, primary, secondary ONLY**, a deliberately fixed *major-road backbone* designed to fit sparse graph resources, **NOT** the entire Overture road network (the previous source ROAD shortest paths also included tertiary, residential, unclassified).
- Actual real Overture eligible segments: **137,207**; all 16/16 source cells with mapped major-road connectors; **120/120** source cell pairs have connected physical major-road paths.
- True physical circuit solve: original source graph reduced by series/deadend elimination to **31,792 nodes**, 35,000-node hard solver bound; numerical solver direct sparse grounded-Laplacian.
- Median \(Q_{\mathrm{road\,backbone}}=1-R_\mathrm{eff}/D_\min=\mathbf{0.867229}\), compared with \(Q_\mathrm{rail}=\mathbf{0.525157}\) in the geographically same broad source block, although endpoints/classes and networks are not directly interchangeable.

**Interpretation:** ROAD infrastructure also has strong alternative-path electrical redundancy. Circuit redundancy is therefore NOT a distinctive indicator of rail-mediated plant dispersal by itself. It would be erroneous to say roads move seeds more than rails from these abstract resistor values. Source maps use different classes, network densities and adjacent anchors, and actual train/vehicle movement/transport frequencies remain unknown.

The French Burgundy independent source run is a separate workflow [major road circuit comparison](../.github/workflows/road-backbone-circuit-fr-v40.yml). No conclusion from its output until a successful live source run and residual/size gates.

Reproduction: [road circuit class-restricted source](../scripts/audit_overture_50km_road_backbone_circuit_v40.py), [offline source tests](../tests/test_road_backbone_circuit_v40.py).

## 3. ACTUAL Julia Circuitscape executable agrees with SciPy small-network reference

[Original official Julia package CI #42775947552](https://github.com/zuizui0223/tetsu/actions/runs/42775947552), SUCCESS, Julia Circuitscape.jl installed and ACTUALLY EXECUTED on three synthetic resistor NETWORK cases. The tetsu adapter exported the 3-column resistor edge file / focal IDs / pairwise INI, the Julia package computed a resistance matrix, and the result was compared with the independently solved sparse SciPy network.

| Genuine executed synthetic network | Independently verified focal-pair values | Max absolute Julia vs SciPy difference |
|---|---:|---:|
| Triangle unit resistors | 3 | \(1\\times10^{-10}\) |
| Diamond two alternative 2-unit paths | 1 | 0 |
| Series 1+2 resistor chain | 1 | 0 |

**Conclusion:** The network format is truly executable by Julia and the two mathematical implementations agree on tested cases. This is **not** a claim that an entire mapped 8,581–41,380-node source railway network has been recomputed in Julia, nor that conductance parameters represent seed dispersal.

Source adapter [SciPy physical network](../src/tetsu/circuit_graph.py), [actual Julia executable input generator and crosscheck](../scripts/verify_actual_julia_circuitscape_v41.py), [Julia driver](../scripts/verify_actual_julia_circuitscape_v41.jl), [real workflow](../.github/workflows/actual-julia-circuitscape-v41.yml). Official [Circuitscape.jl documentation](https://docs.circuitscape.org/Circuitscape.jl/latest/options/).

## 4. Identifiability and next check

Previously, on two broad GRIIS-plants US/FR blocks, 2018–2025 flora similarity partial rank versus log(RAIL/ROAD shortest route) was **−0.264 (US) vs +0.054 (FR)**. Substituting physical RAIL circuit effective resistance changed that to **−0.073 (US) vs +0.016 (FR)**. Independent geography is only two blocks; the direction and magnitude do not demonstrate a replicated positive rail-flora effect.

Before a meaningful biological network analysis, source-equivalent ROAD circuit evidence should be joined as a partial transport confound, accompanied by actual active rail/freight routes, complete DOI-backed GBIF individual event records, observer effort and coordinate accuracy, and independent country×ecoregion geography across the existing eight-country GRIIS source frame. Never treat spatial cell pairs as independent geographic replicates or tune resistances/matching to find significance.

### Research decision

The Circuitscape-style physical network indicator is worth retaining as a **secondary railway network topology axis**, because it has measurable structural variation robust to mapped-parallel-link collapse, AND because independently running Julia Circuitscape validated the exported network format. But it has not earned a causal botanical interpretation or priority over the shortest rail/road route until independent ecological replications and biologically supported conductances exist.
