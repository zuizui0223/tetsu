# tetsu v7.9 — Genuine 2022 FAF6 agricultural commodity RAIL vs TRUCK: macro freight opportunity, NOT observed seed dispersal

**Real data verified 2026-11-12:** [source workflow SUCCESS #43128165864](https://github.com/zuizui0223/tetsu/actions/runs/43128165864), [source JSON artifact #13903049131](https://github.com/zuizui0223/tetsu/actions/runs/43128165864/artifacts/13903049131), [auditable code](../scripts/audit_BTS_FAF6_2022_agricultural_rail_freight_v79.py), [5/5 source tests](../tests/test_BTS_FAF6_2022_agricultural_rail_freight_v79.py).

Original source: US Bureau of Transportation Statistics (BTS), [FAF6.0 state OD 2022](https://www.bts.gov/faf), [2022 State CSV ZIP](https://faf.ornl.gov/faf6/data/Download_Files/ZIP/FAF6.0_State.zip), [FAF6 DOI 10.21949/fhgm-g335](https://doi.org/10.21949/fhgm-g335). Original source ZIP SHA256 is in the machine-readable artifact, as are full source field/row conservation measures.

## Real original data and provenance

- **599,529 original 2022 state-scale commodity × OD × transport mode records**, spanning **51 state-source code categories**. Five targeted tests passed.
- Actual FAF6.0 STATE original columns: **dms_origst, dms_destst, dms_mode, trade_type, sctg2, tons_2022**. A preliminary reader assumed the FAF regional rather than state schema and correctly refused to report results; its source-header bug was repaired before this verified run.
- **RAIL dms_mode=2; TRUCK dms_mode=1.** The source unit is thousand US SHORT TONS of estimated commodity flow, **not** ton-miles, individual rail carloads or viable biological seeds.
- Pre-frozen **SCTG02 (Cereal Grains including seed)** and **SCTG03 (Other Agricultural Products)** as broad potentially propagule-associated commodities before assessing any railway/plants biological effect. SCTG10–12 is another resource category comparator, not a presumed sterile negative control.
- Output includes domestic, domestic portions of imports and domestic portions of exports separately; interstate + within-state domestic OD flow is not counted as distinct plants or observed transit on any specific railroad link.

## Aggregate 2022 potential agricultural commodity transport by mode

| SCTG02+03 RAIL vs TRUCK only | RAIL, million US short tons | TRUCK, million US short tons | RAIL / (RAIL+TRUCK) |
|---|---:|---:|---:|
| **Domestic-only movement** | **128.799** | **1,155.992** | **10.0249%** |
| US domestic transportation portion of imports | 6.052 | 29.104 | 17.2146% |
| US domestic transportation portion of exports | 55.928 | 65.992 | **45.8727%** |

Other domestic source benchmarks: resource code 10–12 rail fraction (among RAIL + TRUCK) 8.7795%, all OTHER codes 11.6194%. These are not biologically independent negative controls.

Domestically, SCTG02 cereal grain freight: **105.832** million US short tons by rail versus **686.194** million by truck. SCTG03 other agricultural freight: **22.967** million rail versus **469.798** million truck. **38** origin-state code categories have positive SCTG02/03 estimated tonnage for BOTH modes; do not interpret these 38 as independent biological replicates.

This is genuine macro commodity movement evidence that RAIL can carry crop-linked material in large estimated tonnage quantities. It is **not** evidence that non-native seeds were carried or deposited. Cereal grain itself is not necessarily an alien plant propagule, and plants observed beside tracks need not have arrived by freight. The 2022 regional economic OD tonnage is NOT a mapped railway segment’s actual train traffic.

## Source scales and next biological identification

Keep **FAF6 2022 state or FAF economic OD zones** separate from the earlier **2026 FRA road–rail CROSSING count sample** and from **2018–2025 GBIF species records at 50km grid cells**. State/region OD freight cannot be assigned to the Overture shortest mapped line as proven path throughput. Such route allocation would require an explicit, independently validated traffic assignment model with uncertainty.

A genuine macroecology test of seed supply should join:
1. The pre-frozen commodity-specific origin/destination rail-vs-road *potential source pressure* from FAF (for the 2022 base year, with historical periods from comparable releases only).
2. Independently documented ACTIVE freight/service corridors, and documented cargo/vehicle seed contamination or dispersal traits, to identify actual carriage rather than a speculative path.
3. Preselected independent country×ecoregion flora, protected from GBIF dataset/provider/observer bias, with real temporal matching and road, urban, climate and vegetation controls.
4. INDEPENDENT mowing, weed-management and gust-wind exposure observations plus actual species dispersal and regeneration traits. GIFT growth form alone is not mowing tolerance or wind resistance. Phenotypic plasticity/evolved adaptation requires phenotypic repeated measures and heritable/common-garden/genetic evidence.

A main transport-by-stress interaction is **not estimable yet**. Current v7.8 apparent negative woody richness correlation survives strict positive-train restriction in only 1 Midwest block of the original 3 US geographic windows. Southeast and Pacific fail sample support; wider freight tonnage availability does NOT convert the old 50km source associations into a replicated biological effect.

**Manuscript-safe outcome:** robust new independent freight vector-OPPORTUNITY source across dozens of economic region codes, not railway-mediated actual seed carriage, plant establishment, grass-cutting/wind adaptation, plasticity or natural selection.
