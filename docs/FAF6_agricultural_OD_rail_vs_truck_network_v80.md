# tetsu v8.0 — Actual FAF6 state OD freight network: agricultural RAIL versus TRUCK, no seed-transport causality

**Source-verified 2026-11-12:** [REAL GitHub Actions SUCCESS #43142647675](https://github.com/zuizui0223/tetsu/actions/runs/43142647675), [original JSON artifact #13910351407](https://github.com/zuizui0223/tetsu/actions/runs/43142647675/artifacts/13910351407), [reproduction code](../scripts/analyze_FAF6_crop_OD_rail_vs_truck_network_v80.py), [six synthetic source/math integrity tests](../tests/test_FAF6_crop_OD_rail_vs_truck_network_v80.py).

## Genuine macro research result

2022 U.S. Bureau of Transportation Statistics FAF6.0 state shipment origin–destination records: 599,529 official commodity/mode/origin/destination rows, 51 state-equivalent source codes, 2022 thousand U.S. short tons (estimated). The original 2022 ZIP checksum, exact 599,529 rows, original 51 states, and nationwide crop and resource RAIL/TRUCK total tonnage were crosschecked **exactly** against the separately verified v7.9 artifacts. This is not a rail-network route assignment.

The prespecified potential plant-propagule commodity groups are **SCTG02 grains (may include seed)**, **SCTG03 other agricultural products**, their SCTG02+03 aggregate, and **SCTG10–12 resources** as a different-category comparison. Neither the commodity label nor movement volume verifies viable exotic seeds. RAIL and TRUCK modes are compared from the same geographic economic source.

### Actual directional state-scale macro network

| 2022 domestic U.S. agricultural commodities SCTG02+03 | RAIL | TRUCK |
|---|---:|---:|
| Estimated domestic tonnage, million U.S. short tons | 128.799 | 1,155.992 |
| **Share of each mode's own tonnage crossing state boundaries** | **63.4285%** | **17.7598%** |
| Positive directed *interstate economic OD* combinations | **155** | **1,347** |
| State source codes with positive interstate outgoing flow | 32 | 49 |
| State code destinations with positive interstate incoming flow | 36 | 49 |
| Median effective destination partner states (Shannon count, origin states with outgoing flow) | **1.789** | **5.484** |
| Nonzero unordered interstate state pairs | 147 | 846 |
| Both opposite directions positive among those unordered pairs | **8** | **501** |
| Weighted interstate directional imbalance | **0.9833** | **0.5190** |

**Same-origin paired comparison:** 32 state codes have positive interstate agricultural outgoing flows in BOTH modes. Rail has smaller Shannon-effective interstate destination number in **29/32**, greater in 3/32; equal-state medians are **1.789 RAIL** versus **5.051 TRUCK**. Rail largest-destination fraction is higher in 25/32, smaller in 7/32; median share **84.10%** RAIL versus **42.05%** TRUCK. These are 32 economic state-origin descriptive comparisons, NOT 32 independent biological experiments or train-seed events.

**Mode overlap:** of 155 RAIL positive interstate OD directions, **130** have nonzero TRUCK freight for identical commodities; **25** have no positive TRUCK modeled flow. Those 25 account for **18.37%** of rail's interstate agricultural tonnage. A zero modelled transport mode is not verified vehicle absence.

**Within- vs interstate:** rail makes up **10.02%** of RAIL+TRUCK domestic agricultural tonnage when intrastate movements are included, but **28.47%** of RAIL+TRUCK **interstate** agricultural tonnage (different denominators). Thus railway network opportunity is unusually interstate-oriented in this economic source.

### Commodity sensitivity: not just grains or only special crop coding

| 2022 domestic fraction of each mode's freight crossing state boundaries | RAIL | TRUCK |
|---|---:|---:|
| SCTG02 cereal grains | **61.92%** | 9.99% |
| SCTG03 other agricultural commodities | **70.38%** | 29.12% |
| SCTG02+03 combined | **63.43%** | 17.76% |
| SCTG10–12 resource comparison | **14.35%** | 5.92% |

The same **rail-greater-than-truck interstate tendency also exists for resource commodities**, so it is NOT specific biological evidence of seed dispersal. The larger 2022 rail-vs-truck interstate gap for crop commodities may instead be economics, cargo type, geography or cost.

### Positive model SOURCE finding, not a verified biological railway effect

This is a real, repeatedly verified **macro transport opportunity structure**: compared with truck freight, 2022 U.S. agricultural RAIL shipments were more heavily interstate, more concentrated into a small number of modeled origin→destination combinations and more directionally asymmetric. Such a pattern generates a **mechanistic hypothesis of focused propagule arrival opportunities at certain freight destinations**, not evidence for actual viable weed propagules.

The v7.4–v7.8 plant results remain source/dataset and metric fragile; they do NOT establish freight-mediated plant introduction, convergence, mowing tolerance, train-gust selection, plasticity or evolved adaptation. Do not portray the v8.0 state economic network as 50km Overture physical railway traffic: this data does not tell which track segment carried freight, when any seed dropped, or what germinated.

## Next source-identified macroecology unit

Obtain independently reproducible state×year country-GRIIS plant recorded species categories via official Census FIPS state names + GBIF stateProvince facets, with explicit observer/source effort and taxonomy QA. Only after that can state-to-state rail and truck agricultural OD links be joined to species assemblage data for exploratory **source-only** association; biological effects require multiple years matched to freight operation, independently observed seed contamination/arrival and trackside/roadside plots.

Freight and filtering hypotheses must stay separated:
- Delivery: origin agricultural species source pool × directional FAF6 RAIL/TRUCK mode-specific OD estimated tonnage, with alternative route allocation and potential cargo contamination measured independently.
- Filtering: destination trackside vegetation × documented cutting/gust/ballast/heat and measured species regenerative/dispersal traits.
- Adaptation: replicated same-species phenotype and genetic/common-garden evidence, not GIFT average growth form.

No cell-pair p-values and no train-carried seed quantities reported.
