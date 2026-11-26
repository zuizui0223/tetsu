# Verified selection sensitivity and manuscript decision — 2026-11-26

## Receipt
Workflow [45758090537](https://github.com/zuizui0223/tetsu/actions/runs/45758090537), SUCCESS; source commit 94e0d6f7c8c797e9d486e1f79fec2695e122fc28. Artifact 15088069603, archive digest sha256:5dfbaeb8c85c9114c471878202c126a316615ffe6b950596558ec0cc07c13b90 (GitHub metadata; archive not independently downloaded here). Previous support and native coefficient replays passed. Results transcribed programmatically from the four successful job-log JSON panels. Native and IPW effects below use the same selected overlap sample; approximate intervals conditional on estimated weights.

| km | land | plant–bird n | three-taxon n | plant slope IPW | insect slope IPW (conditional on plant) | plant coefficient attenuation IPW |
|---:|---:|---:|---:|---|---|---|
| 100 | 0.5 | 1769 | 1576 | 0.21126 [0.10366, 0.31887] | 0.10520 [0.03646, 0.17393] | 0.03430 [0.01170, 0.05690] |
| 100 | 0.8 | 1500 | 1340 | 0.16826 [0.10074, 0.23579] | 0.06752 [-0.00476, 0.13980] | 0.02368 [0.00128, 0.04608] |
| 200 | 0.5 | 933 | 811 | 0.25237 [0.16909, 0.33565] | 0.13906 [0.07053, 0.20759] | 0.04473 [0.01863, 0.07083] |
| 200 | 0.8 | 753 | 653 | 0.26956 [0.17281, 0.36632] | 0.15871 [0.08824, 0.22919] | 0.03977 [0.01417, 0.06536] |

| km | land | plant–bird weighted frame overlap | three-taxon weighted frame overlap | plant–bird ESS | three-taxon ESS | zero-rail plant–bird cells | zero-rail three-taxon cells |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 0.5 | 38.1% | 32.8% | 229.1 | 201.3 | 70 | 35 |
| 100 | 0.8 | 36.6% | 32.7% | 315.2 | 234.5 | 37 | 19 |
| 200 | 0.5 | 61.4% | 54.1% | 342.4 | 280.7 | 42 | 17 |
| 200 | 0.8 | 57.5% | 52.3% | 260.9 | 203.2 | 17 | 9 |

Primary 100km land80, paired weighting change in two-taxon plant slope +.04210 [-.00178,.08598]. In three-taxon model, plant slope change +.05728 [.00984,.10472]; insect -.08937 [-.15684,-.02190]; attenuation -.02302 [-.04523,-.00081]. At 200km land80 insect change -.02968 [-.06600,.00664], attenuation -.01002 [-.02372,.00367], plant slope change +.04506 [.00489,.08524]. Two-taxon spatial association remains positive in every weighted panel. Three-taxon insect association weaker under weighting, 100km land80 interval crosses zero, 200km both land thresholds remain positive. Attenuation is positive in four panels, not mediated causal proportion. No statistically established increase of the two-taxon slope due to weighting. Overlap support is extremely incomplete, especially all-four-zero: only 37 and 19 zero-status cells in primary two- and three-taxon overlap. Not evidence of representative global ecology or railway-specific coupling.

## Evidence hierarchy and stop rules

1. Global mapped railway association with observation support: robust descriptive selection pattern; independent environmental frame, stage analysis [45754547988](https://github.com/zuizui0223/tetsu/actions/runs/45754547988), 100km land80 joint effort active-minus-inactive -.02051 [-.03450,-.00652]; all four panels negative. Only two of 2320 jointly effort-sufficient cells fail species count gate. Inactive status predicts effort support on heldout continents; not causal historical effect, possible mapping/access confounding. Distinguish unadjusted active-versus-zero and conditional continuous status slopes. 
2. Spatial plant–bird observed richness association persists across observer sides, AGB, environment, time-specific panels and measured-selection reweighting, within overlap only. Temporal within-cell coupling [45718095925](https://github.com/zuizui0223/tetsu/actions/runs/45718095925) cross .00891 [-.04326,.06108] vs same .14556 [.08274,.20837], not a longitudinal ecological response. Composition [45716992968](https://github.com/zuizui0223/tetsu/actions/runs/45716992968) retrospective primary uncertain. Do not state direct interactions or true richness.
3. Insect contribution scale/selection dependent; broad fruit+seed dietary contrast unsupported, granivore-only exploratory. Retrospective heldback annual insect-guild test was frozen before earlier insect data in 141391523fc66e957de16a745a288c684ec0ee6a. Run 45758161138 failed at test import due to unmatched brace; repaired only test fixture at 20304cc5f7555fa681ee2e993ae8ad5152072400; rerun [45759346784](https://github.com/zuizui0223/tetsu/actions/runs/45759346784). No numerical result yet. If frozen test unsupported, close this same-source diet-specificity branch rather than further subdivision.
4. Railway-specific ecological connectivity not established. Strict nonrail beta comparator lacks common support; continuous railway moderation not robust; independent eBird primary not replicated (not a source for further work). Do not treat non-significance or inadequate support as evidence of no effect. Regional plot evidence remains separate from global coarse records. 
5. Strongest integrated narrative: mapped infrastructure and selection into citizen-science samples; observed cross-taxon spatial structure within the supported domain. An ecological railway mechanism manuscript needs a genuinely independent observation design, finer-scale surveyed track/control plots or a prospective actual-operation contrast. Current observational dataset cannot by itself identify train dispersal, food-web mediation, or historical railway effects.

This document is a result audit, not a new model. No source thresholds, main branch, or previous analysis files changed. No new eBird usage.
