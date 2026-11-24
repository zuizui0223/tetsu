# v316: strict incidence and trait-consistency audit (2026-11-22)

The previous v315 matched-site preflight has been locally reviewed. Two failure modes were identified: (1) the same species could have different temperature indicators across sites, and (2) a numeric abundance/cover value could be treated as incidence. Either could bias the thermal-affinity contrast.

A new local v316 implementation and eight synthetic tests were executed successfully (8/8): valid paired data; missing-trait coverage; inconsistent taxon trait; abundance-as-incidence; mismatched within-site climate; duplicate taxa; unpaired data; invalid coverage threshold. The v316 code is available as a conversation artifact and **is not yet committed in this PR**. The previous v315 source code remains unchanged.

**No source-backed climate model was fitted.** The scientific bottleneck remains genuine local rail/nonrail matched vegetation data in multiple independent climatic regions. The 19 Slovak stations do not establish independent climatic replication; Slovenian shoulder vs toe is a within-railway comparison.

Implementation specification for source replay: input plot-level binary species incidence (0/1), explicit study/region/site/habitat, fixed temperature indicator and trait provenance for each taxon within study-region, no conversion of 'x' to zero, separate descriptive delta and eligibility flag, sensitivity thresholds 0.50/0.70/0.90, and block/region validation before fitting a climate interaction.
