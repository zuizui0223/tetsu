# Rail-status evidence update (2026-11-25)

Verified completed runs: cross-taxon 45560367041 (artifacts 14998137806, 14998128170); population/local 45561352038 (artifact 14997871993). Both successful, tests passed. Original 288,721,938 rows; eligible plant 7,827,936, bird 4,431,373, insect 4,960,614, fungal 585,925 observations. Exact v369 plant replay/hash gates passed.

## Cross-taxon, 50% land, full status conditional
Relative heldout MSE gains in matched cells (plant, comparison):
- Birds: 100km n=2659 (+1.860%, -1.522%), 200km n=1222 (+2.265%, -0.282%). Plant 0/39 fake maps >= real at both widths. Paired active slope differences 0.00582 [-0.00870,0.02034], 0.00405 [-0.01072,0.01883].
- Insects: n=2519 (+0.0063%, -0.639%), n=1126 (+0.611%, -0.326%). Paired differences -0.00030 [-0.01302,0.01241], 0.00991 [-0.00545,0.02527].
- Fungi: n=1030 (+5.458%, -4.247%), n=475 (+1.576%, -0.674%); only 2 eligible continent folds, not global heldout validation. All pairwise slope difference intervals overlap zero across all panels.
At 80% land, plant vs bird gains remain positive but shrink at 200km (+0.046%); insect-matched plant gains are negative at both widths.

## Population/local plant model
50% land 100km n=3335, 175 blocks: active 0.00534 [-0.00165,0.01233]; inactive -0.00491 [-0.01313,0.00330]; difference 0.01026 [-0.00141,0.02192].
200km n=1370, 159 blocks: active 0.01154 [-0.00006,0.02313]; inactive 0.00419 [-0.00922,0.01760]; difference 0.00735 [-0.01335,0.02805].
Population-adjusted conditional active relative MSE gain 0.770%, 0.906%; at 80% land -0.024%, -0.230%. All primary local status difference intervals include zero; 200km 80% land high-effort and both-positive subset difference estimates are negative.
Conclusion: weak, geographically sensitive predictive evidence for plant-specific rail-status association, not a demonstrated difference in coefficients or a causal operating-train mechanism. Population analysis and cross-taxon analysis were separate. Next combined sensitivity uses fixed common cells, population, local block FE and direct paired double contrast (plant active-minus-inactive slope minus comparison-taxon active-minus-inactive slope). Bird is primary based on prior results, not independent confirmation. No further new raw-data scan.
