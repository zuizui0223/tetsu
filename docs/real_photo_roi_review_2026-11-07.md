# First real-photo ROI feasibility review — 2026-11-07

**This is a provisional single-reviewer visual screen at 280px contact-sheet scale, NOT a calibrated segmentation study.** A second blinded reviewer and full-resolution photo quality assessment are required before publication.

Original source-backed photo cohort: [GitHub Actions #42373519458](https://github.com/zuizui0223/tetsu/actions/runs/42373519458). The artifact includes `contact_sheet.jpg`, 13 CC-BY thumbnails, `photo_attribution.csv` with original photo URLs and credited authors, and `summary.json`. Images are not committed to the source repository. The deterministic sampler is [`scripts/sample_photo_roi_candidates.py`](../scripts/sample_photo_roi_candidates.py).

## Empirical data availability

Original sampling frame: all quality-grade research photos of `Trifolium repens` and `Solidago altissima` within **fixed** Tokyo/Osaka rectangles, 2018–2025, independent of prior flower annotation. 885 plant photographic observation records scanned (570 + 193 + 74 + 48). Only a subset met location accuracy <=100m, unobscured, wild/not explicitly cultivated, and strict **CC0 or CC-BY photo license**. Exactly **13** thumbnails could be selected across the 8 intended city × taxon × near/far strata, with 3 per available stratum, because some strata were empty:
- Osaka: *T. repens* rail <=250m: 3, and none suitable >1000m;
- Tokyo: *T. repens* near: 3, far: 1;
- Tokyo: *S. altissima* near: 3, far: 3;
- Osaka: *S. altissima* no suitable licensed near/far in this first sampling setup.

Only **3/13** of the selected photos carried an explicit iNaturalist Flowers annotation; **0/13** carried an explicit No Flowers or Fruits annotation. This is not an ecological estimate of flowering prevalence.

## Provisional manual photo-geometry screen

[Observation-level provisional screen CSV](photo_2026-11-07_provisional_review.csv), keyed solely to artifact `photo_001.jpg` through `photo_013.jpg`:

| Visibility category | Number of images |
|---|---:|
| Multi-individual/patch, no single whole plant separable | **7** |
| Partial shoots / missing base | **3** |
| Flower close-up | **1** |
| Hand-manipulated branch | **1** |
| Candidate full individual, but still visually ambiguous against other stems | **1** |

- **0/13** unambiguous, complete, isolated individual-plant silhouettes.
- **1/13** *tentative* full-plant candidate, not yet eligible for individual ROI.
- **11/13** showed probable visible open flowers on this review and 2 were classified uncertain, compared with only 3/13 explicitly annotated as Flowers on iNaturalist. These are provisional visual observations, not reliable site flowering probabilities and may reflect how photos were selected.
- **0/13** visibly showed an unmistakable active railway track/trackbed. This does NOT imply the coordinates are far from rails; it means **near-rail distance is not direct evidence of railway verge vegetation or slipstream/cutting exposure**.

These results are **positive for a narrow source-feasibility claim** (licensable photos and railway proximity GIS can be reproducibly assembled) and **negative for a broad morphological inference claim** (uncontrolled public photos do not yield repeatable absolute plant size or single-plant ROI shape without substantial screening).

## Decision: preserve two different ROI targets

1. **Individual plant ROI**: width/height, silhouette solidity, flower relative height; only after complete-plant view, isolated focal individual, same-angle and human-approved masks. No actual records pass this QC yet.
2. **Patch/scene ROI**: apparent floral pixels / *independently annotated conspecific vegetation pixels* could be studied within carefully selected **homogeneous, fixed-view** plant patches, with separate species identity QA and field angle controls. This is a visible floral display fraction, not absolute size, biomass or fitness. Existing photos may suggest this is more feasible but **no patch segmentation is yet validated**, so no ROI values are asserted.

`scripts/extract_photo_roi.py` now requires `target_scope` (`individual` vs `patch` / `closeup` / `unknown`) and `mask_reviewed`. It refuses to mark patches or unreviewed segmentations as **valid for individual-shape comparison**, while retaining pixel diagnostics.

## Stronger life-history tests beyond these records

- **Phenology:** fixed plot repeated surveys recording flowering *and* nonflowering, observer effort, GDD and true photo date, and verified mowing events. Look for secondary bloom peaks after mowing, not simply first photo date.
- **Train wind:** use actual measured gust exposure, trackside distance and shielding, not GIS rail distance alone.
- **Plasticity:** before-after of the same marked individual/plot under known cuts, with matched roadside plots.
- **Adaptation:** common-garden, seeds, or genetic comparisons after site histories are measured.
- **Macro generalization:** multiple independent corridors across regions (rail/road/background); photo counts within one park are not replicated railway systems.

Safety: no unauthorized entry to railway track or rights-of-way, and no unlicensed photo redistribution.
