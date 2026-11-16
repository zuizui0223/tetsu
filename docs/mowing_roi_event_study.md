# Fixed-plot grass cutting × photograph ROI pilot

Status 2026-11-07: **Input-contract and synthetic code tests only. No real cutting/phenology or image-mask response estimated.** This is intended to make useful field photo data immediately analyzable without pretending that uncontrolled iNaturalist closeups give plant height.

## What the paired field photos measure

Two independent observable quantities per **same fixed camera plot**:
- Vegetation projected occupancy: `plant_mask_pixels / fixed_frame_pixels`; a 2D **cover** measure for the same fixed view, NOT height, physical biomass or genotype.
- Floral display fraction: `flower_mask_pixels / reviewed conspecific vegetation mask pixels`; displayed flower pixels, NOT seeds, reproductive investment or phenological onset.

At each independently verified cutting date `t=0`, define:
- Baseline window: two or more complete reviewed fixed-plot photos from **t−28 through t−1 days**.
- Recovery window: two or more equivalent photos from **t+7 through t+42 days**.
- Compare medians of these two indicators; report absolute difference in projected fraction. **Do not divide by small or zero pre-cut fractions or invent missing dates.**
- If another documented cut falls in the window, mark the comparison restricted; don't pretend it is a single cutting event.
- No formal p values or treatment-effect claims, because seasonal growth and other management can also explain pre/post change.

Real causal strengthening requires paired **uncut controls** or a controlled experimental cutting schedule, with a similar camera rig, light and aspect, matched plant taxon, temperature/GDD, soil/urban context and site-year random effects. An operator's *general* grass-cutting policy cannot be substituted for a verified event date.

## Data contract and reproducibility

1. Annotate full-frame conspecific plant ROI and all open flowers in aligned PNG masks; manually review masks. A common plant patch needs `target_scope=patch`, `mask_reviewed=true`, `standardized_view=true`, `fixed_plot_view=true`, `same_taxon_patch=true` and a stable `plot_id`. Patch photos may touch the image border as they occupy a fixed frame.
2. Run `scripts/extract_photo_roi.py` to obtain `quality_status=patch_pass` and `valid_for_fixed_plot_flower_fraction=true`. This does NOT make patches valid for individual-plant silhouette analyses.
3. Supply independent management log `site_id,event_date,event_type,evidence` for real cutting events, not dates estimated from plant photographs.
4. Run `scripts/summarize_mowing_photo_roi.py --roi results/photo_roi.csv --events verified_cut_events.csv --output-dir results/mowing` and inspect `plot_event_roi.csv` and `audit.json`.

Pre-cut / recovery windows and both-photo minimum thresholds are chosen before outcome examination. Future windows require a separate named version and sensitivity comparison.

## Ecological interpretation

A positive cover recovery or rebloom signature may indicate **vegetative regrowth or plastic seasonal response** if repeated identical individual/plot measurements and the cutting event are confirmed. It does not by itself demonstrate that *railway-specific selection* or *heritable adaptation* caused the response. Selection requires comparing survival/reproductive performance of taxon/genotypes, common garden or reciprocal transplantation and separately measured railway conditions.

For wind, use calibrated train-induced peak gusts, shelter and local background wind measurements; never turn 'near passenger railway' into 'under strong train wind'.

## Field safety and licensing

Fixed photos must be acquired in publicly accessible locations with confirmed safe distance/permissions, or through rail operator approval. Never enter tracks/operating railway land without authorization. Record geometry, photo permissions and consistent on-site placement; preserve observer and data-source metadata.
