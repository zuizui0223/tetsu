# v1.0 horsetail spore-cone ROI inputs (not flower ROI)

**Horsetail (*Equisetum arvense*) is a spore plant, not a flowering plant.** Fertile cones (ツクシの胞子嚢穂) and sterile/vegetative shoots have different timing and visibility. Never apply the iNaturalist Flowers/Fruits or a "first flower" variable to horsetail. Source: https://www.qsr.mlit.go.jp/chikugo/siryo/02-kawa/01_chikugo/syokubutu/sugina.html.

## Separate extractor

`scripts/extract_horsetail_spore_roi.py` uses the same validated pixel-mask geometry as the angiosperm `extract_photo_roi.py`, but exposes **different biological names** in the CSV output:

- `spore_cone_mask_annotated`
- `spore_cone_area_px`
- `spore_cone_area_over_plant_area`
- `spore_cone_relative_vertical_position`
- `valid_for_fixed_plot_spore_cone_fraction`
- `not_a_flowering_plant=true`

No `flower_area_px` column is emitted by this script. It does **not** estimate spore production, individual genetics, physical plant height, reproductive fitness, or management response.

Required manifest CSV:
```csv
image_id,site_id,taxon,photo_date,observer_id,plant_mask,spore_cone_mask,mask_reviewed,target_scope,whole_plant,standardized_view,fixed_plot_view,same_taxon_patch,plot_id,spore_phenophase
s1,site01,Equisetum arvense,2026-04-15,obs1,plant.png,cones.png,true,patch,false,true,true,true,plot01,fertile_cone_visible
```

This example is an **input shape only**, not a real rail observation. Both masks must be aligned 8-bit PNGs (background 0, foreground any positive pixel). Fertile cone pixels must be within the target plant foreground. Verified absence requires a supplied, manually annotated **all-black spore_cone_mask** and `fertile_cone_absent_visible`. A blank cone-mask field is **unknown/not annotated**, never spore-cone absence. Accepted phases:
- `fertile_cone_visible`
- `fertile_cone_absent_visible`
- `fertile_cone_senescent`
- `sterile_vegetative_shoots`
- `unknown`

```bash
pip install Pillow
python scripts/extract_horsetail_spore_roi.py --manifest data_manifest.csv --root path/to/masks --output results/horsetail_spore_roi.csv
```

**Individual silhouette:** only if `whole_plant=true`, `standardized_view=true`, `target_scope=individual`, `mask_reviewed=true` and no border truncation or too-small target ROI.

**Repeated fixed quadrat:** requires `target_scope=patch`, `plot_id`, `same_taxon_patch=true`, `fixed_plot_view=true`, `standardized_view=true`, `mask_reviewed=true`, and an explicit cone mask. Then reports `pass_fixed_patch_spore_fraction`. This *is not* a photograph of one genetically identical plant. The plant silhouette occupancy of a fixed plot can be used for repeated cover but never for absolute biomass, height or genotypic plasticity.

**Safety:** Inspect/photograph plants only from safely accessible public places or with rail operator authorization; a projected N02 rail distance does not establish exposure to train turbulence, herbicide, cutting or ballast substrate.

## Strong ecological hypothesis but not yet a modern measured outcome

Abe et al. 2012 published JR East survey gives *E. arvense* **24 of 70 trackbed occurrence reports** vs **22 of 642 slope reports** (relative species-mention share ~10-fold) [source reanalysis](japan_rail_microhabitat_reanalysis_v10.md). This suggests a possible *persistent rhizomatous shoot* filter, but historic maintenance, habitat definition and recorded species are observationally confounded.

To verify within-species plasticity/recovery, follow repeated *same plot* shoot-cover images before/after independently documented rail **cutting or herbicide/growth-retardant** work, with matched roadside and uncut control plots. Report cover and cone fraction separately, and use actual treatment dates. No field photos/masks or cutting logs are yet available in tetsu.
