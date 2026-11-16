#!/usr/bin/env python3
"""Extract plant and floral pixel-ROI descriptors from reviewed PNG masks.

Requires: Pillow. See docs/rail_environment_traits_protocol.md.
No ML segmentation, plant size in cm, phenology timing or train effects inferred.
"""
import argparse
import csv
import sys
from pathlib import Path
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.photo_traits import summarize_masks

REQ = {
    "image_id", "site_id", "taxon", "photo_date", "photographer_id",
    "site_class", "whole_plant", "standardized_view", "phenophase",
    "plant_mask", "flower_mask", "target_scope", "mask_reviewed",
}
ALLOWED_PHASES = {
    "unknown", "vegetative", "bud", "flowering", "fruiting",
    "flowering_fruiting",
}
ALLOWED_SITES = {"rail", "road", "background", "unknown"}
ALLOWED_SCOPES = {"individual", "patch", "closeup", "unknown"}
BOOLEANS = {"true": True, "false": False, "1": True, "0": False}


def convert_bool(value):
    key = str(value).strip().lower()
    if key not in BOOLEANS:
        raise ValueError(f"expected true/false/1/0, got {value!r}")
    return BOOLEANS[key]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True, help="CSV with photo and mask metadata")
    p.add_argument("--root", required=True, help="directory containing masks")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    root = Path(args.root).resolve()
    rows = []
    with open(args.manifest, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        missing = REQ.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"manifest missing: {sorted(missing)}")
        for row in reader:
            meta = {k: (row[k] or "").strip() for k in REQ}
            # Optional additional fixed-quadrat metadata for patch-only floral
            # displays. Never infer these flags from GPS/photo composition.
            for key in ("plot_id", "same_taxon_patch", "fixed_plot_view"):
                meta[key] = (row.get(key) or "").strip()
            record = dict(meta)
            record["quality_status"] = "rejected"
            record["quality_error"] = ""
            try:
                if not meta["image_id"] or not meta["taxon"] or not meta["site_id"]:
                    raise ValueError("missing image_id, site_id or taxon")
                if meta["site_class"] not in ALLOWED_SITES:
                    raise ValueError("invalid site_class")
                if meta["phenophase"] not in ALLOWED_PHASES:
                    raise ValueError("invalid phenophase")
                if meta["taxon"].startswith("Equisetum "):
                    # Horsetails make spores in fertile strobili, not
                    # flowers/fruits. Never allow the flower ROI output or
                    # flowering phenophase to silently mislabel horsetail.
                    if meta["flower_mask"] or meta["phenophase"] in (
                        "flowering", "flowering_fruiting", "bud", "fruiting"
                    ):
                        raise ValueError(
                            "Equisetum has no flowers: use "
                            "extract_horsetail_spore_roi.py for fertile cones"
                        )
                if meta["target_scope"] not in ALLOWED_SCOPES:
                    raise ValueError("target_scope must be individual/patch/closeup/unknown")
                reviewed = convert_bool(meta["mask_reviewed"])
                if not meta["plant_mask"]:
                    raise ValueError("missing plant mask")
                def resolve_mask(path):
                    if not path:
                        return None
                    target = (root / path).resolve()
                    if not target.is_relative_to(root):
                        raise ValueError("mask path escapes --root")
                    if target.suffix.lower() != ".png":
                        raise ValueError("masks must be PNG")
                    return target

                plant_path = resolve_mask(meta["plant_mask"])
                flower_path = resolve_mask(meta["flower_mask"])
                with Image.open(plant_path) as plant:
                    if flower_path is None:
                        features = summarize_masks(
                            plant,
                            whole_plant=convert_bool(meta["whole_plant"]),
                            standardized_view=convert_bool(meta["standardized_view"]),
                        )
                    else:
                        with Image.open(flower_path) as flower:
                            features = summarize_masks(
                                plant, flower,
                                whole_plant=convert_bool(meta["whole_plant"]),
                                standardized_view=convert_bool(meta["standardized_view"]),
                            )
                record.update(features)
                if (features["flower_area_px"] or 0) > 0 and meta["phenophase"] in ("vegetative", "bud", "fruiting"):
                    raise ValueError("visible flower pixels contradict recorded phenophase")
                # A photograph of a patch or a close-up MUST NOT be treated as
                # an individual plant silhouette, even with a filled ROI.
                individual_and_reviewed = (
                    meta["target_scope"] == "individual" and reviewed
                )
                valid = (
                    features["valid_for_cross_image_shape_comparison"]
                    and individual_and_reviewed
                )
                reasons = (
                    [] if features["quality_reasons"] == "pass"
                    else features["quality_reasons"].split(";")
                )
                if meta["target_scope"] != "individual":
                    reasons.append("not_individual_target_scope")
                if not reviewed:
                    reasons.append("mask_not_independently_reviewed")
                record["valid_for_cross_image_shape_comparison"] = valid
                record["quality_reasons"] = ";".join(reasons) if reasons else "pass"
                fixed_view = (
                    convert_bool(meta["fixed_plot_view"]) if meta["fixed_plot_view"]
                    else False
                )
                same_taxon = (
                    convert_bool(meta["same_taxon_patch"]) if meta["same_taxon_patch"]
                    else False
                )
                patch_reasons = []
                if meta["target_scope"] != "patch":
                    patch_reasons.append("not_patch_scope")
                if not reviewed:
                    patch_reasons.append("mask_not_independently_reviewed")
                if not fixed_view or not convert_bool(meta["standardized_view"]):
                    patch_reasons.append("not_fixed_standardized_plot_view")
                if not same_taxon:
                    patch_reasons.append("target_patch_not_verified_conspecific")
                if not meta["plot_id"]:
                    patch_reasons.append("missing_plot_id")
                if not features["flower_mask_annotated"]:
                    patch_reasons.append("flower_mask_unannotated")
                if features["plant_frame_fraction"] < 0.01:
                    patch_reasons.append("plant_roi_too_small")
                # A patch can legitimately touch the image boundary if the
                # quadrat camera frame is fixed; this is NOT individual size.
                valid_patch = not patch_reasons
                record["valid_for_fixed_plot_flower_fraction"] = valid_patch
                record["patch_quality_reasons"] = (
                    ";".join(patch_reasons) if patch_reasons else "pass"
                )
                record["quality_status"] = (
                    "pass" if valid else ("patch_pass" if valid_patch else "restricted")
                )
            except (ValueError, OSError) as exc:
                record["quality_error"] = f"{type(exc).__name__}: {exc}"
            rows.append(record)

    if not rows:
        raise SystemExit("empty manifest")
    if len({x["image_id"] for x in rows}) != len(rows):
        raise SystemExit("duplicate image_id: can produce pseudoreplication")
    cols = list(rows[0])
    cols += [k for row in rows for k in row if k not in cols]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        writer.writerows(rows)
    print(f"images={len(rows)} pass={sum(x['quality_status']=='pass' for x in rows)} "
          f"patch_pass={sum(x['quality_status']=='patch_pass' for x in rows)} "
          f"restricted={sum(x['quality_status']=='restricted' for x in rows)} "
          f"rejected={sum(x['quality_status']=='rejected' for x in rows)}")


if __name__ == "__main__":
    main()
