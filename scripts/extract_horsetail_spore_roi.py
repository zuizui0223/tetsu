#!/usr/bin/env python3
"""Measure visible spore-cone ROI in independently reviewed horsetail photos.

Horsetails (*Equisetum*) reproduce by spores; flowering/floral ROI and
iNaturalist Flowers/Fruits labels are biologically inapplicable.
This script reports projected fertile strobilus fraction only, not spore
yield, population reproduction, individual dry mass or true plant height.
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tetsu.photo_traits import summarize_masks

REQUIRED={
    "image_id", "site_id", "taxon", "photo_date", "observer_id",
    "plant_mask", "spore_cone_mask", "mask_reviewed",
    "target_scope", "whole_plant", "standardized_view",
    "fixed_plot_view", "same_taxon_patch", "plot_id",
    "spore_phenophase",
}
PHASES={"unknown", "fertile_cone_visible", "fertile_cone_absent_visible",
        "fertile_cone_senescent", "sterile_vegetative_shoots"}
TRUE={"true", "1"}
FALSE={"false", "0"}


def boolean(x):
    x=str(x).strip().lower()
    if x in TRUE: return True
    if x in FALSE: return False
    raise ValueError("required boolean must be true/false/1/0")


def resolve(root:Path, value:str):
    if not value: return None
    p=(root/value).resolve()
    if not p.is_relative_to(root) or p.suffix.lower()!=".png":
        raise ValueError("mask must be PNG inside --root")
    return p


def summarize_row(row:dict,root:Path)->dict:
    result=dict(row)
    result["quality_status"]="rejected"
    result["quality_error"]=""
    try:
        if not row["image_id"] or not row["site_id"]:
            raise ValueError("missing photo/site identifier")
        if not row["taxon"].startswith("Equisetum "):
            raise ValueError("spore-cone analysis only accepts Equisetum")
        if row["spore_phenophase"] not in PHASES:
            raise ValueError("invalid horsetail spore phase")
        reviewed=boolean(row["mask_reviewed"])
        whole=boolean(row["whole_plant"])
        standardized=boolean(row["standardized_view"])
        fixed=boolean(row["fixed_plot_view"])
        conspecific=boolean(row["same_taxon_patch"])
        target_scope=row["target_scope"]
        if target_scope not in {"individual", "patch", "unknown", "closeup"}:
            raise ValueError("invalid target_scope")
        plant_path=resolve(root,row["plant_mask"])
        if plant_path is None:
            raise ValueError("plant mask required")
        cone_path=resolve(root,row["spore_cone_mask"])
        with Image.open(plant_path) as plant:
            if cone_path is None:
                measures=summarize_masks(plant,whole_plant=whole,
                                          standardized_view=standardized)
            else:
                with Image.open(cone_path) as cones:
                    measures=summarize_masks(plant,cones,whole_plant=whole,
                                              standardized_view=standardized)
        if row["spore_phenophase"] in (
            "fertile_cone_visible", "fertile_cone_absent_visible"
        ) and not measures["flower_mask_annotated"]:
            raise ValueError("visible/absent cone phase requires independently annotated cone mask")
        if row["spore_phenophase"]=="fertile_cone_visible" and measures["flower_area_px"]==0:
            raise ValueError("cone-present label contradicts reviewed zero-pixel cone mask")
        if row["spore_phenophase"] in (
            "fertile_cone_absent_visible", "sterile_vegetative_shoots"
        ) and (measures["flower_area_px"] or 0)>0:
            raise ValueError("cone-absent phase contradicts visible cone mask")
        # Convert generic mask math into BIOLOGICALLY ACCURATE spore labels.
        measure_names={
            "flower_mask_annotated":"spore_cone_mask_annotated",
            "flower_area_px":"spore_cone_area_px",
            "flower_area_over_plant_area":"spore_cone_area_over_plant_area",
            "flower_relative_vertical_position":"spore_cone_relative_vertical_position",
        }
        for key, value in measures.items():
            result[measure_names.get(key,key)]=value
        # Never label projected patch as a single plant, nor claim
        # physiological stress without exposure.
        individual_valid=(
            measures["valid_for_cross_image_shape_comparison"]
            and reviewed and target_scope=="individual"
        )
        patch_valid=(
            reviewed and target_scope=="patch" and fixed and
            standardized and conspecific and bool(row["plot_id"]) and
            measures["flower_mask_annotated"] and
            measures["plant_frame_fraction"]>=0.01
        )
        result["valid_for_cross_image_shape_comparison"]=individual_valid
        result["valid_for_fixed_plot_spore_cone_fraction"]=patch_valid
        result["quality_status"]=(
            "pass_individual_shape" if individual_valid else
            "pass_fixed_patch_spore_fraction" if patch_valid else
            "restricted"
        )
        result["reproductive_structure"]="spore_strobilus"
        result["not_a_flowering_plant"]=True
        result["metrics_are_projected_pixels_not_spore_yield"]=True
    except (OSError, ValueError) as exc:
        result["quality_error"]=f"{type(exc).__name__}: {exc}"
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",required=True)
    p.add_argument("--root",required=True)
    p.add_argument("--output",required=True)
    args=p.parse_args()
    root=Path(args.root).resolve()
    with open(args.manifest,newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        missing=REQUIRED-set(reader.fieldnames or [])
        if missing: raise SystemExit(f"missing fields: {sorted(missing)}")
        rows=list(reader)
    if not rows:raise SystemExit("no photo metadata")
    if len({r["image_id"] for r in rows})!=len(rows):
        raise SystemExit("duplicate image ID")
    results=[summarize_row(r,root) for r in rows]
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    cols=list(dict.fromkeys(k for r in results for k in r))
    with out.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=cols)
        writer.writeheader();writer.writerows(results)
    print(f"images={len(results)} "
          f"valid_individuals={sum(r['quality_status']=='pass_individual_shape' for r in results)} "
          f"valid_fixed_patches={sum(r['quality_status']=='pass_fixed_patch_spore_fraction' for r in results)} "
          f"restricted={sum(r['quality_status']=='restricted' for r in results)} "
          f"rejected={sum(r['quality_status']=='rejected' for r in results)}")


if __name__=="__main__":
    main()
