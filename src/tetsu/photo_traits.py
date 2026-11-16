"""Explainable two-mask photographic traits for railway vegetation studies.

Input: human-verified/externally segmented full-plant and flower PNG masks.
Does NOT perform segmentation, estimate true physical size, diagnose wind
adaptation or distinguish genetically evolved traits from plasticity.
"""
from __future__ import annotations

from typing import Optional
from PIL import Image


def polygon_area(points):
    if len(points) < 3:
        return 0.0
    return abs(sum(
        a[0] * b[1] - a[1] * b[0]
        for a, b in zip(points, points[1:] + points[:1])
    )) / 2.0


def convex_hull(points):
    """Monotone-chain convex hull in pixel-corner coordinates."""
    pts = sorted(set(points))
    if len(pts) < 3:
        return pts

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def summarize_masks(
    plant_mask: Image.Image,
    flower_mask: Optional[Image.Image] = None,
    *,
    whole_plant: bool = False,
    standardized_view: bool = False,
    min_plant_fraction: float = 0.01,
) -> dict:
    """Measure geometry while retaining defensible quality gates.

    A completely black *supplied* flower mask means annotated no visible flower;
    an absent flower mask means unknown/not annotated and returns None.
    A mask image may be binary 0/255 or any nonzero pixel value.
    Solidity uses pixel-square corner geometry, so full rectangular masks have
    solidity exactly 1.0 (not an artifact > 1).
    """
    if not 0 < min_plant_fraction < 1:
        raise ValueError("min_plant_fraction must be in (0,1)")
    w, h = plant_mask.size
    if w < 3 or h < 3:
        raise ValueError("mask image too small")
    if flower_mask is not None and flower_mask.size != (w, h):
        raise ValueError("plant and flower masks differ in pixel dimensions")
    plant_pixels = list(plant_mask.convert("L").getdata())
    flower_pixels = (
        list(flower_mask.convert("L").getdata())
        if flower_mask is not None else None
    )

    rows = {}
    total_plant = 0
    xmin = w
    xmax = -1
    ymin = h
    ymax = -1
    edge_touches = False
    for i, pixel in enumerate(plant_pixels):
        if pixel == 0:
            continue
        x, y = i % w, i // w
        total_plant += 1
        xmin, xmax = min(xmin, x), max(xmax, x)
        ymin, ymax = min(ymin, y), max(ymax, y)
        if x in (0, w - 1) or y in (0, h - 1):
            edge_touches = True
        if y in rows:
            rows[y][0] = min(rows[y][0], x)
            rows[y][1] = max(rows[y][1], x)
        else:
            rows[y] = [x, x]
    if not total_plant:
        raise ValueError("empty plant mask")

    # Corner expansion makes hull area valid for rasterized foreground pixels.
    corners = []
    for y, (lo, hi) in rows.items():
        corners.extend([
            (lo, y), (lo, y + 1), (hi + 1, y), (hi + 1, y + 1)
        ])
    hull_area = polygon_area(convex_hull(corners))
    solidity = total_plant / hull_area if hull_area > 0 else None
    area_fraction = total_plant / (w * h)
    bbox_w = xmax - xmin + 1
    bbox_h = ymax - ymin + 1

    reasons = []
    if area_fraction < min_plant_fraction:
        reasons.append("target_roi_too_small")
    if edge_touches:
        reasons.append("plant_touches_image_border")
    if not whole_plant:
        reasons.append("not_whole_plant")
    if not standardized_view:
        reasons.append("view_not_standardized")

    result = {
        "frame_width_px": w, "frame_height_px": h,
        "plant_area_px": total_plant,
        "plant_frame_fraction": round(area_fraction, 6),
        "plant_bbox_width_px": bbox_w,
        "plant_bbox_height_px": bbox_h,
        "plant_width_height_ratio": round(bbox_w / bbox_h, 6),
        "plant_silhouette_solidity": round(solidity, 6) if solidity is not None else None,
        "plant_touches_image_border": edge_touches,
        "flower_mask_annotated": flower_pixels is not None,
        "flower_area_px": None,
        "flower_area_over_plant_area": None,
        "flower_relative_vertical_position": None,
        "valid_for_cross_image_shape_comparison": len(reasons) == 0,
        "quality_reasons": ";".join(reasons) if reasons else "pass",
        "note": "All features are pixel geometry, NOT true physical plant dimensions or biomass."
    }

    if flower_pixels is not None:
        flower_count = 0
        outside_count = 0
        ysum = 0
        for i, value in enumerate(flower_pixels):
            if value > 0:
                flower_count += 1
                ysum += (i // w) + 0.5
                outside_count += int(plant_pixels[i] == 0)
        if outside_count:
            raise ValueError(
                f"{outside_count} annotated floral pixels lie outside plant ROI"
            )
        result["flower_area_px"] = flower_count
        result["flower_area_over_plant_area"] = round(
            flower_count / total_plant, 6
        )
        if flower_count:
            # 0 ~ bottom, 1 ~ top of plant projection.
            result["flower_relative_vertical_position"] = round(
                (ymax + 1 - ysum / flower_count) / bbox_h, 6
            )
    return result
