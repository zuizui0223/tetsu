#!/usr/bin/env python3
"""Public-license whole-photo ROI feasibility sampler (actual photos, no masks).

Sample CC0 or CC-BY photos of focal urban taxa, independently of their
flower/fruit annotation, across two geographic metros and near/far official
rail geometry. Produce a contact sheet and attribution manifest for review.

No claims about whole-plant silhouette, stress, mowing, phenology shifts
or genetic adaptation. No raw locations exported.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import time
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit

import requests
from PIL import Image, ImageDraw, ImageOps
from audit_inat_rail_phenophases import (AREAS, URL, N02, RateLimitedClient,
                                          acceptable_observation, rail_index_for_bbox)

TAXA = ("Trifolium repens", "Solidago altissima")
AREAS_PILOT = {"tokyo_bbox": AREAS["tokyo_bbox"], "osaka_bbox": AREAS["osaka_bbox"]}
ALLOWED = {"cc0", "cc-by"}


def allowed_photo_url(url):
    p = urlsplit(str(url or ""))
    host = (p.hostname or "").lower()
    return p.scheme == "https" and (
        host == "static.inaturalist.org"
        or host.endswith(".inaturalist.org")
        or host.startswith("inaturalist-open-data.s3")
    )


def licensed_photo(ob):
    for ph in ob.get("photos") or []:
        license_code = str(ph.get("license_code") or "").lower()
        url = ph.get("url") or ""
        if license_code in ALLOWED and allowed_photo_url(url):
            return ph
    return None


def hash_rank(obs_id):
    return hashlib.sha256(f"tetsu-v04-photo-{obs_id}".encode()).hexdigest()


def image_fetch_url(photo):
    url = photo["url"]
    # iNaturalist API commonly returns the square variant; medium has more
    # visual structure. Use host-preserving substitution only.
    for size in ("square", "thumb", "small"):
        if f"/{size}." in url:
            return url.replace(f"/{size}.", "/medium.")
    return url


def query_all(client, taxon, bbox):
    swlat, swlng, nelat, nelng = bbox
    fixed = {
        "taxon_name": taxon, "quality_grade": "research",
        "photos": "true", "d1": "2018-01-01", "d2": "2025-12-31",
        "swlat": swlat, "swlng": swlng, "nelat": nelat, "nelng": nelng,
        "per_page": 200, "order": "asc", "order_by": "id"
    }
    first = client.get(URL, params={**fixed, "page": 1}).json()
    total = int(first["total_results"])
    if total > 10000:
        raise ValueError(f"Data exceed iNaturalist's 10k search window for {taxon}: use official export.")
    all_results = list(first["results"])
    for page in range(2, math.ceil(total / 200) + 1):
        data = client.get(URL, params={**fixed, "page": page}).json()
        all_results.extend(data["results"])
    if len(all_results) != total:
        raise ValueError(f"iNaturalist total changed while paging: {len(all_results)}/{total}")
    ids = [i["id"] for i in all_results]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate observations across pages")
    return all_results


def group_for_rail(d):
    if d <= 250:
        return "rail_0_250m"
    if d >= 1000:
        return "rail_gt_1000m"
    return "middle_excluded"


def run(out, n_per_group, delay):
    client = RateLimitedClient(delay_seconds=delay)
    archive = client.get(N02, timeout=(15, 120)).content
    if not zipfile.is_zipfile(io.BytesIO(archive)):
        raise ValueError("N02 not ZIP")
    rail_maps = {
        city: rail_index_for_bbox(archive, bounds)[0]
        for city, bounds in AREAS_PILOT.items()
    }
    eligible = defaultdict(list)
    qc = Counter()
    total_by_stratum = {}
    for city, bbox in AREAS_PILOT.items():
        for taxon in TAXA:
            observations = query_all(client, taxon, bbox)
            total_by_stratum[city + "|" + taxon] = len(observations)
            for ob in observations:
                okay, reason = acceptable_observation(ob)
                if not okay:
                    qc[reason] += 1
                    continue
                photo = licensed_photo(ob)
                if photo is None:
                    qc["no_usable_CC0_or_CC_BY_photo"] += 1
                    continue
                lng, lat = ob["geojson"]["coordinates"]
                dist = rail_maps[city](lng, lat)
                group = group_for_rail(dist)
                if group == "middle_excluded":
                    qc[group] += 1
                    continue
                label_codes = sorted({
                    ann.get("controlled_value_id") for ann in (ob.get("annotations") or [])
                    if ann.get("controlled_attribute_id") == 12 and isinstance(ann.get("controlled_value_id"), int)
                })
                if 21 in label_codes and 13 in label_codes:
                    qc["conflicting_annotation_codes"] += 1
                eligible[(city, taxon, group)].append({
                    "obs_id": ob["id"], "obs_url": f"https://www.inaturalist.org/observations/{ob['id']}",
                    "city": city, "taxon": taxon, "rail_proximity": group,
                    "date": ob.get("observed_on"), "photographer_login": (ob.get("user") or {}).get("login"),
                    "photo": photo, "annotation_values": label_codes
                })
    manifest = []
    for (city, taxon, group), items in sorted(eligible.items()):
        items.sort(key=lambda x: hash_rank(x["obs_id"]))
        selected = items[:n_per_group]
        for record in selected:
            photo = record["photo"]
            image_url = image_fetch_url(photo)
            if not allowed_photo_url(image_url):
                qc["image_url_not_allowed"] += 1
                continue
            try:
                request = client.get(image_url, timeout=(15, 45))
                if len(request.content) > 5_000_000:
                    raise ValueError("image too large")
                image = Image.open(io.BytesIO(request.content))
                image = ImageOps.exif_transpose(image).convert("RGB")
                if min(image.size) < 120:
                    raise ValueError("image too small for structure review")
            except (OSError, ValueError, requests.RequestException) as exc:
                qc["photo_download_or_decode_failed"] += 1
                continue
            number = len(manifest) + 1
            filename = f"photo_{number:03d}.jpg"
            image.thumbnail((700, 700), Image.Resampling.LANCZOS)
            image.save(out / filename, quality=85)
            manifest.append({
                "photo_no": number, "file": filename, "city": city, "taxon": taxon,
                "rail_proximity": group, "obs_date": record["date"],
                "annotations_controlled_value_ids": ",".join(map(str, record["annotation_values"])),
                "attribution": photo.get("attribution") or "",
                "photographer_login": record["photographer_login"] or "",
                "license": photo.get("license_code"),
                "photo_source": image_url, "obs_url": record["obs_url"]
            })
    if not manifest:
        raise RuntimeError("No licensed photos passed quality gates")
    # Compact contact sheet for manual full-plant/framing QC.
    cols, tile_w, tile_h = 4, 280, 315
    nrows = math.ceil(len(manifest) / cols)
    sheet = Image.new("RGB", (tile_w * cols, tile_h * nrows), "white")
    draw = ImageDraw.Draw(sheet)
    for i, row in enumerate(manifest):
        with Image.open(out / row["file"]) as photo:
            im = photo.copy()
        im.thumbnail((tile_w - 12, 258), Image.Resampling.LANCZOS)
        x0, y0 = (i % cols) * tile_w, (i // cols) * tile_h
        sheet.paste(im, (x0 + (tile_w - im.width)//2, y0 + (258 - im.height)//2))
        draw.text((x0 + 7, y0 + 264), f"{row['photo_no']:03d} {row['taxon'][:18]}", fill="black")
        draw.text((x0 + 7, y0 + 283), f"{row['city'][:5]} {row['rail_proximity']}", fill="black")
        draw.text((x0 + 7, y0 + 300), f"{row['license']} | image for review ONLY", fill="black")
    sheet.save(out / "contact_sheet.jpg", quality=86)
    import csv
    with (out / "photo_attribution.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
    counts = Counter(f"{x['city']}|{x['taxon']}|{x['rail_proximity']}" for x in manifest)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": URL, "official_rail_source": N02,
        "intention": "photos for independent full-plant / flower ROI QC, NOT rail biological effects",
        "obs_scanned_by_city_and_taxon": total_by_stratum,
        "n_selected_licensed_photo_thumbnails": len(manifest),
        "selected_per_stratum": dict(counts),
        "qc_losses": dict(qc),
        "n_with_explicit_flowers_annotation": sum("13" in x["annotations_controlled_value_ids"].split(",") for x in manifest),
        "n_with_explicit_no_flower_annotation": sum("21" in x["annotations_controlled_value_ids"].split(",") for x in manifest),
        "notes": [
            "Sampled deterministically from all species photos independent of flowering annotation; selection is still conditional on photo and license.",
            "Images are licensed CC0 or CC BY; see photo_attribution.csv and original iNaturalist links.",
            "Nearest passenger rail at <=250m does not make this a trackside vegetation image.",
            "No foreground masks have been annotated or validated. Review each photo for whole plant, target identity, clipping and occlusion before measuring shape.",
            "No exact observation coordinates exported. Photos themselves may imply location/context.",
            "Do not infer flowering onset, response to cutting, wind stress or genetic adaptation."
        ]
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output-dir", default="results/inat/roi_photos")
    p.add_argument("--per-group", type=int, default=3)
    p.add_argument("--delay", type=float, default=1.10)
    a=p.parse_args()
    if not 1 <= a.per_group <= 20 or a.delay < 1:
        p.error("per-group must be 1..20; delay >=1.0s")
    out=Path(a.output_dir)
    out.mkdir(parents=True,exist_ok=True)
    run(out,a.per_group,a.delay)


if __name__=="__main__":
    main()
