#!/usr/bin/env python3
"""v339 C4 preflight: source-verified TAXON NAME INVENTORY ONLY.
No field plot, site, coordinate or source abundance values are exported.
"""
import argparse, hashlib, json
from pathlib import Path
from collections import Counter
from analyze_Kashmir_rail_interior_EIVE_thermal_v327 import species_names_from_publisher
from analyze_Kashmir_same_sites_reference_year_v254 import DRYAD_SHA256

def run(path):
    digest=hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if digest != DRYAD_SHA256:
        raise ValueError("Authentic original Dryad checksum mismatch")
    names=species_names_from_publisher(path)
    valid=sorted(set(n for n in names if n))
    return {
       "version":"v339-taxonomy-preflight",
       "original_source_sha256":digest,
       "original_source_species_headers":len(names),
       "parsed_binomials":len(valid),
       "genera":len({x.split()[0] for x in valid}),
       "species_names_only_no_field_records":valid,
       "genus_histogram":dict(sorted(Counter(x.split()[0] for x in valid).items()))
    }

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--kashmir-original",required=True)
    ap.add_argument("--output",required=True)
    args=ap.parse_args()
    out=run(args.kashmir_original)
    dest=Path(args.output);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False),flush=True)
