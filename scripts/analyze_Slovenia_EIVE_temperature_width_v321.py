#!/usr/bin/env python3
"""v321 original-Slovenia vegetation joined to independent EIVE 1.0 thermal traits.

Phenomenon tested: within active railway shoulder (BS) versus railway toe
(BT) species temperature-indicator centre and niche-width proxy. Both zones
are railway. This does NOT identify climate offset, microclimate causality,
population adaptation, train dispersal or independently replicated climates.
"""
import argparse
import hashlib
import json
import math
import re
import unicodedata
import xml.etree.ElementTree as ET
import zipfile
from collections import defaultdict
from pathlib import Path

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
from analyze_Slovenia_thermal_generalist_section_v319 import (
    make_pairs, section_test, section_drop)
from audit_EIVE_temperature_niche_width_source_v320 import (
    EXPECTED_MD5, S, cell_value, colnum)

EIVE_COLUMNS=("TaxonConcept","TaxonRank","EIVEres-T",
              "EIVEres-T.nw3","EIVEres-T.n")
SPECIES_RE=re.compile(r"^[A-Z][a-zA-Z-]+ [a-z][a-zA-Z-]+$")
NON_SPECIES_EPITHETS={"sect","section","agg","aggregate","sp","spp",
    "subsp","var","cf","aff","group","complex","spec","species"}

def normal_name(v):
    return " ".join(unicodedata.normalize("NFKC",str(v)).split()).casefold()

def latin_map(original):
    codebook=original.get("species_codebook",[])
    if len(codebook)!=312: raise ValueError("Unexpected original codebook source size")
    names={}
    for row in codebook:
        code=str(row["EPPO_code"])
        matches=set()
        for key,v in row.items():
            if not any(k in str(key).lower() for k in
                       ("scientific","latin","species","taxon","name")):
                continue
            if key=="EPPO_code" or not isinstance(v,str):continue
            m=re.match(r"^([A-Z][a-zA-Z-]+ [a-z][a-zA-Z-]+)(?:\b|$)",v.strip())
            if m and SPECIES_RE.fullmatch(m.group(1)):
                if m.group(1).split()[1].lower() not in NON_SPECIES_EPITHETS:
                    matches.add(m.group(1))
        if len(matches)==1:
            names[code]=normal_name(next(iter(matches)))
    if len(names)<100:
        raise ValueError("Too few explicit original Latin binomials")
    return names

def float_finite(v):
    try:
        x=float(v)
    except (ValueError,TypeError):
        return None
    return x if math.isfinite(x) else None

def eive_table(path):
    data=Path(path).read_bytes()
    if hashlib.md5(data).hexdigest()!=EXPECTED_MD5:
        raise ValueError("EIVE source MD5 mismatch")
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:raise ValueError("EIVE ZIP CRC fail")
        shared=ET.fromstring(z.read("xl/sharedStrings.xml"))
        strings=["".join(t.text or "" for t in item.iter(S+"t"))
                 for item in shared.findall(S+"si")]
        sheet=ET.fromstring(z.read("xl/worksheets/sheet2.xml"))
        rows=sheet.findall(".//"+S+"sheetData/"+S+"row")
        if len(rows)!=14836:
            raise ValueError("EIVE published mainTable 14836 rows expected")
        def unpack(row):
            return {colnum(c.attrib["r"]):cell_value(c,strings)
                    for c in row.findall(S+"c")}
        header=unpack(rows[0])
        mapping={str(v):k for k,v in header.items()}
        if any(c not in mapping for c in EIVE_COLUMNS):
            raise ValueError("Missing EIVE temperature position/width/n columns")
        out={}
        duplicates=0
        species_rows=0
        for row in rows[1:]:
            cells=unpack(row)
            if str(cells.get(mapping["TaxonRank"],"")).casefold()!="species":
                continue
            species_rows+=1
            raw=cells.get(mapping["TaxonConcept"],"")
            name=normal_name(raw)
            # Do not coerce hybrid names, infra-specific names or authors.
            if not re.fullmatch(r"[a-z-]+ [a-z-]+",name):
                continue
            centre=float_finite(cells.get(mapping["EIVEres-T"]))
            breadth=float_finite(cells.get(mapping["EIVEres-T.nw3"]))
            systems=float_finite(cells.get(mapping["EIVEres-T.n"]))
            if centre is None or breadth is None or systems is None:
                continue
            if breadth<=0 or centre<=0 or systems<1 or not systems.is_integer():
                continue
            item={"position":centre,"width":breadth,"systems":int(systems)}
            if name in out:
                duplicates+=1
                if out[name]!=item:
                    raise ValueError("Duplicate EIVE accepted names with conflicting traits")
            out[name]=item
    return out,{"EIVE_total_mainTable_data_rows":len(rows)-1,
                "EIVE_taxon_rank_species_rows":species_rows,
                "EIVE_strict_binomials_with_temperature_position_and_width":len(out),
                "duplicate_same_name_same_value_rows":duplicates,
                "EIVE_original_md5":EXPECTED_MD5}

def score(codes, crosswalk, eive, metric, min_systems):
    vals=[]
    for c in codes:
        name=crosswalk.get(c)
        tr=eive.get(name) if name else None
        if tr is None or tr["systems"]<min_systems:continue
        vals.append(tr[metric])
    return (sum(vals)/len(vals) if vals else None,
            len(vals)/len(codes) if codes else 0)

def compare(pairs,names,eive,metric,min_systems,coverage):
    bysite=defaultdict(list)
    observed_cover={"BS":[],"BT":[]}
    for p in pairs:
        a,ca=score(p["BS"],names,eive,metric,min_systems)
        b,cb=score(p["BT"],names,eive,metric,min_systems)
        observed_cover["BS"].append(ca);observed_cover["BT"].append(cb)
        if a is None or b is None or min(ca,cb)<coverage:
            continue
        bysite[(p["section"],p["site"])].append(a-b)
    sections=defaultdict(list)
    for (sec,site),diffs in bysite.items():
        sections[sec].append(sum(diffs)/len(diffs))
    location_diffs=[x for v in sections.values() for x in v]
    return {"metric":metric,"min_source_EIV_systems":min_systems,
            "minimum_coverage_each_BS_BT_zone":coverage,
            "numeric_trait_coverage_across_all_paired_sides":{
                z:sum(v)/len(v) for z,v in observed_cover.items()},
            "eligible_original_track_sides":sum(len(x) for x in bysite.values()),
            "eligible_distinct_field_locations":len(bysite),
            "eligible_line_sections":len(sections),
            "field_location_mean_BS_minus_BT":(
                sum(location_diffs)/len(location_diffs) if location_diffs else None),
            "BS_greater_sites":sum(x>0 for x in location_diffs),
            "BT_greater_sites":sum(x<0 for x in location_diffs),
            "exact_section_signed_direction_diagnostic":section_test(sections),
            "leave_one_section_out":section_drop(sections)}

def analyze(original, path):
    pairs,traits,_,_ = make_pairs(original)
    names=latin_map(original)
    eive,source=eive_table(path)
    observed=set().union(*(p["BS"]|p["BT"] for p in pairs))
    matched={c for c in observed if names.get(c) in eive}
    strict3={c for c in matched if eive[names[c]]["systems"]>=3}
    if len(matched)<30:
        raise ValueError("Insufficient exact original Latin binomial EIVE trait match")
    panels=[]
    for metric in ("position","width"):
        for minsystems in (1,2,3):
            for coverage in (0,0.5):
                panels.append(compare(pairs,names,eive,metric,minsystems,coverage))
    return {"version":"v321","independent_eive_source_doi":"10.5281/zenodo.7534792",
      "slovenia_source_doi":"10.1007/s44473-026-00155-z",
      "scope":"same-site BS vs BT within railway; NOT regional climate offset",
      "source":source,"source_original_name_codebook_mapped":len(names),
      "n_rail_EPPO_codes_observed_in_original_BS_or_BT":len(observed),
      "n_exact_mapped_BS_BT_taxa_with_EIVE_temperature_traits":len(matched),
      "n_exact_mapped_BS_BT_taxa_with_ge3_EIVE_systems":len(strict3),
      "taxon_crosswalk":"exact binomial from original publisher codebook to EIVE TaxonConcept; no fuzzy or inferred synonym",
      "panels":panels,
      "limits":[
        "EIVE consensus width is an indicator-model niche-width proxy, not a measured thermal tolerance breadth",
        "n=1 system estimates may have model-imposed baseline width; evaluate >=2 and >=3 separately",
        "Original source nomenclature mismatches are excluded rather than guessed",
        "Two zones are railway; no same-locality natural vegetation control",
        "P-values are exploratory section sign reversal conditional on assumed section-direction exchangeability; not independent climate regions",
        "BS and BT original plot areas and measured temperature exposures not proven equal"]}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-slovenia-xlsx",required=True)
    ap.add_argument("--verified-v102",required=True)
    ap.add_argument("--eive-xlsx",required=True)
    ap.add_argument("--output",required=True)
    args=ap.parse_args()
    original=extract_original_xlsx(args.source_slovenia_xlsx,args.verified_v102)
    out=analyze(original,args.eive_xlsx)
    p=Path(args.output)
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({k:out[k] for k in (
        "version","source","source_original_name_codebook_mapped",
        "n_rail_EPPO_codes_observed_in_original_BS_or_BT",
        "n_exact_mapped_BS_BT_taxa_with_EIVE_temperature_traits",
        "n_exact_mapped_BS_BT_taxa_with_ge3_EIVE_systems","panels")},ensure_ascii=False),flush=True)

if __name__=="__main__":main()
