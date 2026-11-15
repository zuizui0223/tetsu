#!/usr/bin/env python3
"""v18.7: extract Parma 2020 complete author 438 taxa × 28 station incidence.

Publisher original Adorni & Ghillani 2020, Annali Museo Civico di Rovereto
36:149–200. This parser reads PDF-native text and faithfully retains printed
station codes, a separate printed [n] count, biological lifeforms and alien
status. It does not equate whole stations to standardized ballast plots.
Survey only 28 stations; Citerna Taro (CIT) was NOT surveyed.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,re
from collections import Counter
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_Parma_2020_28rail_station_flora_original_v186 import (
    STATION_CODES, ORIGINAL_URL, ORIGINAL_BIBLIO, EXCLUDED_UNSURVEYED_STATION)
START="Abies alba Mill."
END_MARK="4.2 Diversità floristica"
CLOSING_RE=re.compile(r"\[(\d{1,2})\]")
STATION_RE=re.compile(r"\b(?:"+"|".join(re.escape(x) for x in STATION_CODES)+r")\b")
FAMILY_RE=re.compile(r"^\s*(.+?)\s+\(([A-Za-zÀ-ž-]+aceae)\),\s*(.*)$",re.S)
ALIEN_RE=re.compile(r"\b(Autoct\.|Archeof\.\s+\w+\.|Neof\.\s+\w+\.|Xenof\.\s+\w+\.)")
LIFE_RE=re.compile(r",\s*(T|H|P|G|Ch|He|NP|I)(?:\s+([A-Za-zÀ-ž-]+))?\s*,")
EXPECTED_TAXA=438
KNOWN_BENCHMARKS={"SAL":157,"BVT":151,"PAR":145,"CHI":32}
EXTRA_EFFORT_STATIONS=set("BUS PAR FID MEZ NOC".split())

def original_lines(source):
    import fitz
    raw=Path(source).read_bytes()
    if not raw.startswith(b"%PDF-") or len(raw)<300000:
        raise ValueError("Publisher PDF missing / malformed")
    doc=fitz.open(stream=raw,filetype="pdf")
    if len(doc)!=54:raise ValueError("Publisher original source PDF not 54 pages")
    text_lines=[];opened=False;closed=False;pages=[]
    for pageno,page in enumerate(doc):
        text=page.get_text("text",sort=True)
        if not opened:
            i=text.find(START)
            if i<0:continue
            opened=True
            text=text[i:]
        if "4.2 Diversità floristica" in text or re.search(r"\b4\.2\s+(?:Diversit|Diversità)",text):
            text=re.split(r"4\.2\s+(?:Diversit|Diversità)",text,maxsplit=1)[0]
            closed=True
        for line in text.splitlines():
            line=line.strip()
            if not line:continue
            if re.fullmatch(r"\d{3}",line):continue
            if "Annali Museo Civico di Rovereto" in line:continue
            text_lines.append(line)
            pages.append(pageno)
        if closed:break
    if not opened or not closed:raise ValueError("First taxon or 4.2 terminal section missing in source PDF")
    return text_lines,pages,dict(sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),pages=len(doc))

def concatenate(lines,pages):
    # Keep original page at the beginning of each record for provenance.
    cleaned=[]
    for line in lines:
        cleaned.append(line)
    joined="\n".join(cleaned)
    # Only original PDF word-wrapped hyphen immediately followed by lower-case
    # continuation is repaired. Do NOT strip taxonomic or chorological hyphens
    # where the next word begins upper-case.
    joined=re.sub(r"([A-Za-zÀ-ž])-\s*\n\s*([a-zà-ž])",r"\1\2",joined)
    joined=re.sub(r"\s+"," ",joined)
    return joined.strip()

def parse_entry(raw,printed_count):
    entry=re.sub(r"\s+"," ",raw).strip()
    family=FAMILY_RE.match(entry)
    if not family:
        return None,dict(cause="no_publisher_family_pattern",fragment=entry[:180])
    taxon_original=family.group(1).strip()
    after=family.group(3).strip()
    observed=sorted(set(STATION_RE.findall(after)))
    life_matches=list(LIFE_RE.finditer(","+after))
    life_full=None
    life_category=None
    if life_matches:
        m=life_matches[-1]
        life_full=m.group(1)+(" "+m.group(2) if m.group(2) else "")
        life_category=m.group(1)
        if life_category=="NP":life_category="P"  # nanophanerophytes
    status=ALIEN_RE.search(after)
    # 'Autoct.' may be absent from some source records; preserve missingness.
    if not re.match(r"^[A-Z][A-Za-zÀ-ž-]+\s+[a-z][A-Za-zÀ-ž-]+",taxon_original):
        return None,dict(cause="no_species_binomial",fragment=entry[:180])
    return {
        "taxon_original":taxon_original,
        "family_original":family.group(2),
        "native_status_original":status.group(1) if status else None,
        "lifeform_original":life_full,
        "lifeform_coarse":life_category,
        "station_codes":observed,
        "printed_station_total":printed_count,
        "parsed_station_total":len(observed),
        "printed_minus_codes":printed_count-len(observed),
        "original_line_text":entry,
    },None

def parse_pdf(pdf):
    lines,pages,source_meta=original_lines(pdf)
    joined=concatenate(lines,pages)
    starts=joined.startswith(START)
    if not starts:raise ValueError("Parsed section no longer begins with Abies alba")
    records=[];errors=[]
    pos=0
    for m in CLOSING_RE.finditer(joined):
        part=joined[pos:m.start()].strip()
        pos=m.end()
        if not part:continue
        row,error=parse_entry(part,int(m.group(1)))
        if error:errors.append(error)
        else:records.append(row)
    tail=joined[pos:].strip()
    return records,errors,{
        "source_original_url":ORIGINAL_URL,"source_bibliography":ORIGINAL_BIBLIO,
        "source_pdf_sha256":source_meta["sha256"],
        "source_pdf_bytes":source_meta["bytes"],
        "source_pdf_pages":source_meta["pages"],
        "n_original_source_rows_parsed":len(records),
        "n_publisher_reports_missing_from_parser":EXPECTED_TAXA-len(records),
        "n_parse_failures":len(errors),
        "parse_errors_preview":errors[:35],
        "section_trailing_unparsed_excerpt":tail[:250],
        "n_taxon_name_unique":len(set(x["taxon_original"] for x in records)),
        "native_status_missing_rows":sum(x["native_status_original"] is None for x in records),
        "lifeform_missing_rows":sum(x["lifeform_coarse"] is None for x in records),
        "publisher_frequency_disagreements":[
            {"taxon":x["taxon_original"],"printed_n":x["printed_station_total"],
             "reconstructed_station_codes_n":x["parsed_station_total"],
             "diff":x["printed_minus_codes"]}
            for x in records if x["printed_station_total"]!=x["parsed_station_total"]
        ][:70],
        "publisher_frequency_disagreements_count":sum(
            x["printed_station_total"]!=x["parsed_station_total"] for x in records),
        "diagnostic_raw_record_examples":records[:4]+records[-4:],
        "expect_publisher_taxa":EXPECTED_TAXA,
        "source_station_code_count":len(STATION_CODES),
        "CIT_not_surveyed_not_zero":True
    }

def station_stats(rows):
    richness={code:sum(code in x["station_codes"] for x in rows) for code in STATION_CODES}
    family_count=Counter(x["family_original"] for x in rows)
    life=Counter(x["lifeform_coarse"] for x in rows)
    natives=Counter(x["native_status_original"] for x in rows)
    return richness,family_count,life,natives

def save(rows,diag,outdir):
    root=Path(outdir);root.mkdir(parents=True,exist_ok=True)
    richness,families,life,natives=station_stats(rows)
    diag["original_reported_sample_station_richness"]=KNOWN_BENCHMARKS
    diag["observed_station_richness"]=richness
    diag["original_reported_station_richness_reproduced"]={key:richness.get(key)==value for key,value in KNOWN_BENCHMARKS.items()}
    diag["source_family_count"]=len(families)
    diag["lifeform_counts"]=dict(life)
    diag["status_counts"]=dict(natives)
    diag["all_28_stations_with_positive_occurrences"]=all(richness.values())
    diag["station_richness_min"]=min(richness.values()) if richness else None
    diag["station_richness_max"]=max(richness.values()) if richness else None
    diag["mean_station_richness"]=sum(richness.values())/28 if richness else None
    diag["passed_taxa_438"]=len(rows)==EXPECTED_TAXA and diag["n_taxon_name_unique"]==438
    diag["passed_field_stations_28"]=len(richness)==28 and diag["all_28_stations_with_positive_occurrences"]
    diag["passed_original_4_richness_benchmarks"]=all(diag["original_reported_station_richness_reproduced"].values())
    diag["source_complete_matrix_gate"]=bool(
        diag["passed_taxa_438"] and diag["passed_field_stations_28"]
        and diag["passed_original_4_richness_benchmarks"]
        and not diag["n_parse_failures"]
        and not diag["lifeform_missing_rows"]
    )
    audit_path=root/"Parma_2020_original_438taxa_28station_QA.json"
    audit_path.write_text(json.dumps(diag,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({k:v for k,v in diag.items()
                      if k not in ("diagnostic_raw_record_examples","parse_errors_preview",
                                   "publisher_frequency_disagreements")},ensure_ascii=False,indent=2),flush=True)
    if not diag["source_complete_matrix_gate"]:
        raise ValueError("Publisher original 438×28 original sources do not pass identity/richness gates; do not claim a genuine incidence matrix")
    with (root/"Parma_2020_28stations_438taxa_original_incidence.csv").open("w",newline="",encoding="utf8") as out:
        w=csv.writer(out)
        w.writerow(["original_study","country","station_code","taxon_original","presence",
                    "lifeform_original","lifeform_coarse","status_original","family_original",
                    "extra_observation_effort_known","whole_station_mixed_habitat"])
        for r in rows:
            for code in STATION_CODES:
                w.writerow(["Adorni_Ghillani_2020","IT",code,r["taxon_original"],
                            int(code in r["station_codes"]),r["lifeform_original"],
                            r["lifeform_coarse"],r["native_status_original"],r["family_original"],
                            int(code in EXTRA_EFFORT_STATIONS),True])
    with (root/"Parma_2020_438taxa_original_metadata.csv").open("w",newline="",encoding="utf8") as out:
        w=csv.writer(out)
        w.writerow(["taxon_original","family_original","lifeform_original",
                    "lifeform_coarse","native_status_original","station_codes","source_frequency_original","reconstructed_frequency"])
        for r in rows:
            w.writerow([r["taxon_original"],r["family_original"],r["lifeform_original"],
                        r["lifeform_coarse"],r["native_status_original"],";".join(r["station_codes"]),
                        r["printed_station_total"],r["parsed_station_total"]])

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--source-pdf",required=True)
    a.add_argument("--output-dir",required=True)
    p=a.parse_args()
    rows,errors,diag=parse_pdf(p.source_pdf)
    save(rows,diag,p.output_dir)

if __name__=="__main__":
    main()
