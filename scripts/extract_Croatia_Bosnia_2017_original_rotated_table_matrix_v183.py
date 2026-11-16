#!/usr/bin/env python3
"""Original-source 2017 Croatia/Bosnia railway stations: rotated vector Tab2.

Jasprica et al., Natura Croatica 26:271-303; doi:10.20302/NC.2017.26.21.
Recover only PDF-native source text (NOT OCR). Publisher survey of live track
areas at 15 stations. Source includes station numbers 1-15, life form and 
regional status; this is not experimental ballast equivalence across nations.
"""
import argparse,csv,hashlib,json,re
from pathlib import Path
PDF_URL="https://hrcak.srce.hr/file/282429"
SPECIES_RE=re.compile(r"^\*?[A-Z][A-Za-zÀ-ž\-]+\s+[a-z][a-zA-ZÀ-ž\-]+")
RECORD_RE=re.compile(r"^\s*\d{1,2}(?:\s*,\s*\d{1,2})*\s*,?\s*$")
LF={"T","H","P","Ch","G"}
FREQUENCIES={round(n*100/15) for n in range(1,16)}

def rows_from_original(pdf):
    import fitz
    raw=Path(pdf).read_bytes()
    if not raw.startswith(b"%PDF-"):raise ValueError("Not a PDF")
    doc=fitz.open(stream=raw,filetype="pdf")
    if len(doc)!=33:raise ValueError("Publisher PDF page count changed")
    out=[];page_reports=[];discarded=[]
    for pi in range(6,24):
        page=doc[pi]
        lines=[]
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines",[]):
                di=line.get("dir",(1,0))
                if abs(di[0])>.05 or abs(di[1]+1)>.05:continue
                txt="".join(sp.get("text","") for sp in line.get("spans",[])).strip()
                if not txt: continue
                x0,y0,x1,y1=line["bbox"]
                if x0<102 or x0>440:continue
                lines.append(dict(text=txt,x=(x0+x1)/2,y0=y0,y1=y1,x0=x0))
        taxon_anchors=[l for l in lines if l["y1"]>=560 and SPECIES_RE.match(l["text"])]
        page_reports.append({"page":pi+1,"n_taxon_anchor":len(taxon_anchors),
              "anchor_examples":[{"name":x["text"][:80],"x":round(x["x"],1)} for x in taxon_anchors[:4]]})
        for cand in taxon_anchors:
            # A row is a constant PDF x: all columns are vertically rotated.
            # Long species authority lines may use one extra glyph line.
            same=[l for l in lines if abs(l["x"]-cand["x"])<4.5]
            species=cand["text"].lstrip("*").strip()
            stationparts=[l["text"] for l in same if (400<=l["y0"] and l["y1"]<=507
                           and RECORD_RE.fullmatch(l["text"]))]
            rawstations=" ".join(stationparts).strip()
            if rawstations and not RECORD_RE.fullmatch(rawstations):
                rawstations=rawstations.replace(", ",",").replace(" ,",",").strip(",")
            vals=[int(v) for v in re.findall(r"\d{1,2}",rawstations)] if rawstations else []
            life=[l["text"] for l in same if 306<=l["y0"]<=333 and l["text"] in LF]
            family=[l["text"] for l in same if 344<=l["y0"]<=412 and re.fullmatch(r"[A-Za-z]+",l["text"])]
            frequency=[int(l["text"]) for l in same if 70<=l["y0"]<=110 and l["text"].isdigit()
                       and int(l["text"]) in FREQUENCIES]
            (out if vals and life else discarded).append({"taxon_original":species,"page":pi+1,
                        "source_station_numbers":sorted(set(vals)),
                        "source_station_list_text":rawstations,
                        "life_form":life[0] if life else "",
                        "family":family[0] if family else "",
                        "frequency_percent_source":frequency[0] if frequency else None,
                        "source_row_x":round(cand["x"],2)})
    metadata={"original_DOI":"10.20302/NC.2017.26.21","source_pdf_url":PDF_URL,
       "publisher_original_pdf_sha256":hashlib.sha256(raw).hexdigest(),
       "publisher_original_pdf_bytes":len(raw),"source_n_pdf_pages":len(doc),
       "n_record_taxon_candidates":len(out)+len(discarded),
       "n_discarded_unmatched_table_lines":len(discarded),
       "discarded_unmatched_table_lines":discarded[:12],
       "n_unique_candidate_taxa":len(set(o["taxon_original"] for o in out)),
       "n_station_lists_nonempty":sum(bool(o["source_station_numbers"]) for o in out),
       "n_life_form_nonempty":sum(bool(o["life_form"]) for o in out),
       "n_family_nonempty":sum(bool(o["family"]) for o in out),
       "n_frequency_nonempty":sum(o["frequency_percent_source"] is not None for o in out),
       "first_source_row_geometry":[{"text":l["text"],"x":round(l["x"],1),"y0":round(l["y0"],1),"y1":round(l["y1"],1)} for l in lines if abs(l["x"]-taxon_anchors[0]["x"])<5][:18] if taxon_anchors else [],
       "n_frequency_disagreements":sum(o["frequency_percent_source"] is not None
             and round(len(o["source_station_numbers"])*100/15)!=o["frequency_percent_source"]
             for o in out),
       "frequency_disagreements_detail":[{"taxon":o["taxon_original"],"page":o["page"],"printed_station_codes":o["source_station_numbers"],"printed_frequency":o["frequency_percent_source"],"calculated_frequency":round(len(o["source_station_numbers"])*100/15)} for o in out if o["frequency_percent_source"] is not None and round(len(o["source_station_numbers"])*100/15)!=o["frequency_percent_source"]],
       "malformed_station_ids":[{"taxon":o["taxon_original"],"records":o["source_station_numbers"]}
              for o in out if any(v<1 or v>15 for v in o["source_station_numbers"])][:15],
       "taxa_missing_station_preview":[{"taxon":o["taxon_original"],"page":o["page"],"x":o["source_row_x"]}
              for o in out if not o["source_station_numbers"]][:30],
       "page_reports":page_reports,
       "publication_expected_taxa":359,"publication_expected_stations":15,
       "source_count_gate":len(out)==359 and len(set(o["taxon_original"] for o in out))==359,
       "full_matrix_verified":False}
    # The original printed table itself has one detectable internal
    # inconsistency: Canna indica L. is listed at station 12 only, whereas its
    # printed frequency is 53% (8/15). Site incidence is taken from the explicit
    # station column; preserve the printed frequency as conflicting metadata.
    expected_original_anomaly = [{"taxon":"Canna indica L.", "page":9,
                                  "printed_station_codes":[12],
                                  "printed_frequency":53,
                                  "calculated_frequency":7}]
    metadata["source_frequency_anomaly_is_exactly_publisher_Canna_indica"] = (
        metadata["frequency_disagreements_detail"] == expected_original_anomaly)
    metadata["site_frequency_agreement_other_358_taxa"] = (
        metadata["n_frequency_nonempty"] == 359
        and metadata["n_frequency_disagreements"] == 1
        and metadata["source_frequency_anomaly_is_exactly_publisher_Canna_indica"])
    # Independent publisher text benchmarks (Natura Croatica, p. 276).
    station_occurrence_frequencies = {rr["taxon_original"]:len(rr["source_station_numbers"]) for rr in out}
    singletons = sum(c==1 for c in station_occurrence_frequencies.values())
    abundant = sum(c>=8 for c in station_occurrence_frequencies.values())
    universal=[k for k,c in station_occurrence_frequencies.items() if c==15]
    metadata["n_single_stations_original"]=singletons
    metadata["n_taxa_in_published_259_only_analysis"]=len(out)-singletons
    metadata["n_taxa_occurring_over_half_15_stations"]=abundant
    metadata["universal_station_taxa"]=universal
    # The publisher's prose says 100 taxa in one station / 259 elsewhere,
    # but its fully extracted table of 359 taxa has 99 / 260. Never modify
    # original station lists to force agreement with inconsistent prose.
    metadata["publisher_text_singletons_claim"]=100
    metadata["publisher_text_multisite_analysis_claim"]=259
    metadata["publisher_table_singletons_prose_discrepancy"]=(
        singletons==99 and len(out)-singletons==260)
    metadata["independent_publisher_crosscheck_72_frequent_Chondrilla_and_99_table_singletons"]=(
        metadata["publisher_table_singletons_prose_discrepancy"]
        and abundant==72 and universal == ["Chondrilla juncea L."])
    station_taxa={i:sum(i in rr["source_station_numbers"] for rr in out) for i in range(1,16)}
    form_taxa={}
    for rr in out:
        form_taxa[rr["life_form"]]=form_taxa.get(rr["life_form"],0)+1
    pair_jaccard=[]
    pair_function_profile_overlap=[]
    per_station_form_counts={}
    for i in range(1,16):
        per_station_form_counts[i]={
           form:sum(i in rr["source_station_numbers"] and rr["life_form"]==form for rr in out)
           for form in LF}
    for i in range(1,16):
        for j in range(i+1,16):
            shared=sum(i in rr["source_station_numbers"] and j in rr["source_station_numbers"] for rr in out)
            union=sum(i in rr["source_station_numbers"] or j in rr["source_station_numbers"] for rr in out)
            pair_jaccard.append(shared/union if union else None)
            ni=station_taxa[i];nj=station_taxa[j]
            pair_function_profile_overlap.append(sum(
                min(per_station_form_counts[i][form]/ni,per_station_form_counts[j][form]/nj)
                for form in LF))
    metadata["station_taxa_counts"] = station_taxa
    metadata["life_form_taxa_counts"] = form_taxa
    metadata["n_pairwise_station_jaccard_descriptive_nonindependent"] = len(pair_jaccard)
    metadata["station_jaccard_min"]=round(min(pair_jaccard),4)
    metadata["station_jaccard_median"]=round(sorted(pair_jaccard)[len(pair_jaccard)//2],4)
    metadata["station_jaccard_max"]=round(max(pair_jaccard),4)
    metadata["per_station_original_lifeform_counts"]=per_station_form_counts
    metadata["station_lifeform_profile_overlap_min"]=round(min(pair_function_profile_overlap),4)
    metadata["station_lifeform_profile_overlap_median"]=round(sorted(pair_function_profile_overlap)[len(pair_function_profile_overlap)//2],4)
    metadata["station_lifeform_profile_overlap_max"]=round(max(pair_function_profile_overlap),4)
    metadata["lifeform_overlap_not_comparable_directly_to_Jaccard_or_regional_pool_null"]=True
    metadata["station_jaccard_is_within_ONE_study_not_cross_region_homogenization"]=True
    metadata["full_matrix_verified"]=(metadata["source_count_gate"]
         and metadata["n_station_lists_nonempty"]==359
         and metadata["n_life_form_nonempty"]==359
         and metadata["site_frequency_agreement_other_358_taxa"]
         and metadata["independent_publisher_crosscheck_72_frequent_Chondrilla_and_99_table_singletons"]
         and not metadata["malformed_station_ids"])
    return out,metadata

def save(pdf,outdir):
    rows,a=rows_from_original(pdf)
    p=Path(outdir);p.mkdir(parents=True,exist_ok=True)
    (p/"Croatia_Bosnia_2017_native_rotated_table_source_QA.json").write_text(
         json.dumps(a,ensure_ascii=False,indent=2),encoding="utf8")
    print(json.dumps({k:v for k,v in a.items() if k!="page_reports"},ensure_ascii=False,indent=2))
    print("PAGE REPORTS",json.dumps(a["page_reports"],ensure_ascii=False))
    if not a["full_matrix_verified"]:
         raise ValueError("Original vector PDF parsing did not pass 359-taxon station-list/trait frequency source gates")
    with (p/"Croatia_Bosnia_2017_15station_359taxa_source_traits.csv").open("w",newline="",encoding="utf8") as f:
         w=csv.writer(f);w.writerow(["taxon_original","station_codes","life_form","family","page"])
         for r in rows:w.writerow([r["taxon_original"],";".join(map(str,r["source_station_numbers"])),
                                    r["life_form"],r["family"],r["page"]])
    with (p/"Croatia_Bosnia_2017_15station_species_incidence.csv").open("w",newline="",encoding="utf8") as f:
         w=csv.writer(f);w.writerow(["study","rail_station_number","taxon_original","presence",
                                      "life_form","family","pdf_page"])
         for r in rows:
             for station in range(1,16):
                 w.writerow(["Jasprica2017",station,r["taxon_original"],
                             int(station in r["source_station_numbers"]),
                             r["life_form"],r["family"],r["page"]])

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--original-pdf",required=True);p.add_argument("--output-dir",required=True)
    z=p.parse_args();save(z.original_pdf,z.output_dir)
