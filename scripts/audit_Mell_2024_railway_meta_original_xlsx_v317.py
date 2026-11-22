#!/usr/bin/env python3
"""First-party public Mell et al. 2024 meta-dataset source eligibility audit.

Read XLSX workbook via stdlib ZIP+ElementTree, not an opaque spreadsheet
package. Do not treat source case rows as independent studies. The output
only contains structural metadata, column headings and a bounded set of rows
with *rail* literals; never assign effect directions before source checks.
"""
from __future__ import annotations
import argparse,hashlib,json,re,zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import Counter,defaultdict

X_NS="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
R_NS="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
P_NS="http://schemas.openxmlformats.org/package/2006/relationships"
EXTRACT_SOURCE="10.1186/s13750-024-00328-3"
META_SOURCE="https://link.springer.com/article/10.1186/s13750-024-00328-3"
MELL8="13750_2024_328_MOESM8_ESM.xlsx"
MELL5="13750_2024_328_MOESM5_ESM.xlsx"

def _col(a):
    x=re.match(r"^([A-Z]+)[1-9][0-9]*$",a)
    if not x:raise ValueError("Source XLSX cell coordinate unsupported: "+repr(a))
    k=0
    for c in x.group(1):k=k*26+ord(c)-64
    return k-1

def source_original_workbook(path:Path):
    raw=path.read_bytes()
    if len(raw)<15000 or not raw.startswith(b"PK\x03\x04"):
        raise ValueError("Original public source is not an XLSX file")
    out={}
    with zipfile.ZipFile(path) as z:
        names=set(z.namelist())
        if not ("xl/workbook.xml" in names and "xl/_rels/workbook.xml.rels" in names):
            raise ValueError("Genuine XLSX lacks required workbook relationship info")
        shared=[]
        if "xl/sharedStrings.xml" in names:
            root=ET.fromstring(z.read("xl/sharedStrings.xml"))
            for x in root.findall(f"{{{X_NS}}}si"):
                shared.append("".join(v.text or "" for v in x.iter(f"{{{X_NS}}}t")))
        workbook=ET.fromstring(z.read("xl/workbook.xml"))
        sh=workbook.findall(f"{{{X_NS}}}sheets/{{{X_NS}}}sheet")
        rels=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rmap={e.get("Id"):e.get("Target") for e in rels.findall(f"{{{P_NS}}}Relationship")}
        for sheet in sh:
            title=sheet.get("name")
            loc=rmap.get(sheet.get(f"{{{R_NS}}}id"),"")
            loc="xl/"+loc.lstrip("/").removeprefix("xl/")
            if loc not in names:raise ValueError("Worksheet original target is not readable: "+title)
            root=ET.fromstring(z.read(loc))
            rows=[]
            for row in root.findall(f"{{{X_NS}}}sheetData/{{{X_NS}}}row"):
                cells={}
                for c in row.findall(f"{{{X_NS}}}c"):
                    k=_col(c.get("r",""))
                    typ=c.get("t","")
                    val=c.find(f"{{{X_NS}}}v")
                    v=val.text if val is not None else None
                    if typ=="s":
                        if v is not None:v=shared[int(v)]
                    elif typ=="inlineStr":
                        v="".join(t.text or "" for t in c.iter(f"{{{X_NS}}}t"))
                    elif typ in ("n",""):
                        if v is not None:
                            try:v=float(v)
                            except ValueError:pass
                    elif typ=="b":
                        v=v=="1"
                    elif typ in ("e","str"):
                        pass
                    else:raise ValueError("Unrecognized original source XLSX type "+str(typ))
                    if v is not None:cells[k]=v
                rows.append(cells)
            out[title]=rows
    return out,{"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest(),
                "zip_members":len(names),"n_original_sheets":len(out)}


def audit(source:Path):
    sheets,qa=source_original_workbook(source)
    profiles=[]
    for name,records in sheets.items():
        candidates=[]
        for index,record in enumerate(records):
            matches=[(col,str(value)[:110]) for col,value in record.items()
                     if isinstance(value,str) and re.search(r"\brail(?:way|ways|road|track)?\b",value,re.I)]
            if matches:
                candidates.append({"source_row_1based":index+1,
                                   "matches":matches[:10],
                                   "source_publication_candidate":str(record.get(0,""))[:100],
                                   "selected_ambient_columns":[[i,str(record[i])[:75]] for i in range(min(8,max(record,default=0)+1)) if i in record]})
        profiles.append({
            "sheet":name,
            "source_nonempty_rows":len(records),
            "max_source_column_index":max((max(r) for r in records if r),default=-1),
            "source_initial_3_rows":[[[col,str(value)[:90]] for col,value in list(sorted(row.items()))[:30]] for row in records[:3]],
            "n_rows_with_literal_rail_word_in_any_column":len(candidates),
            "bounded_literal_rail_rows_NOT_assumed_effect_size":candidates[:25],
            "source_first_cell_role_caveat":"Some spreadsheets contain notes/multi-header rows; do not classify by any single unvalidated keyword."
        })
    if not profiles:raise ValueError("No source meta-analysis worksheets found")
    if source.name!=MELL8:
        raise ValueError("Wrong supplementary original effect size file")
    if qa["sha256"]!="0fd8ff9922c1138d5abd5289e743058b10be3b52c73d4d338cc9b9d2505def17":
        raise ValueError("Previously authenticated original publisher Supplement 8 SHA changed")
    if "File" not in sheets or len(sheets["File"])!=206:
        raise ValueError("Original source 205 case rows changed")
    records=sheets["File"]
    header=records[0]
    if (header.get(0),header.get(1),header.get(3),header.get(14))!=(
            "Study.ID","Author_date","Case.ID","LTI"):
        raise ValueError("Unexpected original publication/table reference headers")
    category=Counter()
    rail=[]
    for no,row in enumerate(records[1:],start=2):
        lti=row.get(14)
        category[str(lti)]+=1
        if str(lti).lower()=="railway":
            effect={str(header.get(col,f"source_column_{col}")):(
                None if row.get(col) in (None,"NA","") else row.get(col))
                    for col in range(max(header)+1)}
            effect["original_excel_row_1based"]=no
            rail.append(effect)
    if len(rail)!=3 or category["Railway"]!=3:
        raise ValueError("Original publication exact three Railway Q2 cases not reproduced")
    rail_studies=sorted({str(r["Author_date"]) for r in rail})
    if len(rail_studies)!=2:
        raise ValueError("Original railway cases do not form two primary publication IDs")
    return {
        "n_actual_source_Q2_effect_rows":len(records)-1,
        "source_LTI_case_counts":dict(sorted(category.items())),
        "source_unique_primary_author_date_IDs_all_LTI":len({str(r.get(1)) for r in records[1:]}),
        "n_railway_effect_rows":len(rail),
        "n_railway_distinct_publications":len(rail_studies),
        "railway_original_publication_identifiers":rail_studies,
        "source_original_railway_effect_metadata_PUBLIC":rail,
        "publisher_railway_cases_are_not_independent_marginals":True,
        "cannot_claim_rail_specific_climate_interaction_from_two_original_papers":True,

        "version":"v317","official_publication_doi":EXTRACT_SOURCE,
        "original_publisher_record":META_SOURCE,
        "original_source_xlsx":source.name,
        "source_source_sha_and_workbook":qa,
        "sheets":profiles,
        "study_cases_reported_by_2024_review":205,
        "independent_primary_studies_reported_by_review":47,
        "railway_comparisons_reported_rounded_percent":1,
        "rail_only_effect_aggregation_permitted":False,
        "status":"SOURCE_STRUCTURAL_GATE_ONLY_NO_EFFECT_INFERENCE",
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--xlsx",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    x=p.parse_args()
    data=audit(x.xlsx)
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({
        "status":data["status"],
        "original_source":data["original_source_xlsx"],
        "sha":data["source_source_sha_and_workbook"]["sha256"],
        "sheets": [{k:q[k] for k in ("sheet","source_nonempty_rows",
                    "max_source_column_index","n_rows_with_literal_rail_word_in_any_column")}
                    for q in data["sheets"]]
    },ensure_ascii=False,indent=2))
if __name__=="__main__":main()
