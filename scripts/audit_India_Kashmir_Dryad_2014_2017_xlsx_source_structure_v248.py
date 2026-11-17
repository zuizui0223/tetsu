#!/usr/bin/env python3
"""v24.8. Read ONLY the authentic 122,187-byte Rashid et al. Kashmir
Himalaya 2014/2017 31 railway T transects original Dryad Matrix_plot.xlsx.

First-party Dryad dataset: DOI 10.5061/dryad.0p2ngf20h
(dataset 60988, version 99735, file 546003), file SHA256
529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474.
This probe introspects only workbook XML structure/first original row
labels, instead of automatically attributing species to plot types.

Do not mistake source author aggregate quadrat cover vs species by
site, or 31 transects×4 zones×2 years as automatically 248 complete
independent 50×2m plot records. Do not impute absent cells as absence.
No ecological effect or network connectivity estimator at this step.
"""
from __future__ import annotations
import argparse,hashlib,json,re,zipfile,xml.etree.ElementTree as ET
from pathlib import Path
from collections import Counter

EXPECTED_SHA="529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474"
NAMESPACES={"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main",
            "r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships",
            "p":"http://schemas.openxmlformats.org/package/2006/relationships"}
NS=NAMESPACES["m"]

def ss(z):
    f="xl/sharedStrings.xml"
    if f not in z.namelist():return []
    root=ET.fromstring(z.read(f))
    return ["".join(c.text or "" for c in si.iter(f"{{{NS}}}t")) for si in root.findall("m:si",NAMESPACES)]

def display_row(node,shared,max_cells=30):
    result=[]
    for cell in list(node.findall("m:c",NAMESPACES))[:max_cells]:
        val=cell.find("m:v",NAMESPACES)
        text=val.text if val is not None else ""
        typ=cell.get("t","")
        if typ=="s" and text:
            idx=int(text)
            if idx>=len(shared):raise ValueError("Original shared string index exceeds source table")
            text=shared[idx]
        elif typ=="inlineStr":
            text="".join(x.text or "" for x in cell.iter(f"{{{NS}}}t"))
        result.append({"source_cell":cell.get("r"),"text":str(text)[:115],"source_type":typ})
    return result

def audit(src):
    path=Path(src);raw=path.read_bytes()
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=EXPECTED_SHA:
        raise ValueError("Authentic Kashmir author original Dryad Matrix_plot.xlsx SHA changed")
    if not raw.startswith(b"PK\x03\x04"):raise ValueError("Original Dryad XLSX is not a ZIP XML workbook")
    result={"version":"v24.8","source_sha256":digest,"source_original_bytes":len(raw),
            "source_dryad_file_id":546003,"source_first_party_DOI":"10.5061/dryad.0p2ngf20h",
            "source_Dryad_original_version_id":99735,"source_original_workbook_sheets":[]}
    with zipfile.ZipFile(path) as z:
        workbook=ET.fromstring(z.read("xl/workbook.xml"))
        rels=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rid_to_file={}
        for r in rels.findall("p:Relationship",NAMESPACES):
            target=r.get("Target") or ""
            if "worksheet" in (r.get("Type") or ""):
                rid_to_file[r.get("Id")]=("xl/"+target.lstrip("/").removeprefix("xl/"))
        shared=ss(z)
        result["source_original_shared_strings_n"]=len(shared)
        for sh in workbook.findall("m:sheets/m:sheet",NAMESPACES):
            rid=sh.get("{"+NAMESPACES["r"]+"}id")
            filename=rid_to_file.get(rid)
            if filename is None or filename not in z.namelist():
                raise ValueError("Dryad source worksheet relative URI absent")
            sheet=ET.fromstring(z.read(filename))
            dim=sheet.find("m:dimension",NAMESPACES)
            data=sheet.find("m:sheetData",NAMESPACES)
            rows=data.findall("m:row",NAMESPACES) if data is not None else []
            preview=[{"original_source_row":r.get("r"),
                "source_nonempty_cells":len(r.findall("m:c",NAMESPACES)),
                "actual_source_cells_first35":display_row(r,shared,35)}
                 for r in rows[:9]]
            allstrings=Counter()
            for row in rows:
                for c in row.findall("m:c",NAMESPACES):
                    typ=c.get("t")
                    if typ in ("s","inlineStr"):
                        allstrings[typ]+=1
            result["source_original_workbook_sheets"].append({
                "original_sheet_name":sh.get("name"),
                "source_dimension_ref":dim.get("ref") if dim is not None else None,
                "source_stored_XML_row_nodes":len(rows),
                "source_stored_XML_cells_total":sum(len(r.findall("m:c",NAMESPACES)) for r in rows),
                "source_text_cell_types":dict(allstrings),
                "original_first_rows_column_header_preview_not_reconstructed_incidence_matrix":preview,
                "source_original_XML_path":filename})
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument("--source",required=True);p.add_argument("--output",required=True)
    a=p.parse_args();report=audit(a.source)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf8")
    print("V248_KASHMIR_ORIGINAL_XLSX_STRUCTURE="+json.dumps(report,ensure_ascii=False,separators=(",",":")),flush=True)
if __name__=="__main__":main()
