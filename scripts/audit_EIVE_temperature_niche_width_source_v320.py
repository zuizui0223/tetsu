#!/usr/bin/env python3
"""v320 inspect genuine open EIVE 1.0 SM08 niche-position/width XLSX.

Does not infer taxonomy, guess column meanings or republish the source rows.
Expected publisher Zenodo md5:0f0ce19cd2a781eeed027d8b0acadbcb.
Read-only OOXML inspection via Python stdlib; result is schema-only metadata.
"""
import argparse
import hashlib
import json
import posixpath
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

EXPECTED_MD5="0f0ce19cd2a781eeed027d8b0acadbcb"
FILE="EIVE_Paper_1.0_SM_08.xlsx"
NS={"s":"http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
RELS="{http://schemas.openxmlformats.org/package/2006/relationships}"
S="{"+NS["s"]+"}"

def colnum(cell_ref):
    m=re.match(r"^([A-Z]+)[0-9]+$",str(cell_ref))
    if not m: raise ValueError("Invalid OOXML cell address")
    n=0
    for c in m.group(1): n=n*26+ord(c)-64
    return n-1

def cell_value(cell,strings):
    t=cell.attrib.get("t")
    if t=="inlineStr":
        return "".join(x.text or "" for x in cell.iter(S+"t"))
    v=cell.find(S+"v")
    if v is None: return None
    raw=v.text
    if t=="s":
        i=int(raw)
        if not 0<=i<len(strings): raise ValueError("Bad shared-string index")
        return strings[i]
    if t=="e": raise ValueError("Workbook error cell")
    return raw

def audit(path):
    data=Path(path).read_bytes()
    md5=hashlib.md5(data).hexdigest()
    if md5!=EXPECTED_MD5:
        raise ValueError("Not the exact publisher EIVE SM08 file")
    if not zipfile.is_zipfile(path): raise ValueError("Not a valid OOXML workbook")
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None: raise ValueError("ZIP CRC failed")
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        relroot=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rels={r.attrib["Id"]:r.attrib["Target"] for r in relroot
              if r.tag.startswith(RELS)}
        if "xl/sharedStrings.xml" in z.namelist():
            shared=ET.fromstring(z.read("xl/sharedStrings.xml"))
            strings=["".join(t.text or "" for t in row.iter(S+"t"))
                     for row in shared.findall(S+"si")]
        else:
            strings=[]
        sheets=[]
        for sh in wb.findall(".//s:sheet",NS):
            key=sh.attrib["{"+NS["r"]+"}id"]
            target=rels[key]
            file=(target.lstrip("/") if target.startswith("/")
                  else posixpath.normpath(posixpath.join("xl",target)))
            if not file.startswith("xl/worksheets/"):
                raise ValueError("Workbook sheet target outside expected namespace")
            root=ET.fromstring(z.read(file))
            rows=root.findall(".//s:sheetData/s:row",NS)
            preview=[]
            for row in rows[:9]:
                cells={str(colnum(c.attrib["r"])):str(cell_value(c,strings))[:110]
                       for c in row.findall(S+"c") if cell_value(c,strings) is not None}
                preview.append({"sheet_row":row.attrib.get("r"),"cells":cells})
            sheets.append({"name":sh.attrib["name"],"xml":file,
                           "sheet_xml_rows":len(rows),
                           "first_9_row_schema_preview":preview})
    return {"version":"v320","source":"Dengler et al. 2023 EIVE 1.0 supplemental material 8",
            "source_doi":"10.5281/zenodo.7534792",
            "original_filename":FILE,
            "original_md5":md5,
            "original_sha256":hashlib.sha256(data).hexdigest(),
            "original_bytes":len(data),
            "sheets":sheets,"source_species_rows_not_redistributed":True,
            "position_and_niche_width_impact_not_estimated":True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True)
    ap.add_argument("--output",required=True)
    args=ap.parse_args()
    result=audit(args.input)
    dest=Path(args.output)
    dest.parent.mkdir(exist_ok=True,parents=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
