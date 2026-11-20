"""Extract published matrix cells and plot metadata; no inference from blank prose."""
import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
import fitz


def lines(page):
    out = []
    for block in page.get_text('dict')['blocks']:
        for line in block.get('lines', []):
            text = ''.join(s['text'] for s in line['spans']).strip()
            box = fitz.Rect(line['bbox']) * page.rotation_matrix
            if text:
                out.append((box.x0, box.y0, box.x1, box.y1, text))
    return out


def extract(path):
    doc = fitz.open(path)
    records, taxa, sites, metadata = [], set(), set(), {}
    for page_no in range(387, len(doc)):
        page = doc[page_no]
        if 'Matrix of Species Found Within Each Site' not in page.get_text():
            continue
        rows = lines(page)
        headers = [r for r in rows if re.fullmatch(r'[A-Z]{2,3}\d+[A-Za-z]?', r[4])]
        assert headers, page_no
        left = min(r[0] for r in headers)
        header_y = min(r[1] for r in headers)
        names = [r for r in rows if r[0] < left - 20 and r[1] > header_y + 4
                 and re.match(r'^[A-Z][a-z]+\s', r[4])]
        for name in names:
            taxa.add(name[4])
            for h in headers:
                sites.add(h[4])
            for cell in rows:
                if abs(cell[1] - name[1]) < 2 and cell[0] >= left - 20:
                    if re.fullmatch(r'\+|[\d.]+%', cell[4]):
                        h = min(headers, key=lambda h: abs((h[0]+h[2])/2-(cell[0]+cell[2])/2))
                        records.append(dict(site=h[4], taxon=name[4], cover_raw=cell[4], pdf_page=page_no+1))
    for page_no in range(193, 380):
        text = doc[page_no].get_text()
        m = re.search(r'Site\s+([A-Z]{2,3}\d+[A-Za-z]?)', text)
        if not m:
            continue
        site = m[1]
        def grab(pattern):
            m = re.search(pattern, text, re.S)
            return m.group(1).strip() if m else None
        coords = re.search(r'MGA\s+(\d+)\s+(\d+)\s*mE\s+(\d+)\s*mN', text)
        dims = re.search(r'Type\s+(?:[QR]\s+)?(\d+)\s*x\s*(\d+)\s*m?', text)
        metadata[site] = dict(site=site, pdf_page=page_no+1,
            zone=int(coords[1]) if coords else None,
            easting=int(coords[2]) if coords else None,
            northing=int(coords[3]) if coords else None,
            type=grab(r'Type\s+([QR])'),
            area_m2=int(dims[1])*int(dims[2]) if dims else None,
            habitat=grab(r'Habitat\s+(.*?)\s+Soil'),
            condition=grab(r'\nVeg\s+(.*?)\s+Fire'),
            location=grab(r'Location\s+(.*?)\s+MGA'))
    sheets = []
    site = None
    for page_no in range(193, 380):
        page = doc[page_no]
        text = page.get_text()
        if 'Appendix G' not in text:
            continue
        found = re.search(r'Site\s+([A-Z]{2,3}\d+[A-Za-z]?)', text)
        if found:
            site = found[1]
        rr = lines(page)
        for cover in rr:
            if not (320 < cover[0] < 360 and re.fullmatch(r'\+|[\d.]+%', cover[4])):
                continue
            row = [r for r in rr if abs(r[1]-cover[1]) < 1.5]
            names = [r[4] for r in row if r[0] < 120 and re.match(r'^[A-Z][a-z]+\s',r[4])]
            if len(names) != 1:
                raise ValueError(('Unresolved sheet row',page_no+1,cover,row))
            notes = ' '.join(r[4] for r in rr if r[0] > 480 and abs((r[1]+r[3]-cover[1]-cover[3])/2)<4)
            height = ' '.join(r[4] for r in row if 390 < r[0] < 435)
            sheets.append(dict(site=site,taxon=names[0],cover_raw=cover[4],height_raw=height,
                notes=notes,out_of_plot=bool(re.search(r'\bOP\b',notes)),pdf_page=page_no+1))
    pairs = [(r['site'], r['taxon']) for r in records]
    duplicates = [p for p,c in Counter(pairs).items() if c>1]
    summary = dict(matrix_sites=len(sites), matrix_taxa=len(taxa), incidences=len(records),
        duplicate_cells=duplicates, metadata_sites=len(metadata),
        matrix_without_metadata=sorted(sites-set(metadata)),
        metadata_without_matrix=sorted(set(metadata)-sites),
        area_counts=dict(Counter(str(m['area_m2']) for m in metadata.values())),
        type_counts=dict(Counter(str(m['type']) for m in metadata.values())),
        missing_coords=sum(m['easting'] is None for m in metadata.values()))
    return dict(source_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        summary=summary, matrix=records, metadata=list(metadata.values()), taxa=sorted(taxa),sheet_records=sheets)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--pdf', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    result = extract(a.pdf)
    Path(a.output).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(result['summary'], indent=2))
