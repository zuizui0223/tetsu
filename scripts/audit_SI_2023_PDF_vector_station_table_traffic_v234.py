#!/usr/bin/env python3
"""v23.4 original-2023 operator PDF cell-based station/km recovery, no invented sites.

True 2023 Slovenian infrastructure-manager Annex 2B PDF, 16 pages, SHA pinned.
Re-read the publisher table using PDF VECTOR LINE GEOMETRY (pdfplumber
find_tables -> extract) to repair multi-line station labels ("Ljubljana
Zalog", "Ljubljana Šiška"), which the original PDF extract_text layout
can split or falsely prefix. Then apply EXACT original v23.1 national
rail line + two endpoint station km crosswalk to Eurostat 2020 real
directionally observed goods/passenger trains; collapse opposing
direction counts to one unique physical interval. Map original 89
botanical sites only when line+source chainage uniquely contains site.

The true railway vegetation was field-surveyed in 2022/2023; 2020 train
counts are HISTORICAL and are not assigned field-year contemporaneous
causal status. No ecological regression is performed in this source
coverage audit. No actual source 2022/23 freight counts are fabricated.
"""
from __future__ import annotations
import argparse,hashlib,json,math,re
from collections import Counter,defaultdict
from pathlib import Path
from audit_SI_2023_official_station_chainage_2020_eurostat_v231 import (
    extract_authoritative_stations, normalize,PDF_SOURCE,RE_LINE,
    source_vegetation_sites,eurostat_2020_segments,segment_route_station_join,
    collapse_exact_reverse_direction_segments,assign_segment
)
from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx

PDF_SHA="83f12dfdb91784a653b9a39a8d927a027e7eeb5d8af18ce17b850b5318796425"
RE_FLOAT3=re.compile(r"^\d{1,4}\.\d{3}$")
RE_INT=re.compile(r"^\d{2,5}$")
SOURCE_LINES=frozenset(("10","20","21","31","32","34","41","42","44","50","64","70","72","80","81"))

def text_value(v):
    return " ".join(str(v or "").split())

def classify_cell_row(row,current):
    """Only accept strict four first official station columns; never
    convert track length / platform metres into railway station km."""
    cells=[text_value(q) for q in row]
    joined=" ".join(cells)
    header=RE_LINE.search(joined)
    if header:
        return None,header.group(1),"source_line_heading"
    if current is None:return None,current,"no_source_national_line"
    if len(cells)<4:return None,current,"fewer_than_four_station_columns"
    name,low,high,length=cells[:4]
    if not name or not RE_FLOAT3.fullmatch(low) or not RE_FLOAT3.fullmatch(high):
        return None,current,"no_four_cell_station_record"
    length=length.replace(" ","")
    if not RE_INT.fullmatch(length):
        return None,current,"original_station_length_not_int"
    lo,hi,ln=float(low),float(high),int(length)
    if not(0<=lo<hi<1000 and ln>=20 and abs(1000*(hi-lo)-ln)<=15):
        return None,current,"length_km_not_source_consistent"
    if len(name)>90 or re.search(r'\b(STATION|INITIAL|FINAL|LENGTH|EQUIPMENT)\b',name,re.I):
        return None,current,"source_name_not_station"
    # A new station must not be accepted based on a fractional platform
    # field, station name must occur in the first original table column.
    return {'national_line':current,'name':name,
            'name_normalized':normalize(name),
            'station_start_km':lo,'station_end_km':hi,
            'station_midpoint_km':(lo+hi)/2},current,"strict_official_station_table_row"

def pdf_station_tables(pdfpath):
    import pdfplumber
    data=[];sections={};page_qa=[];confusions=Counter()
    current=None
    with pdfplumber.open(pdfpath) as pdf:
        if len(pdf.pages)!=16:raise ValueError("Real 2023 government source was not 16p")
        for n,p in enumerate(pdf.pages):
            if n<2:continue
            tables=p.find_tables(table_settings={
                "vertical_strategy":"lines",
                "horizontal_strategy":"lines",
                "intersection_tolerance":5,
                "snap_tolerance":3,
                "join_tolerance":3,
                "edge_min_length":8,
            })
            previews=[];records=0
            for tb in sorted(tables,key=lambda t:t.bbox[1]):
                cells=tb.extract(x_tolerance=1.5,y_tolerance=3)
                if len(previews)<3:
                    previews.append({'table_rows':len(cells),
                      'table_columns_max':max((len(row) for row in cells),default=0),
                      'first_four_original_vector_rows':[
                         [text_value(q)[:100] for q in row[:6]] for row in cells[:4]]})
                for row in cells:
                    station,current,why=classify_cell_row(row,current)
                    confusions[why]+=1
                    if station is None:continue
                    station['source_pdf_page_zero_index']=n
                    records+=1;data.append(station)
            page_qa.append({'source_page_zero_index':n,
                'vector_tables_detected':len(tables),
                'exact_four_column_station_records':records,
                'example_original_pdf_table_cells':previews})
    return data,{'source_vector_table_station_records':len(data),
          'source_line_ids_covered':sorted({q['national_line'] for q in data},key=int),
          'original_sampled_15_lines_covered':sorted({q['national_line'] for q in data}&SOURCE_LINES,key=int),
          'per_line_source_vector_records':dict(Counter(q['national_line'] for q in data)),
          'vector_cells_parse_rejection_reason_counts':dict(confusions),
          'original_page_vector_cell_QA':page_qa}

def merge_old_new(old,new):
    """Same national line and EXACT original low/high station km are required.
    Conflicting non-nested station names produce unknown, not fake match."""
    by=defaultdict(list)
    for s in old:by[s['national_line'],round(s['station_start_km'],3),
                    round(s['station_end_km'],3)].append(('line_text',s))
    for s in new:by[s['national_line'],round(s['station_start_km'],3),
                    round(s['station_end_km'],3)].append(('table_vector',s))
    result=[];qa=Counter();bad=[]
    for k,entries in sorted(by.items()):
        byname=defaultdict(list)
        for source,s in entries:byname[s['name_normalized']].append((source,s))
        source_names=sorted(byname)
        if len(source_names)==1:
            if any(source=="table_vector" for source,s in entries):
                result.append(next(s for source,s in entries if source=="table_vector"))
                qa['station_positions_from_source_table_exact']+=1
            else:
                result.append(entries[0][1]);qa['station_positions_text_only']+=1
            continue
        # Wrapped label is a provable extension of the originally cropped
        # label; reject arbitrary same-km but different source names.
        extension=[(small,big) for small in source_names for big in source_names
                   if small!=big and (big.endswith(" "+small) or big.startswith(small+" "))]
        if len(source_names)==2 and extension:
            longest=max(source_names,key=len)
            chosen=byname[longest]
            if len(chosen)>=1:
                result.append(chosen[0][1]);qa['verified_multiline_wrapped_station_name_repaired']+=1
                continue
        qa['conflicted_source_table_text_label_excluded']+=1
        if len(bad)<25:bad.append({'line':k[0],'low_km':k[1],
              'high_km':k[2],'conflicting_source_names':source_names})
    return result,dict(qa)|{'source_conflicting_name_examples_no_join':bad}

def analysis(pdfpath,original):
    raw=Path(pdfpath).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PDF_SHA:
        raise ValueError("True 2023 operator station PDF changed")
    baseline,sourceqa=extract_authoritative_stations(pdfpath)
    table,tableqa=pdf_station_tables(pdfpath)
    combined,mergeqa=merge_old_new(baseline,table)
    vegetation,botanicalqa=source_vegetation_sites(original)
    segment,euroqa=eurostat_2020_segments()
    # Exact same source biology and traffic input, only station label parsing
    # changes; do not force a larger matching cohort in the face of ambiguity.
    base_matches,basejoinqa=segment_route_station_join(baseline,segment)
    basephysical,basepqa=collapse_exact_reverse_direction_segments(base_matches)
    basemapped,basecoverage=assign_segment(vegetation,basephysical)
    if len(basemapped)!=27 or basecoverage['n_unique_2020_traffic_segments_mapped']!=9:
        raise ValueError("Original certified v23.1 baseline exact 27 sites and nine physical intervals drifted")
    matches,joinqa=segment_route_station_join(combined,segment)
    physical,phyqa=collapse_exact_reverse_direction_segments(matches)
    mapped,newcoverage=assign_segment(vegetation,physical)
    old={(v['source_site'],v['eurostat_segment_code'],
          v['historical_2020_goods_trains_per_year']) for v in basemapped}
    new={(v['source_site'],v['eurostat_segment_code'],
          v['historical_2020_goods_trains_per_year']) for v in mapped}
    lost=old-new;added=new-old
    if lost:
        raise ValueError('Source vector table changed an already verified true station-to-traffic link; no silent replacement: '+str(len(lost)))
    if len(table)<80:
        raise ValueError('2023 PDF vector table parsed fewer than 80 authentic named stations')
    return {'version':'v23.4',
      'status':'SOURCE_AUTHENTICATED_2023_STATION_PDF_VECTOR_TABLE_RECOVERY',
      'operator_2023_annex_pdf_sha256':PDF_SHA,
      'source_original_2020_eurostat_json_sha':euroqa['source']['sha256_original_json_response'],
      'original_source_89_vegetation_2022_2023_sites':len(vegetation),
      'original_text_parser_station_positions':len(baseline),
      'new_true_vector_table_station_positions':len(table),
      'combined_nonconflicted_true_operator_station_positions':len(combined),
      'source_vector_table_diagnostics':tableqa,
      'combined_old_new_original_source_name_conflicts':mergeqa,
      'baseline_source_27_flora_site_crosswalk':basecoverage,
      'recovered_site_crosswalk_after_real_PDF_table_parsing':newcoverage,
      'new_certified_source_site_to_2020_train_segment_pairs':len(added),
      'old_certified_source_site_to_2020_train_segment_pairs_lost':len(lost),
      'source_unique_two_direction_physical_segments_before':len(basephysical),
      'source_unique_two_direction_physical_segments_after':len(physical),
      'original_2020_freight_train_directions_are_not_independent_lines':True,
      '2020_TRAIN_COUNT_NOT_2022_2023_EXPOSURE':True,
      'no_fitted_ecological_traffic_effect_claimed':True,
      'unmatched_sites_unknown_traffic_not_zero':True}
def main():
    p=argparse.ArgumentParser()
    for n in ('official_station_pdf','original_xlsx','source_audit','output'):
        p.add_argument('--'+n.replace('_','-'),required=True)
    args=p.parse_args()
    org=extract_original_xlsx(args.original_xlsx,args.source_audit)
    result=analysis(args.official_station_pdf,org)
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    slim={k:result[k] for k in ('status','original_text_parser_station_positions',
      'new_true_vector_table_station_positions',
      'combined_nonconflicted_true_operator_station_positions',
      'combined_old_new_original_source_name_conflicts',
      'baseline_source_27_flora_site_crosswalk',
      'recovered_site_crosswalk_after_real_PDF_table_parsing',
      'new_certified_source_site_to_2020_train_segment_pairs',
      'source_unique_two_direction_physical_segments_before',
      'source_unique_two_direction_physical_segments_after')}
    # Diagnostics shows original pdf table cells and labels but never the
    # original plant by site incidence matrix or field site GPS.
    print('V234_TRUE_PDF_TABLE_CROSSWALK='+json.dumps(slim,ensure_ascii=False,separators=(',',':')))
if __name__=='__main__':main()
