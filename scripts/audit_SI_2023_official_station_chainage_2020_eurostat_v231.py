#!/usr/bin/env python3
"""v23.1: Genuine 2023 Slovenian railway-station chainage vs Eurostat 2020 traffic join opportunity.

Source A: Slovenian rail infrastructure official Network Statement 2023
ANNEX 2B (16 page station table):
https://infrastruktura.sz.si/wp-content/uploads/2023/02/NS_2023_1.2_Annex_2B_Stations.pdf
Source B: original 2022–23 source vegetation by Šabić et al. (2026) SHA-locked
publisher XLSX (89 true ballast shoulder+toe sites and true national rail
route number and 100m transect start/end km).
Source C: Eurostat 2015/2020 official goods/passenger trains by network segment.

Requires physical source station name AND source national route kilometer
interval. Eurostat segment ID is NOT the national 2-digit route number.
No geodesic interpolation of Eurostat station endpoints, and no copying of
fake site-level passenger or freight train counts.
No causal estimates from temporally mismatched 2020 traffic and 2022/23
flora. 2020 pandemic specifically excluded from contemporary exposure.
"""
from __future__ import annotations
import argparse,hashlib,json,math,re,unicodedata
from collections import Counter,defaultdict
from pathlib import Path
from statistics import median

from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx,validate_original_keys
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import source_sites_and_guilds
from analyze_Slovenia_rail_spatial_matched_taxon_identity_v142 import original_gps_coordinates
from audit_Eurostat_Slovenia_2015_2020_actual_trains_by_segment_v230 import (
    query_series,unpack_jsonstat,source_segment_count_rows
)

PDF_SOURCE="https://infrastruktura.sz.si/wp-content/uploads/2023/02/NS_2023_1.2_Annex_2B_Stations.pdf"
RE_LINE=re.compile(r"^\s*LINE\s+(\d{2})\b",re.I)
# Source official station label + initial km + final km; no text from
# timetable or unnamed rail track is a substitute for station identity.
RE_STATION=re.compile(r"^\s*([A-Za-zÀ-ÿČčĆćŠšŽžĐđ\s.\-]+?)\s+(\d{1,4}\.\d{3})\s+(\d{1,4}\.\d{3})\s+(\d{2,5})\b",re.U)
N_ORIGINAL_SITES=89
ORIGINAL_SECTIONS=frozenset(("10","20","21","31","32","34","41","42","44","50","64","70","72","80","81"))
NO_CONSTRUCTED_SITE_NETWORK=False

def normalize(s):
    s=s.casefold().replace('–','-').replace('—','-')
    s=unicodedata.normalize('NFKD',s)
    s=''.join(c for c in s if not unicodedata.combining(c))
    s=re.sub(r"\bmeja\b|\bd\.?\s*m\.?\b", "",s)
    s=s.replace('siska','siska').replace('poljcane','poljcane')
    return ' '.join(re.findall(r'[a-z0-9]+',s))

def rows_from_text(lines,initial_code=None):
    current=initial_code;got=[];pending=None
    for raw in lines:
        line=' '.join(raw.split())
        if not line:continue
        m=RE_LINE.match(line)
        if m:
            current=m.group(1)
            pending=None
            continue
        if current is None:continue
        mm=RE_STATION.match(line)
        if mm:
            name=mm.group(1).strip()
            low=float(mm.group(2));high=float(mm.group(3))
            stationlen=int(mm.group(4))
            # The official table reports metres station length, subject to
            # unit consistency and slight source roundings.
            if not(0<=low<high<1000 and stationlen>=20 and
                   abs((high-low)*1000-stationlen)<15):
                pending=None
                continue
            # In source PDF, a multi-line station name may put the first
            # word on its own line before the numeric source row.
            if pending and len(name.split())<=2 and len(pending.split())<=2:
                name=(pending+' '+name).strip()
            pending=None
            got.append({'national_line':current,'name':name,
                        'name_normalized':normalize(name),
                        'station_start_km':low,'station_end_km':high,
                        'station_midpoint_km':(low+high)/2})
            continue
        # Treat isolated candidate prefix text cautiously. Do NOT
        # join source headers / platform data to station names.
        if (len(line)<28 and line[0].isalpha() and
            not any(char.isdigit() for char in line) and
            not re.search(r'\b(electro|central|main|line|equipment|station|network|useful|track|platform)\b',line,re.I)
            and 'see data at line' not in line.lower()):
            pending=line
        elif line.startswith(('1 ','2 ','3 ','4 ','5 ')):
            # Track details need not clear a wrapped name preamble.
            pass
        else:
            pending=None
    return got,current

def extract_authoritative_stations(pdfpath):
    raw=Path(pdfpath).read_bytes()
    if not raw.startswith(b'%PDF-') or len(raw)<50000:
        raise ValueError('Original 2023 national operator PDF invalid/truncated')
    import pdfplumber
    stations=[];preview=[];page_qc=[]
    current=None
    with pdfplumber.open(pdfpath) as pdf:
        if len(pdf.pages)!=16:raise ValueError('Source original annex PDF page count changed')
        for n,page in enumerate(pdf.pages):
            if n<2:continue
            text=page.extract_text(x_tolerance=1.5,y_tolerance=3) or ""
            if len(text)<200:raise ValueError('Official PDF text not extractable, no OCR guessing')
            these,current=rows_from_text(text.splitlines(),current)
            for x in these:x['source_pdf_page_zero_index']=n
            stations.extend(these)
            page_qc.append({'page_zero_based':n,'found_source_station_km_records':len(these),
                            'source_original_heading_count':len(RE_LINE.findall(text))})
            if n in (2,3,4):preview.append({'page':n,'first_lines':text.splitlines()[:18],
                                      'station_labels_first':[(q['national_line'],q['name'],q['station_midpoint_km'])
                                                              for q in these[:10]]})
    # Failing to parse major stations must halt, rather than silently report
    # zero mapped botanical sites as evidence of no exposure.
    if len(stations)<75:raise ValueError('Only %d parsed station positions; source parser too incomplete' % len(stations))
    src_lines={q['national_line'] for q in stations}
    if len(src_lines & ORIGINAL_SECTIONS)<10:raise ValueError('Fewer than ten actual survey railway lines present in official PDF')
    unique=defaultdict(list)
    for s in stations:unique[(s['national_line'],s['name_normalized'])].append(s)
    duplicates={f"{key[0]}:{key[1]}":len(v) for key,v in unique.items() if len(v)>1}
    qa={'sha256_original_pdf_source':hashlib.sha256(raw).hexdigest(),
        'original_pdf_bytes':len(raw),'original_official_pdf_url':PDF_SOURCE,
        'n_original_pages':16,'source_station_records_parsed':len(stations),
        'source_national_line_codes':sorted(src_lines,key=lambda z:int(z)),
        'source_survey_line_coverage':sorted(src_lines & ORIGINAL_SECTIONS,key=int),
        'duplicate_national_line_and_normalized_station_name':duplicates,
        'n_parsed_stations_by_national_line':dict(Counter(q['national_line'] for q in stations)),
        'page_QA':page_qc,
        'first_pages_parser_examples':preview}
    return stations,qa

def eurostat_2020_segments():
    original,source=query_series('rail_tf_ns20_si')
    rows,qa=unpack_jsonstat(original,'2020')
    segment_counts=source_segment_count_rows(rows)
    result=[]
    for (code,label),v in segment_counts.items():
        sep=re.split(r"\s+[-–—]\s+",label)
        if len(sep)!=2:continue
        if 'TRN_GD' not in v:continue
        result.append({'source_network_segment_code':code,
                       'source_original_label':label,
                       'start_station':sep[0].strip(),'end_station':sep[1].strip(),
                       'start_normalized':normalize(sep[0]),
                       'end_normalized':normalize(sep[1]),
                       'actual_2020_freight_trains_per_year':v['TRN_GD'],
                       'actual_2020_passenger_trains_per_year':v.get('TRN_PAS')})
    return result,{'source':source,'segment_value_qa':qa,
                   'n_source_segments_with_parseable_two_endpoints':len(result)}

def segment_route_station_join(stations,segments):
    # Unique names on the SAME national route are required at BOTH ends.
    ix=defaultdict(list)
    for s in stations:ix[(s['national_line'],s['name_normalized'])].append(s)
    resolved=[];reasons=Counter();examples=[]
    for row in segments:
        options=[]
        for line in ORIGINAL_SECTIONS:
            a=ix.get((line,row['start_normalized']),[])
            b=ix.get((line,row['end_normalized']),[])
            if len(a)==1 and len(b)==1 and a[0]['name_normalized']!=b[0]['name_normalized']:
                L=a[0]['station_midpoint_km'];R=b[0]['station_midpoint_km']
                if abs(L-R)>0.1:
                    options.append({'national_line':line,'start_km':L,'end_km':R})
        # Endpoint "Meja" is not a true railway station; do not imply
        # that a border marker equals named station in source table.
        if len(options)==1:
            q={**row,**options[0]}
            resolved.append(q)
        else:
            reasons['none_both_endpoints_unique_same_line' if not options else 'ambiguous_multiple_national_lines']+=1
            if len(examples)<25:examples.append({'source_segment':row['source_original_label'],
              'source_code':row['source_network_segment_code'],'national_line_candidates':len(options)})
    return resolved,{'all_eurostat_2020_source_segments':len(segments),
                     'exact_same_line_unique_two_named_endpoints':len(resolved),
                     'unresolved_reasons':dict(reasons),
                     'examples_25_unresolved':examples,
                     'n_national_lines_matched':len({x['national_line'] for x in resolved}),
                     'matched_by_source_national_line':dict(Counter(x['national_line'] for x in resolved))}

def source_vegetation_sites(original):
    fields,known,meta=source_sites_and_guilds(original)
    withgps,qa=original_gps_coordinates(original,fields)
    if len(withgps)!=N_ORIGINAL_SITES or {x['section'] for x in withgps}!=ORIGINAL_SECTIONS:
        raise ValueError('Genuine publisher author survey site universe drifted')
    if any(x['mid_source_route_chainage_km'] is None for x in withgps):
        raise ValueError('Original real railway route kilometrics not present, cannot map sites')
    return withgps,qa

def assign_segment(sites,matched,margin_km=0.05):
    byline=defaultdict(list)
    for x in matched:byline[x['national_line']].append(x)
    out=[];reason=Counter()
    for site in sites:
        k=site['mid_source_route_chainage_km']
        options=[s for s in byline[site['section']]
                 if min(s['start_km'],s['end_km'])+margin_km <=k <=
                    max(s['start_km'],s['end_km'])-margin_km]
        if len(options)==1:
            x=options[0]
            out.append({'source_site':site['site'],'national_line':site['section'],
                        'eurostat_segment_code':x['source_network_segment_code'],
                        'eurostat_source_segment_name':x['source_original_label'],
                        'historical_2020_goods_trains_per_year':
                              x['actual_2020_freight_trains_per_year'],
                        'historical_2020_passenger_trains_per_year':
                              x['actual_2020_passenger_trains_per_year'],
                        'source_taxa':site['known']})
        else:
            reason['no_two_endpoints_matched_on_line' if len(byline[site['section']])==0 else
                   'no_segment_contains_original_site_midpoint' if not options else
                   'multiple_matching_network_segments']+=1
    cov={'source_89_genuine_sites':len(sites),
         'mapped_to_unique_2020_segment_n':len(out),
         'source_original_lines_with_any_site_mapped':len({x['national_line'] for x in out}),
         'source_line_counts':dict(Counter(x['national_line'] for x in out)),
         'n_unique_2020_traffic_segments_mapped':len({x['eurostat_segment_code'] for x in out}),
         'unmapped_genuine_site_reasons':dict(reason),
         'mapped_genuine_sources_not_zero_imputed':True}
    # No source site coordinates/taxon-by-site assemblage is exported.
    return out,cov

def main():
    p=argparse.ArgumentParser()
    for a in ('official_station_pdf','original_xlsx','source_audit','output'):
        p.add_argument('--'+a.replace('_','-'),required=True)
    a=p.parse_args()
    stations,stationqa=extract_authoritative_stations(a.official_station_pdf)
    if stationqa['sha256_original_pdf_source']!='83f12dfdb91784a653b9a39a8d927a027e7eeb5d8af18ce17b850b5318796425':
        raise ValueError('2023 genuine station source PDF SHA changed: stop')
    original=extract_original_xlsx(a.original_xlsx,a.source_audit)
    sites,gps=source_vegetation_sites(original)
    euro,eq=eurostat_2020_segments()
    matches,mq=segment_route_station_join(stations,euro)
    mapped,coverage=assign_segment(sites,matches)
    result={'version':'v23.1','status':'AUTHENTIC_SOURCE_CROSSWALK_FEASIBILITY_ONLY',
       '2023_rail_operator_station_km_source':stationqa,
       'original_2022_2023_vegetation_geographical_source':gps,
       '2023_official_station_pdf_EXPECTED_sha256':'83f12dfdb91784a653b9a39a8d927a027e7eeb5d8af18ce17b850b5318796425',
       '2020_historical_eurostat_traffic_source':eq,
       'national_route_code_not_assumed_equal_eurostat_segment_code':True,
       '2020_eurostat_named_endpoints_matched_to_authoritative_station_positions':mq,
       'real_original_vegetation_site_to_historical_traffic_coverage':coverage,
       'mapped_goods_2020_summary':{
           'min':min((x['historical_2020_goods_trains_per_year'] for x in mapped),default=None),
           'median':median([x['historical_2020_goods_trains_per_year'] for x in mapped]) if mapped else None,
           'max':max((x['historical_2020_goods_trains_per_year'] for x in mapped),default=None)},
       'censoring_and_causality_cautions':[
           '2020 COVID-affected actual train counts are not 2022/23 field-year freight volumes',
           '2015–2020 code rank stability does not establish matching 2022/2023 traffic rank',
           'Eurostat network segment ≠ Slovenian national route label; join requires exactly matching endpoint station km',
           'unmapped source study sites are unknown, not zero traffic and not absence of railway',
           'station kilometre references and Eurostat timetable segment boundaries can differ and must be verified by geography',
           'no seed arrival counts, inter-station freight OD, local vegetation management or independent regional species pool',
           'do not claim route-level train traffic caused railway traits or station pair identity']}
    # No misdated traffic ecological regression while 2022/23 actual is absent.
    target=Path(a.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf8')
    print('V231_TRUE_OFFICIAL_CHAINAGE_CROSSWALK='+json.dumps(result,ensure_ascii=False,separators=(',',':')),flush=True)
if __name__=='__main__':main()
