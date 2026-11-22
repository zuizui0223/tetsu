#!/usr/bin/env python3
"""Inspect original Dryad 2014/2017 Kashmir site metadata prior to temperature analysis.

Fail closed if source does not provide genuine altitude, without substituting
operator station elevations or latitude. Aggregate only, no site names or GPS.
"""
from __future__ import annotations
from pathlib import Path
import argparse,hashlib,json,math
from collections import defaultdict
from analyze_Kashmir_same_sites_reference_year_v254 import parse_xlsx_matrix,DRYAD_SHA256

def audit(path:Path)->dict:
    import pandas as pd
    plots,qa=parse_xlsx_matrix(path)
    frame=pd.read_excel(path)
    if tuple(frame.shape)!=(248,129):
        raise ValueError("Published Kashmir source matrix shape changed")
    headers=[str(x) for x in frame.columns[:6]]
    if str(frame.columns[0])!=str(frame.columns[0]):
        raise ValueError("Corrupt source header")
    if len(set(headers))!=6:
        raise ValueError("Ambiguous original metadata labels")
    originals=frame.copy()
    originals.iloc[:,0]=originals.iloc[:,0].astype(str).str.strip()
    originals.iloc[:,5]=originals.iloc[:,5].astype(str).str.strip()
    originals.iloc[:,4]=originals.iloc[:,4].astype(int)
    if len(originals[originals.iloc[:,4].eq(2014)])!=124:
        raise ValueError("Original 2014 plot row count invalid")
    audit_columns=[]
    for col in range(1,4):
        name=headers[col]
        vals=pd.to_numeric(originals.iloc[:,col],errors="coerce")
        complete=int(vals.notna().sum())
        cells_by_site=defaultdict(set)
        for site,v in zip(originals.iloc[:,0],vals):
            if pd.notna(v):cells_by_site[site].add(float(v))
        one=sum(len(v)==1 for v in cells_by_site.values())
        varying=sum(len(v)>1 for v in cells_by_site.values())
        sampled=[next(iter(v)) for v in cells_by_site.values() if len(v)==1]
        audit_columns.append({
            "original_source_metadata_column_index_zero_based":col,
            "original_source_column_name":name,
            "count_numeric_field_plot_cells":complete,
            "n_original_sites_with_exact_single_numeric_value_across_all_years_and_zones":one,
            "n_original_sites_with_multiple_different_numeric_values_across_years_zones":varying,
            "unique_stable_site_value_minmax":[min(sampled),max(sampled)] if sampled else None,
            "plausible_original_altitude_m":one==31 and all(1000<=x<=2500 for x in sampled)
                                         and max(sampled)-min(sampled)>40,
        })
    plausible=[x for x in audit_columns if x["plausible_original_altitude_m"]]
    return {
      "original_workbook_sha256":DRYAD_SHA256,
      "source_replay_checks":qa,
      "n_original_site_ids":31,
      "original_site_name_frequencies_not_exported":True,
      "original_first_six_metadata_column_headers":headers,
      "source_auxiliary_column_gates":audit_columns,
      "plausible_original_altitude_columns":len(plausible),
      "altitude_identifiable_from_original_Dryad_data":len(plausible)==1 and ("elev" in plausible[0]["original_source_column_name"].lower()
                                      or "alt" in plausible[0]["original_source_column_name"].lower()),
      "do_not_infer_altitude_from_lat_or_station_name":True,
      "raw_names_gps_and_taxon_matrix_exported":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--xlsx",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    r=audit(a.xlsx)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(r,indent=2,ensure_ascii=False))
if __name__=="__main__": main()
