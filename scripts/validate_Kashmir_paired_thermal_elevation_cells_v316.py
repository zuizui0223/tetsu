#!/usr/bin/env python3
"""Independent numerical arithmetic for Kashmir v316 using openpyxl direct cells.

Does not import v316 code or JAE XML parser. It DOES reuse the published
v315 trait resolver and v314 strict nomenclatural unavailable-key policy, so
it is an arithmetic/source-cell check, not an independent taxonomic authority.
No original locations or species-by-transect matrix leave the runner.
"""
import argparse,hashlib,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np
from openpyxl import load_workbook

from analyze_matched_foreign_rail_road_v314 import normalize_name
from source_thermal_affinity_v315 import load_traits
from source_matched_rail_road_v314 import audit_names,unresolved_keys

SFILE="529f3b97521933018e12623b61f531f012c67241cae2602a6675dbd534192474"
ZONES=("Track","Embankment","Interior Plot 1","Interior Plot 2")

def validate(original:Path,traitdir:Path,derived:dict):
    if hashlib.sha256(original.read_bytes()).hexdigest()!=SFILE:
        raise ValueError("Second direct-cell validator received non-original publisher matrix")
    b=load_workbook(original,read_only=True,data_only=True)
    try:
        w=b.worksheets[0]
        if (w.max_row,w.max_column)!=(249,129):
            raise ValueError("Wrong genuine publisher original shape")
        rows=w.iter_rows(values_only=True)
        headers=next(rows)
        names=headers[6:]
        raw={}
        alt=defaultdict(set)
        for cells in rows:
            site=" ".join(str(cells[0]).split())
            year=int(cells[4])
            zone=" ".join(str(cells[5]).split())
            if zone not in ZONES or year not in (2014,2017):
                raise ValueError("Unexpected original site-year-habitat cell")
            a=float(cells[3]);alt[site].add(a)
            source_abundance=cells[6:]
            raw[(site,year,zone)]=(None if any(x is None for x in source_abundance)
                                 else {names[i] for i,v in enumerate(source_abundance) if v>0})
        if len(raw)!=248 or len(alt)!=31 or any(len(v)!=1 for v in alt.values()):
            raise ValueError("Original 2014/2017 source missing plot or unstable altitude")
        # Independent data walk, but the same predeclared taxon-mask rule.
        seen=set().union(*(v for v in raw.values() if v is not None))
        name_map=audit_names(seen,normalize_name)
        excluded=unresolved_keys(name_map,seen,normalize_name)
        traits=load_traits(traitdir,normalize_name)
        results={}
        for year in (2014,2017):
            results[str(year)]={}
            for indicator,table in (("tichy_numeric",traits["tichy"]),("eive_all",traits["eive"])):
                samples=[]
                source_pairs=0
                for site in sorted(alt):
                    groups=[raw.get((site,year,z)) for z in ZONES]
                    if any(v is None for v in groups):continue
                    source_pairs+=1
                    rail={name_map[s] for s in groups[0]|groups[1]}-{None}-excluded
                    interior={name_map[s] for s in groups[2]|groups[3]}-{None}-excluded
                    score=table["values"]
                    rail_values=[score[s] for s in rail if s in score]
                    inner_values=[score[s] for s in interior if s in score]
                    if min(len(rail_values),len(inner_values))<3:continue
                    paired=float(math.fsum(rail_values)/len(rail_values)-math.fsum(inner_values)/len(inner_values))
                    samples.append((next(iter(alt[site])),int(not site.startswith("Between ")),paired))
                n=len(samples)
                if n<20:raise ValueError("Second original direct-cell check has too few scored paired sites")
                e=np.array([t[0] for t in samples],float)
                st=np.array([t[1] for t in samples],float)
                yy=np.array([t[2] for t in samples],float)
                X=np.column_stack((np.ones(n),(e-e.mean())/100,st))
                independent=np.linalg.lstsq(X,yy,rcond=None)[0]
                old=derived["main_2014" if year==2014 else "year_2017_sensitivity"][indicator]
                if n!=old["model"]["n_transects"] or source_pairs!=old["original_rows_after_require_all_4_source_plots"]:
                    raise ValueError("Independent original-cell taxon support differs from v316")
                if abs(independent[1]-old["model"]["altitude_slope_per_100m"])>1e-10:
                    raise ValueError("Independent original-cell altitude slope differs from v316")
                if abs(yy.mean()-old["model"]["mean_paired_thermal_difference"])>1e-10:
                    raise ValueError("Independent original-cell paired thermal mean differs")
                results[str(year)][indicator]={
                    "n_actual_source_transects":n,
                    "n_four_plot_source_pairs":source_pairs,
                    "independently_computed_altitude_slope_per_100m":float(independent[1]),
                    "independently_computed_pair_mean":float(yy.mean()),
                    "absolute_slope_difference_from_v316":float(abs(independent[1]-old["model"]["altitude_slope_per_100m"])),
                    "independent_original_openpyxl_cell_arithmetic_match":True
                }
        return {"source":"independently reread original XLSX cells with openpyxl",
                "shared_dependencies":"v315 Zenodo trait resolver and v314 scientific-name ambiguity masks",
                "n_source_rows_verified":248,"n_source_physical_transects":31,
                "source_original_matrix_sha256":SFILE,"results":results,
                "site_taxon_and_gps_republished":False}
    finally:b.close()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--traits",required=True,type=Path)
    p.add_argument("--result",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    r=validate(a.source,a.traits,json.loads(a.result.read_text()))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"source_original_field_cells_checked":r["n_source_rows_verified"],
                      "recomputed":r["results"],
                      "all_original_numeric_matches":True},indent=2,ensure_ascii=False))
if __name__=="__main__":main()
