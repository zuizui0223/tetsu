"""Exploratory height-only geographic support extension, declared after core run began."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from source_region_rail_elevation_v302 import read_sources, src, FROZEN_3_STUDY_LOG_MEANS, FROZEN_3_STUDY_LOG_SDS
from audit_MIREN_original_native_elevation_source_v299 import original_paired_site_elevation
from analyze_region_rail_elevation_transfer_v302 import region_holdout, ANCHORS, ZONES, fit_foreign, rail_analysis
from source_region_rail_elevation_v302 import rail_records


def run(folder,output):
    source=read_sources(folder)
    height={t:(math.log(v)-FROZEN_3_STUDY_LOG_MEANS[0])/FROZEN_3_STUDY_LOG_SDS[0]
            for t,v in source['height'].items() if v>0 and math.isfinite(v)}
    result={}
    for minimum,share,name in [(5,.25,'primary'),(10,.5,'strict')]:
        regions,qa=src.valid_region_panels(source['flora'],height,minimum,share,2007)
        arrays={z:{} for z in ZONES}
        support={}
        for region,roads in sorted(regions.items()):
            chosen=[]
            for road,rows in roads:
                paired=[(r,original_paired_site_elevation(r,region,road,source['elevations'])) for r in rows]
                pairs=[(r,e) for r,e in paired if e is not None]
                if len(pairs)<6:
                    continue
                x=np.array([e/1000 for r,e in pairs])
                if np.ptp(x)<.15-1e-12:
                    continue
                chosen.append((x,{z:np.array([np.mean([height[t] for t in sorted(r[z])]) for r,e in pairs]) for z in ZONES}))
            original_n=sum(len(r) for _,r in roads)
            n=sum(len(x) for x,y in chosen)
            passed=len(chosen)>=2 and n>=.7*original_n
            support[region]=dict(original_sites=original_n,retained_sites=n,roads=len(chosen),eligible=passed)
            if passed:
                for z in ZONES:
                    arrays[z][region]=[(x,y[z]) for x,y in chosen]
        folds={z:region_holdout(a) for z,a in arrays.items()}
        macro={}
        for z in ZONES:
            macro[z]={}
            for anchor in ANCHORS:
                scores=[r['scores'][anchor]['gain'] for r in folds[z].values()]
                macro[z][anchor]=dict(gain=float(np.mean(scores)),positive=sum(v>0 for v in scores),n=len(scores))
        # Europe-trained slope alone tests transfer across a continent boundary.
        europe={'AUT','CHE','ITA','NOR'}
        europe_training={r:v for r,v in arrays['roadside'].items() if r in europe}
        from analyze_region_rail_elevation_transfer_v302 import predict_score,average_scores
        external={r:{a:average_scores([predict_score(x,y,fit_foreign(europe_training),a) for x,y in roads])
                     for a in ANCHORS} for r,roads in arrays['roadside'].items() if r not in europe}
        rc,ra,_=rail_records(folder,source,None,minimum,share,height_only=True)
        rail=rail_analysis(rc,ra,europe_training)
        result[name]=dict(source_support=qa,retained_support=support,regions=folds,macro=macro,
                          Europe_trained_nonEurope_road=external,Europe_trained_Kashmir_rail=rail)
    output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--data',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    run(a.data,a.output)
