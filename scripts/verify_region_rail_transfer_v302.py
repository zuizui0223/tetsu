"""Independent aggregate arithmetic and publisher-XLSX railway gradient checks."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from openpyxl import load_workbook
from audit_SI_CZ_original_continuous_trait_coverage_v220 import binomial, read_leda
from audit_Europe_4independent_railstudies_fixed_community_trait_label_null_v285 import FROZEN_3_STUDY_LOG_MEANS, FROZEN_3_STUDY_LOG_SDS


def check_close(actual, reference, label):
    if isinstance(reference, dict):
        assert actual.keys()==reference.keys(),label
        for k in reference:
            check_close(actual[k],reference[k],label+'/'+k)
    elif isinstance(reference,list):
        assert len(actual)==len(reference),label
        for i,(a,b) in enumerate(zip(actual,reference)):
            check_close(a,b,label+'/'+str(i))
    elif isinstance(reference,float):
        assert math.isfinite(actual) and math.isclose(actual,reference,rel_tol=1e-8,abs_tol=1e-10),(label,actual,reference)
    else:
        assert actual==reference,(label,actual,reference)


def verify(data, result):
    checks=0
    for panel, item in result['panels'].items():
        road=item['road']
        for zone,regions in road['regions'].items():
            for region,row in regions.items():
                assert region not in row['training_regions']
                for anchor,score in row['scores'].items():
                    assert abs(score['gain']-(score['baseline_rmse']-score['rmse']))<1e-10
                    assert abs(score['excess_gain']-(score['gain']-score['null_mean_gain']))<1e-10
                    checks+=3
            for anchor,macro in road['macro'][zone].items():
                gains=[r['scores'][anchor]['gain'] for r in regions.values()]
                assert abs(macro['gain']-np.mean(gains))<1e-10
                assert macro['positive_regions']==sum(v>0 for v in gains)
                checks+=2
        for key, row in item['rail'].items():
            for field in ('result','same_sites_both_years'):
                r=row[field]
                if r['status']=='eligible':
                    assert r['n']>=6 and r['elevation_span_m']>=150-1e-9
                    for anchor,s in r['scores'].items():
                        assert s['evaluated_sites']==r['n']-(anchor!='none')
                        assert abs(s['gain']-(s['baseline_rmse']-s['rmse']))<1e-10
                        checks+=2
                else:
                    assert 'scores' not in r
                    checks+=1
    # A fresh openpyxl reconstruction, including complete-plot handling, uses
    # an explicit 2-column least-squares fit rather than the production slope.
    height,_=read_leda(data/'canopy_height.txt','height_m')
    wb=load_workbook(data/'Matrix_plot.xlsx',read_only=True,data_only=True)
    rows=list(wb.active.values);wb.close()
    species=[binomial(v) for v in rows[0][6:]]
    for minimum,share,panel in [(5,.25,'primary_5species_25percent'),(10,.5,'strict_10species_50percent')]:
        outputs=result['panels']['LCVP/'+panel]['rail_height_only_all_LEDA_sensitivity']
        for year in (2014,2017):
            stations={}
            for row in rows[1:]:
                if str(row[0]).startswith('Between') or row[4]!=year or row[5] not in ('Track','Embankment'):
                    continue
                key=' '.join(row[0].split())
                stations.setdefault(key,[]).append(row)
            x,y=[],[]
            for key,parts in sorted(stations.items()):
                assert len(parts)==2
                if any(v is None for row in parts for v in row[6:]):
                    continue
                taxa={species[i] for row in parts for i,v in enumerate(row[6:]) if v>0}-{None}
                known=taxa&height.keys()
                if len(known)<minimum or len(known)/len(taxa)<share:
                    continue
                x.append(parts[0][3]/1000)
                y.append(np.mean([(math.log(height[t])-FROZEN_3_STUDY_LOG_MEANS[0])/FROZEN_3_STUDY_LOG_SDS[0] for t in sorted(known)]))
            observed=outputs[f'station/{year}/Track+Embankment']['result']
            fit=np.linalg.lstsq(np.column_stack([np.ones(len(x)),x]),y,rcond=None)[0]
            assert observed['n']==len(x)
            assert abs(observed['observed_slope']-fit[1])<1e-10
            checks+=2
    return checks

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--data',type=Path,required=True)
    p.add_argument('--result',type=Path,required=True)
    p.add_argument('--reference',type=Path)
    p.add_argument('--height-result',type=Path)
    p.add_argument('--height-reference',type=Path)
    a=p.parse_args()
    result=json.loads(a.result.read_text(encoding='utf-8'))
    n=verify(a.data,result)
    if a.reference:
        check_close(result,json.loads(a.reference.read_text(encoding='utf-8')),'core')
    if a.height_reference:
        check_close(json.loads(a.height_result.read_text(encoding='utf-8')),json.loads(a.height_reference.read_text(encoding='utf-8')),'height')
    print(f'{n} independent arithmetic and raw-XLSX checks PASS; requested frozen comparisons PASS')
