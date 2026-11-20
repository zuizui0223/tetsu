"""Original European railway altitude tables and Kashmir rail external transfer."""
import argparse
import json
import re
from pathlib import Path
import numpy as np
import pymupdf
import source_region_rail_elevation_v302 as sources
import validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 as eu
from analyze_region_rail_elevation_transfer_v302 import slope,fit_foreign,predict_score,rail_analysis,ANCHORS


def european_tables(folder):
    italy,_=eu.publisher_italy(folder/'parma.pdf')
    adriatic,_=eu.publisher_hr_ba(folder/'hrba.pdf')
    italy_alt={}
    with pymupdf.open(folder/'parma.pdf') as doc:
        for page in (doc[7],doc[8]):
            words=page.get_text('words')
            codes=[w for w in words if re.fullmatch('[A-Z]{3}',w[4]) and 65<w[0]<90]
            for code in codes:
                numbers=[w for w in words if 312<w[0]<330 and re.fullmatch('[0-9]+',w[4]) and abs(w[1]-code[1])<1]
                if len(numbers)!=1:
                    raise ValueError('Ambiguous Italian published altitude-column cell')
                italy_alt[code[4]]=int(numbers[0][4])
    assert len(italy_alt)==29 and 'CIT' in italy_alt
    with pymupdf.open(folder/'hrba.pdf') as doc:
        text=doc[3].get_text()
    rows=re.findall(r'(?m)^(\d+)\s*\n([^\n]+)\n[^\n]*°[^\n]*\n(\d+)\s*\n(TSS|TWS|MED)\s*$',text)
    assert len(rows)==15 and [int(r[0]) for r in rows]==list(range(1,16))
    adriatic_alt={'source_station_'+r[0]:int(r[2]) for r in rows}
    assert min(italy_alt.values())==32 and max(italy_alt.values())==412
    assert min(adriatic_alt.values())==2 and max(adriatic_alt.values())==348
    return {'IT28':(italy,italy_alt),'HR_BA15':(adriatic,adriatic_alt)}


def run(folder,output):
    source=sources.read_sources(folder)
    studies=european_tables(folder)
    result={}
    choices=[(t,False) for t in sources.src.TREES]+[('height_only',True)]
    for tree_name,height_only in choices:
        tree=None if height_only else sources.tree_context(folder,tree_name)
        for minimum,share,panel in [(5,.25,'primary'),(10,.5,'strict')]:
            values={t:(np.log(h)-sources.FROZEN_3_STUDY_LOG_MEANS[0])/sources.FROZEN_3_STUDY_LOG_SDS[0]
                    for t,h in source['height'].items() if h>0}
            if not height_only:
                values={t:v for t,v in values.items() if t in source['seed'] and t in tree['lookup']}
            training={};support={};slopes={}
            for name,(stations,altitudes) in studies.items():
                x=[];y=[];coverage=[]
                for s in stations:
                    taxa={eu.binomial(t) for t in s['source_names']}-{None}
                    known=taxa&values.keys()
                    if len(known)<minimum or len(known)/len(taxa)<share:
                        continue
                    x.append(altitudes[s['station']]/1000)
                    y.append(np.mean([values[t] for t in sorted(known)]))
                    coverage.append(len(known)/len(taxa))
                x,y=np.array(x),np.array(y)
                assert len(x)>=6 and np.ptp(x)>=.15
                training[name]=[(x,y)]
                deleted=[slope(np.delete(x,i),np.delete(y,i)) for i in range(len(x))]
                slopes[name]=dict(n=len(x),slope=slope(x,y),delete_one_range=[min(deleted),max(deleted)])
                support[name]=dict(n=len(x),min_coverage=min(coverage),median_coverage=float(np.median(coverage)),
                                   altitude_min_max_m=[float(x.min()*1000),float(x.max()*1000)])
            cohorts,audit,_=sources.rail_records(folder,source,tree,minimum,share,height_only)
            external=rail_analysis(cohorts,audit,training)
            # Study-heldout EU check has ONE independent source study for training.
            eu_holdout={}
            for heldout,roads in training.items():
                model=fit_foreign({r:v for r,v in training.items() if r!=heldout},minimum_sources=1)
                x,y=roads[0]
                eu_holdout[heldout]={a:predict_score(x,y,model,a) for a in ANCHORS}
            result[tree_name+'/'+panel]=dict(studies=support,observed_slopes=slopes,
                                            european_study_holdout=eu_holdout,Kashmir=external)
            print(tree_name,panel,slopes,flush=True)
    output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.data,a.output)
