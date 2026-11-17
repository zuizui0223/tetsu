import hashlib
import re
from collections import Counter
import pymupdf
import pandas as pd
import validate_Europe_3studies_82stations_original_LEDA_2D_trait_beta_v245 as eu
from analyze_Slovenia_source_rail_functional_taxonomic_identity_v141 import source_sites_and_guilds
from select_two_original_rail_20_species_for_inat_v164 import original_slovenian_EPPO_codebook_binomial_mapping

def names(raw):
    return {eu.binomial(t) for t in raw if t} - {None}

def frequencies(sites):
    count=Counter(t for site in sites for t in site)
    return {t:n/len(sites) for t,n in count.items()}

def load_rail(data):
    studies={}
    for tag, result in [('CZ39',eu.publisher_czech(data/'cz1.xlsx',data/'cz2.xlsx')),
                        ('IT28',eu.publisher_italy(data/'parma.pdf')),('HR_BA15',eu.publisher_hr_ba(data/'hrba.pdf'))]:
        studies[tag]=frequencies([names(x['source_names']) for x in result[0]])
    si=next((data/'si').rglob('*.xlsx'))
    assert hashlib.sha256(si.read_bytes()).hexdigest()=='26be7ca9513dd54136d4ec499d8413c1cb288045f87fe3cf8814731f2edac99f'
    from analyze_Slovenia_rail_microhabitats_v103 import extract_original_xlsx
    original=extract_original_xlsx(si,next((data/'si').rglob('*.json')))
    sites,_,_=source_sites_and_guilds(original)
    code,_=original_slovenian_EPPO_codebook_binomial_mapping(original)
    studies['SI89']=frequencies([names([code.get(t) for t in site['all_original_taxa']]) for site in sites])
    assert hashlib.sha256((data/'sk30_2.pdf').read_bytes()).hexdigest()=='ff8271a396ecda01d92bba9817721577bcdd44fdac1309cc4884240c9e3c9838'
    with pymupdf.open(data/'sk30_2.pdf') as doc:
        txt='\n'.join(p.get_text() for p in doc)
    entries=re.findall(r'(?m)^([^\n]+)\n(\d+)\s*%\s*$',txt)
    assert len(entries)==413
    sk={}
    for raw,freq in entries:
        key=eu.binomial(raw.strip())
        if key: sk.setdefault(key,[]).append(int(freq)/100)
    # Species aggregates with ambiguous duplicate marginals cannot be unioned.
    studies['SK30']={t:v[0] for t,v in sk.items() if len(v)==1}
    k=pd.read_excel(data/'Matrix_plot.xlsx')
    railrich={}
    kashmir_floras={}
    for year,sub in k.groupby('Year of sampling'):
        st=sub[~sub.iloc[:,0].str.startswith('Between')].copy()
        st['site']=st.iloc[:,0].str.strip()
        track=st[st.Location=='Track'].set_index('site').sort_index().iloc[:,6:]
        bank=st[st.Location=='Embankment'].set_index('site').sort_index().iloc[:,6:]
        assert track.index.equals(bank.index) and len(track)==16 and not (track.isna().any().any() or bank.isna().any().any())
        counts=(bank>0).sum(axis=1)-(track>0).sum(axis=1)
        railrich[str(year)]=dict(mean=float(counts.mean()),positive=int((counts>0).sum()),ties=int((counts==0).sum()),n=len(counts))
        for label,frame in [('track',track),('bank',bank)]:
            kashmir_floras[f'{year}_{label}']=frequencies([names([c for c in frame if v[c]>0]) for _,v in frame.iterrows()])
        if year==2014:
            studies['Kashmir16']=frequencies([names([c for c in track if v[c]>0]) for _,v in track.iterrows()])

    return studies, railrich, kashmir_floras
