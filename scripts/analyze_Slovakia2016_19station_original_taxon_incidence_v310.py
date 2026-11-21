#!/usr/bin/env python3
"""v310 Slovak 2016 original author TRUE railway 19 stations × 318 taxa.

First recover author's printed 318 separate numbered station-incidence clauses
from the accessible 9-page original PDF (NO OCR). Every result is blocked
unless five original station richness values and nine author-declared
all-station species are exactly reproduced from the ORIGINAL integer codes.
No true 2012–14 train counts: this is original railway community structure,
not independent freight-traffic replication.

Original original paper PDF stays transient in CI; aggregate-only output does
not republish author station coordinates or species-by-station incidence table.
"""
from __future__ import annotations
import argparse,hashlib,json,re,math,random,unicodedata
from collections import Counter
from pathlib import Path
import numpy as np
import fitz
from audit_Slovakia2016_original_19station_318taxa_pdf_preflight_v310 import preflight,SOURCE_URL,clean_text

KNOWN_STATION_RICHNESS={1:48,6:146,12:134,15:171,19:202}
PUBLISHER_UNIVERSAL={
    "Achillea millefolium agg.","Arrhenatherum elatius","Artemisia vulgaris",
    "Chenopodium album agg.","Medicago lupulina","Pastinaca sativa",
    "Taraxacum sect. Ruderalia","Trifolium repens","Viola arvensis"
}
EXPECTED_SOURCE_TAXA=318
EXPECTED_ORIGINAL_STATIONS=19
DRAWS_PER_CHAIN=199
SEEDS=(310901,310902)
BURN=18000
THIN=650

def clean_original_page(t):
    # Only source text for published page41–43, after per-page headers.
    t=clean_text(t)
    t=re.sub(r"(?m)^.*Bull\.\s*Slov\.\s*Bot\.\s*Spolo[čc]n\..*$"," ",t)
    t=re.sub(r"(?m)^.*J\.\s*M[áa]jekov[áa]\s*&\s*J\.\s*Lim[áa]nek.*$"," ",t)
    t=re.sub(r"(?m)^\s*(?:3[7-9]|4[0-5])\s*$"," ",t)
    # The original's line-ending hyphenation is broken across visual lines:
    # 'B. hor-\\ndeaceus' is a single taxon 'B. hordeaceus'.
    t=re.sub(r"(?<=[A-Za-zÀ-ž])-\s*\n\s*(?=[A-Za-zÀ-ž])","",t)
    return t

def source_list_clauses(pdf:Path)->tuple[list[str],dict]:
    d=fitz.open(str(pdf))
    if len(d)!=9:raise ValueError("2016 source printed paper must have 9 pages")
    selected=[clean_original_page(d[i].get_text("text",sort=True)) for i in (4,5,6)]
    d.close()
    raw="\n".join(selected)
    if "Acer platanoides:" not in raw or "Poďakovanie" not in raw:
        raise ValueError("Cannot locate original source first and last taxon list boundaries")
    body="Acer platanoides:"+raw.split("Acer platanoides:",1)[1].split("Poďakovanie",1)[0]
    body=re.sub(r"\s+"," ",body).strip()
    clauses=[p.strip() for p in body.split(";") if p.strip()]
    # Last source author entry 'Zea mays: 17.' is period-terminated, not semicolon.
    return clauses,{
        "original_taxon_list_page_numbers_in_author_paper":[41,42,43],
        "n_semicolon_delimited_species_incidence_clauses":len(clauses),
        "does_source_last_clause_end_with":clauses[-1][-12:] if clauses else None,
        "author_Slovakia2016_primary_taxon_list_reconstructed_with_no_OCR":True,
    }

def expand_station_codes(code:str)->set[int]:
    val=code.strip().rstrip(".").strip()
    if not val or re.search(r"[^0-9\s,–\-]",val):
        raise ValueError("Original taxon station incidence codes contain non-number, non-range content: "+repr(val[:50]))
    result=set()
    for term in val.split(","):
        term=term.strip()
        if not term:raise ValueError("Original station code contains empty token")
        if "–" in term or "-" in term:
            pieces=re.split("[–-]",term)
            if len(pieces)!=2 or not all(p.strip().isdigit() for p in pieces):
                raise ValueError("Malformed author original station code range")
            lo,hi=map(lambda z:int(z.strip()),pieces)
            if not (1<=lo<=hi<=EXPECTED_ORIGINAL_STATIONS):
                raise ValueError("Original station incidence interval outside 1..19")
            obs=set(range(lo,hi+1))
        else:
            if not term.isdigit():raise ValueError("Malformed original integer station number")
            n=int(term)
            if not (1<=n<=EXPECTED_ORIGINAL_STATIONS):raise ValueError("Source station ID outside range 1..19")
            obs={n}
        if result&obs:raise ValueError("Duplicate source station ID inside one taxon list")
        result|=obs
    return result

def parse_clauses(clauses:list[str])->tuple[dict[str,set[int]],dict]:
    if len(clauses)!=EXPECTED_SOURCE_TAXA:
        raise ValueError(f"Only {len(clauses)} source taxon clauses; original publisher reported 318; STOP")
    obs={}
    genus=None
    abbreviated=0
    ambiguous=[]
    for i,clause in enumerate(clauses):
        if clause.count(":")!=1:
            raise ValueError(f"Taxon clause {i+1} has {clause.count(':')} colons, unexpected non-source data")
        name,rest=(x.strip() for x in clause.split(":",1))
        if re.search(r"\d",name):
            raise ValueError("Source scientific name is contaminated with location numbers: "+name)
        m=re.match(r"^([A-Z])\.\s+(.+)$",name)
        if m:
            prefix,epithet=m.groups()
            if genus is None or genus[0]!=prefix:
                ambiguous.append({"original_index":i+1,"abbreviation":name})
                raise ValueError("Do not invent abbreviated source genus "+repr(name))
            name=genus+" "+epithet
            abbreviated+=1
        else:
            candidate=name.split()[0]
            if not re.match(r"^[A-Z][a-zà-ž]+$",candidate):
                raise ValueError("Unexpected nonbotanical full genus source taxon label: "+repr(name))
            genus=candidate
        if name in obs:
            raise ValueError("Duplicate original source scientific name after normalization: "+repr(name))
        obs[name]=expand_station_codes(rest)
    return obs,{
        "n_original_full_taxon_labels":len(obs),
        "n_original_genus_abbreviation_expansions_from_prior_source_full_genus":abbreviated,
        "n_ambiguous_source_species_labels_silently_corrected":0,
        "n_original_species_with_empty_source_station_list":sum(not x for x in obs.values()),
        "n_original_site_taxon_positive_cells":sum(map(len,obs.values()))
    }

def validate_original_source(obs:dict[str,set[int]])->dict:
    taxa=list(obs)
    if len(taxa)!=EXPECTED_SOURCE_TAXA:raise ValueError("Original taxon count wrong")
    site={k:set() for k in range(1,20)}
    for taxon,sites in obs.items():
        for k in sites:site[k].add(taxon)
    n={k:len(t) for k,t in site.items()}
    fails={k:{"extracted":n[k],"original_published":v} for k,v in KNOWN_STATION_RICHNESS.items() if n[k]!=v}
    if fails:
        raise ValueError("Original 2016 author independent five-station richness QA FAIL: "+json.dumps(fails))
    universal={t for t,s in obs.items() if len(s)==19}
    if universal!=PUBLISHER_UNIVERSAL:
        raise ValueError("Original 9 author universal source taxa misparsed: "+
                         str(sorted(universal^PUBLISHER_UNIVERSAL)[:15]))
    if min(n.values())<10 or len(site)!=19:raise ValueError("Unreasonable original source station subset")
    return {
       "n_pub_source_total_taxa":len(obs),
       "n_author_numbered_rail_stations":len(site),
       "n_original_station_taxon_presence":sum(n.values()),
       "five_original_pub_station_richness_values_reproduced_exactly":dict(KNOWN_STATION_RICHNESS),
       "nine_original_universal_taxa_strict_source_list_verified":True,
       "n_true_author_taxa_at_all_19_stations":len(universal),
       "source_station_richness_min":min(n.values()),
       "source_station_richness_max":max(n.values()),
       "source_station_richness_mean":float(np.mean(list(n.values()))),
       "source_taxon_occupancy_fraction_singleton":float(np.mean([len(s)==1 for s in obs.values()])),
       "no_unobserved_plant_absence_imputed":True
    }

def station_sets(obs):
    return {k:{t for t,x in obs.items() if k in x} for k in range(1,20)}

def jaccard(a,b):
    if not a|b:return None
    return len(a&b)/len(a|b)

def beta_slope(sets,dist):
    pairs=[(i,j) for i in range(1,20) for j in range(i+1,20)]
    vals=np.array([jaccard(sets[i],sets[j]) for i,j in pairs])
    d=np.array([dist(i,j) for i,j in pairs],float)
    if not np.isfinite(vals).all():raise ValueError("Source similarity undefined")
    X=np.stack([np.ones(len(d)),d],axis=1)
    return float(np.linalg.lstsq(X,vals,rcond=None)[0][1]), float(np.mean(vals))

def swap_states(sets,rng,nsteps):
    accepted=0
    for _ in range(nsteps):
        i,j=rng.sample(range(1,20),2)
        a,b=sets[i]-sets[j],sets[j]-sets[i]
        if not a or not b:continue
        x=rng.choice(sorted(a));y=rng.choice(sorted(b))
        sets[i].remove(x);sets[i].add(y);sets[j].remove(y);sets[j].add(x)
        accepted+=1
    return accepted

def fixed_margin_null(obs,draws=DRAWS_PER_CHAIN):
    if draws!=199:raise ValueError("Source design frozen at two×199 autocorrelated technical draws")
    original=station_sets(obs)
    oc=Counter(x for ss in original.values() for x in ss)
    counts={k:len(v) for k,v in original.items()}
    original_slope,original_mean=beta_slope(original,lambda i,j:j-i)
    technical=[]
    chain=[]
    for ichain,seed in enumerate(SEEDS):
        sets={k:set(v) for k,v in original.items()}
        rng=random.Random(seed)
        accepted=swap_states(sets,rng,BURN)
        y=[]
        for _ in range(draws):
            accepted+=swap_states(sets,rng,THIN)
            if Counter(x for xx in sets.values() for x in xx)!=oc or any(len(sets[k])!=v for k,v in counts.items()):
                raise ValueError("Site richness and original taxon occupancy NOT conserved")
            v,_=beta_slope(sets,lambda i,j:j-i)
            technical.append(v);y.append(v)
        chain.append({"chain":ichain+1,
                      "fraction_accepted_source_checkerboard_proposals":accepted/(BURN+draws*THIN),
                      "mean_null_station_order_slope":float(np.mean(y))})
    null=np.array(technical)
    return {
      "response":"Jaccard of original botanical taxon identities at 19 source railway stations",
      "main_predictor":"absolute original station index gap along a single line, NOT physical rail travel OD",
      "original_station_ordinal_gap_Jaccard_slope":original_slope,
      "original_mean_station_pair_Jaccard":original_mean,
      "n_dependent_original_station_pairs_NOT_independent":19*18//2,
      "source_margin_null_mean_ordinal_gap_slope":float(null.mean()),
      "original_minus_null_source_margin_slope":float(original_slope-null.mean()),
      "technical_null_range_2p5_97p5_NOT_field_CI":list(map(float,np.quantile(null,[.025,.975]))),
      "technical_lower_rank_NOT_independent_ECOLOGICAL_p":float((1+sum(null<=original_slope))/(1+len(null))),
      "technical_upper_rank_NOT_independent_ECOLOGICAL_p":float((1+sum(null>=original_slope))/(1+len(null))),
      "original_site_count_and_taxon_incidence_margins_held_exact":True,
      "source_replicated_chain_diagnostics":chain,
      "n_source_technical_autocorrelated_null_draws_NOT_independent_field_studies":len(technical),
    }

def run(path:Path)->dict:
    raw=preflight(path)
    if raw["source_incidence_list_first_page_index"] is None:
        # The preflight may fail to locate the LIST HEADER when the original
        # author's glyph text breaks mid-word, but direct taxon boundary is
        # independently verified below. Do not call missing header a success.
        pass
    clauses,caq=source_list_clauses(path)
    obs,tqa=parse_clauses(clauses)
    sourceqa=validate_original_source(obs)
    patterns=fixed_margin_null(obs)
    return {
      "version":"v310", "source":"publisher-authored Slovakia 2016 Poprad–Plaveč original paper",
      "source_original_pdf_sha256":raw["pdf_sha256"],
      "source_pdf_url":SOURCE_URL,
      "source_sha_authenticated_by_original_live_fetch":True,
      "original_source_taxon_list_parser_qa":{**caq,**tqa},
      "author_published_numbers_validated_from_original_rows":sourceqa,
      "source_rail_single_line_taxon_identity_spatial_structure":patterns,
      "n_real_original_freight_or_passenger_train_exposure_values":0,
      "source_raw_species_by_station_not_republished":True,
      "notes":["one author source study, 19 sites and 171 dependent pairs; not independent field n=171",
               "source railway station-wide vegetation, not Slovenia BS/BT 100-m paired microhabitats",
               "original 2012-14 train traffic, field-level management and local nonrail vegetation not measured",
               "original source taxa status cannot be inferred from paper-wide alien counts",
               "2016 Slovakia and 2025 Slovakia may share physical sampling stations, not automatically independent publications"]
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--pdf",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    data=run(args.pdf)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps(data,ensure_ascii=False,indent=2))
if __name__=="__main__":main()
