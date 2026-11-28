"""Bounded, outcome-blind inventory of public archived bird tracking data.
Metadata candidates only; no tracking data download or movement inference.
"""
import argparse,collections,csv,hashlib,html,io,json,math,pathlib,re,time
from urllib.parse import urlencode,urlsplit
from urllib.request import Request,urlopen

BASE="https://api.datacite.org/dois"
PREFIX="10.5441/001/1."
AV="https://raw.githubusercontent.com/BrentPease1/alan/34baf9b9241d3c2d55b92287412294cc128b0f72/data/traits/avonet.csv"
SHA="a8f0e0c889a05485f306029433cd45068836324e2bc6cf5b1ed8e34417f3d5fe"
def get(url,cap):
    req=Request(url,headers={"User-Agent":"tetsu-research-catalogue/1.0","Accept":"application/json, text/plain, */*"})
    with urlopen(req,timeout=60) as r:
        b=r.read(cap+1)
    if len(b)>cap:raise ValueError("response exceeds frozen byte limit")
    return b,hashlib.sha256(b).hexdigest()
def normalized(s):
    return " ".join(html.unescape(str(s or "")).lower().replace("_"," ").split())
def avonet():
    raw,sha=get(AV,10_000_000)
    if sha!=SHA:raise ValueError("trait source changed")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(rows)!=11009:raise ValueError("AVONET row count changed")
    out={}
    for r in rows:
        n=normalized(r["Species1"])
        if n in out:raise ValueError("duplicate AVONET name")
        try:v=float(r["Hand-Wing.Index"])
        except (ValueError,TypeError):continue
        if math.isfinite(v) and 0<=v<=100:out[n]=v
    return out
def textfield(a,key,sub):
    return " ".join(normalized(v.get(sub)) for v in a.get(key,[]) if isinstance(v,dict))
def names(text,traits):
    # Exact binomial tokens; explicit title/subject evidence separate from abstract.
    return {m for m in re.findall(r"(?=(\b[a-z][a-z-]+ [a-z][a-z-]+\b))",text) if m in traits}
def license_label(a):
    raw=" ".join(str(v.get(k,"")) for v in a.get("rightsList",[]) for k in ("rights","rightsIdentifier","rightsUri")).lower()
    if re.search(r"\bcc0\b|creativecommons\.org/publicdomain/zero",raw):return "CC0"
    if re.search(r"cc-by-nc|cc-by-nd|cc by-nc|cc by-nd|noncommercial|non-commercial|no derivatives|no-derivatives|creativecommons\.org/licenses/by-(nc|nd)",raw):return "restricted_or_other"
    if re.search(r"\bcc-by(?:-4\.0|-3\.0)?\b|creativecommons\.org/licenses/by/|creative commons attribution [34]",raw):return "CC-BY"
    return "unverified"
def run():
    traits=avonet();url=BASE+"?"+urlencode({"prefix":"10.5441","page[size]":100,"page[cursor]":"1","disable-facets":"true","sort":"name"})
    seen=set();records=[];receipts=[];n_all=None
    while url:
        if len(receipts)>=100:raise ValueError("pagination bound exceeded; not a complete census")
        u=urlsplit(url)
        if u.scheme!="https" or u.hostname!="api.datacite.org" or u.path!="/dois":raise ValueError("unexpected next link")
        raw,sha=get(url,20_000_000);j=json.loads(raw)
        n=int(j["meta"]["total"])
        if n_all is None:n_all=n
        elif n_all!=n:raise ValueError("catalogue count changed during pagination")
        data=j["data"]
        receipts.append({"page":len(receipts)+1,"n":len(data),"sha256":sha})
        for row in data:
            a=row["attributes"];doi=str(a.get("doi",row.get("id",""))).lower()
            if doi in seen:raise ValueError("duplicate DOI")
            seen.add(doi)
            if not doi.startswith(PREFIX) or a.get("types",{}).get("resourceTypeGeneral")!="Dataset":continue
            title=textfield(a,"titles","title")
            subjects=textfield(a,"subjects","subject")
            abstract=textfield(a,"descriptions","description")
            ts=names(title+" "+subjects,traits)
            ab=names(abstract,traits)
            # Mention in text does not verify tracked species or access.
            candidates=ts|ab
            gps=bool(re.search(r"\bgps\b|global positioning system",title+" "+subjects+" "+abstract))
            lic=license_label(a)
            records.append({"doi":doi,"year":a.get("publicationYear"),"gps_mention":gps,"license":lic,
                 "title_subject_match_count":len(ts),"abstract_only_match_count":len(ab-ts),
                 "title_subject_HWI_groups":sorted({"low" if traits[x]<21 else "high" for x in ts}),
                 "title_subject_species":sorted(ts),"abstract_only_species":sorted(ab-ts),
                 "candidate_bird":bool(candidates)})
        print(json.dumps({"catalogue_progress":{"pages":len(receipts),"prefix_records":len(seen),"repository_dataset_candidates":len(records)}}),flush=True)
        url=j.get("links",{}).get("next")
        if len(seen)==n_all:break
        if not url:raise ValueError("pagination ended before total")
        time.sleep(.5)
    if len(seen)!=n_all:raise ValueError("incomplete prefix catalogue")
    qualified=[r for r in records if r["gps_mention"] and r["license"] in ("CC0","CC-BY") and r["title_subject_match_count"]]
    low={x for r in qualified for x in r["title_subject_species"] if traits[x]<21}
    high={x for r in qualified for x in r["title_subject_species"] if traits[x]>=21}
    summary={"prefix_DOIs":len(seen),"repository_dataset_DOIs":len(records),"bird_metadata_candidates":sum(r["candidate_bird"] for r in records),
       "GPS_mention_and_CC0_or_CC_BY_and_title_subject_bird_name":len(qualified),
       "low_HWI_candidate_dataset_DOIs":sum("low" in r["title_subject_HWI_groups"] for r in qualified),
       "high_HWI_candidate_dataset_DOIs":sum("high" in r["title_subject_HWI_groups"] for r in qualified),
       "distinct_low_HWI_candidate_species":len(low),"distinct_high_HWI_candidate_species":len(high),
       "candidate_license_counts":dict(collections.Counter(r["license"] for r in records if r["candidate_bird"]))}
    return {"version":"movebank-datacite-HWI-opportunity-v1","status":"complete prefix census","summary":summary,
        "receipts":receipts,"avonet_sha256":SHA,"qualified_candidates":qualified,
        "limits":["DOI records not independent studies; title/subject and abstract mentions not validated tracked species",
                  "GPS mention does not confirm fix cadence, terrestrial availability, access or suitable steps",
                  "CC0/CC-BY screening only; verify each dataset landing page and author terms before reuse",
                  "no raw coordinates, tracks or ecological effect estimates; not representative of global birds",
                  "repository prefix 10.5441/001/1; later or other prefix records not necessarily included"]}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",default="results/railstatus/movebank_hwi_catalogue.json");p.add_argument("--self-test",action="store_true");a=p.parse_args()
    assert names("bird parus major and corvus corax",{"parus major":14,"corvus corax":30})=={"parus major","corvus corax"}
    assert license_label({"rightsList":[{"rightsIdentifier":"cc-by-nc-4.0"}]})=="restricted_or_other"
    if a.self_test:return
    try:z=run()
    except Exception as e:z={"version":"movebank-datacite-HWI-opportunity-v1","status":"blocked or incomplete","error":type(e).__name__,"reason":str(e)[:300]}
    out=pathlib.Path(a.out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(z,indent=2,ensure_ascii=False,allow_nan=False))
    print(json.dumps({k:v for k,v in z.items() if k!="qualified_candidates"},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
