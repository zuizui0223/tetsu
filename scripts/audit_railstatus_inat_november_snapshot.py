"""November monthly source lock; no changes to frozen ecological estimator."""
import csv,gzip,io,hashlib,math
from collections import Counter,defaultdict
from pathlib import Path
from botocore.exceptions import ClientError
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
    S3,HashReader,source_client,source_ancestry_ids,OBS_ETAG,OBS_BYTES,TAXA_ETAG,TAXA_BYTES)
from scan_world_railstatus_bird_avonet_guilds import norm,num

FILES=("taxa.csv.gz","observations.csv.gz")
ARCHIVE="metadata/inaturalist-open-data-20261027.tar.gz"
# Fix the November release window before any outcome, not arbitrary latest.
WINDOW=("2026-11-26","2026-11-29")
def receipt(h):
    return {"etag":h["ETag"].strip('"'),"bytes":int(h["ContentLength"]),
            "last_modified":h["LastModified"].isoformat()}
class Snapshot:
    def __init__(self,client):
        self.client=client
        self.files={k:receipt(client.head_object(Bucket=S3,Key=k)) for k in FILES}
        dates=[x["last_modified"][:10] for x in self.files.values()]
        if not all(WINDOW[0]<=d<WINDOW[1] for d in dates):
            raise ValueError("root metadata not in frozen November 2026 release window")
        if any(v["bytes"]<=0 or not v["etag"] for v in self.files.values()):
            raise ValueError("invalid root receipt")
        if self.files["taxa.csv.gz"]["etag"]==TAXA_ETAG or self.files["observations.csv.gz"]["etag"]==OBS_ETAG:
            raise ValueError("mixed old/new monthly root files")
        self.old_archive={"key":ARCHIVE}
        try:
            self.old_archive.update(status="public HEAD available",**receipt(client.head_object(Bucket=S3,Key=ARCHIVE)))
        except ClientError as e:
            self.old_archive.update(status="HEAD unavailable",error_code=str(e.response.get("Error",{}).get("Code","unknown")))
        # Archive is only inspected via HEAD here. Old files not mixed with
        # new ones. Whole archived tar is not downloaded in this workflow.
    def verify(self,key):
        if receipt(self.client.head_object(Bucket=S3,Key=key))!=self.files[key]:
            raise ValueError("monthly source changed during analysis: "+key)
    def open(self,key):
        self.verify(key)
        h=self.files[key]
        r=self.client.get_object(Bucket=S3,Key=key,IfMatch='"'+h["etag"]+'"')
        if int(r["ContentLength"])!=h["bytes"] or r["ETag"].strip('"')!=h["etag"]:
            r["Body"].close()
            raise ValueError("conditional GET mismatch")
        raw=HashReader(r["Body"])
        stream=io.TextIOWrapper(gzip.GzipFile(fileobj=raw),encoding="utf-8",newline="",errors="strict")
        rows=csv.DictReader(stream,delimiter="\t",quoting=csv.QUOTE_NONE)
        return rows,raw,stream

def taxonomy(snapshot):
    rows,raw,stream=snapshot.open("taxa.csv.gz")
    seen=set();birds={};qa=Counter();names=defaultdict(list)
    root=False
    try:
        if not {"taxon_id","ancestry","rank","active","name"}<=set(rows.fieldnames or []):
            raise ValueError("taxon schema changed")
        for r in rows:
            tid=int(r["taxon_id"])
            if tid in seen:raise ValueError("duplicate taxonomy ID")
            seen.add(tid);qa["all_taxon_rows"]+=1
            if tid==3:root=r["name"].strip().lower()=="aves"
            anc=tuple(map(int,source_ancestry_ids(r["ancestry"])))
            if 3 not in anc:continue
            active=r["active"].strip().lower() in ("t","true","1")
            rank=r["rank"].strip().lower()
            if active and rank in ("species","subspecies","variety","form","subvariety"):
                birds[tid]=(rank,anc)
                if rank=="species":names[norm(r["name"])].append(tid)
    finally:stream.close()
    if raw.bytes!=snapshot.files["taxa.csv.gz"]["bytes"] or not root:
        raise ValueError("incomplete taxonomy or Aves root mismatch")
    snapshot.verify("taxa.csv.gz")
    canonical={t:("birds",t) for t,(rank,anc) in birds.items() if rank=="species"}
    qa["active_bird_species"]=len(canonical)
    if qa["active_bird_species"]<8000:raise ValueError("inadequate active bird taxonomy")
    for tid,(rank,anc) in birds.items():
        if rank=="species":continue
        hit=next((canonical[x] for x in reversed(anc) if x in canonical and birds[x][0]=="species"),None)
        if hit is not None:canonical[tid]=hit
    return canonical,names,dict(qa),raw.hash.hexdigest()

def traits(rows,names):
    # Preserve original absolute HWI cutoff 21.0, not updated median.
    out={};qa=Counter()
    for r in rows:
        ids=names.get(norm(r["Species1"]),[])
        if len(ids)!=1:continue
        tid=ids[0]
        qa["unique_exact_matches"]+=1
        h,m,v=map(num,(r["Habitat.Density"],r["Migration"],r["Hand-Wing.Index"]))
        if v is None or not 0<=v<=100:continue
        qa["HWI_known"]+=1
        if h not in (1,2,3) or m not in (1,2,3):continue
        if tid in out:raise ValueError("ambiguous matched trait ID")
        low=int(v<21.)
        out[tid]=(low,f"h{int(h)}_m{int(m)}")
        qa["low_HWI" if low else "high_HWI"]+=1
    if len(out)<8000 or min(qa["low_HWI"],qa["high_HWI"])<3000:
        raise ValueError("insufficient new-snapshot exact trait support")
    return out,dict(qa)
