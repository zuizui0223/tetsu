"""Source-only closure-history feasibility. Never estimate a biological effect."""
import argparse
import collections
import hashlib
import html
import io
import json
import pathlib
import re
import urllib.parse
import urllib.request
import zipfile

JAPAN = "https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N05-2025.html"
OHIO = "https://tims.dot.state.oh.us/ags/rest/services/Assets/Rail_Lines/FeatureServer/0"
FILENAME = "N05-25_GML.zip"

def get(url, max_bytes):
    req = urllib.request.Request(url, headers={"User-Agent": "tetsu-noncommercial-research-source-audit/1.0"})
    with urllib.request.urlopen(req, timeout=45) as r:
        final = urllib.parse.urlsplit(r.geturl())
        if final.hostname != urllib.parse.urlsplit(url).hostname or final.scheme != "https":
            raise ValueError("unexpected redirect host or scheme")
        body = r.read(max_bytes + 1)
        if len(body) > max_bytes:
            raise ValueError("byte budget exceeded")
        return body, {"url": r.geturl(), "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest(),
                      "etag": r.headers.get("ETag"), "last_modified": r.headers.get("Last-Modified")}

def extract_japan_link(page):
    candidates = set()
    for m in re.finditer(r"""(?:https?://|(?:\.\./)+|/)[^"'<>\\\s(),;]*?N05-25_GML\.zip""", html.unescape(page), re.I):
        url = urllib.parse.urljoin(JAPAN, m.group(0))
        parts = urllib.parse.urlsplit(url)
        if parts.scheme == "https" and parts.hostname == "nlftp.mlit.go.jp" and parts.path.endswith("/" + FILENAME):
            candidates.add(url)
    if len(candidates) != 1:
        raise ValueError("official page does not expose one unambiguous exact archive URL")
    return candidates.pop()

def year(x):
    s = str(x or "").strip()
    return int(s) if re.fullmatch(r"\d{4}", s) and 1800 <= int(s) <= 2025 else None

def group_audit(rows):
    groups = collections.defaultdict(list)
    ends = collections.Counter()
    missing = 0
    for r in rows:
        rid = str(r["N05_006"] or "").strip()
        end = str(r["N05_005e"] or "").strip()
        y = year(end)
        ends["present_9999" if end == "9999" else "finite" if y is not None else "unknown"] += 1
        if not rid:
            missing += 1
            continue
        groups[rid].append((end, y, str(r["N05_010"] or "").strip()))
    counts = collections.Counter()
    decades = collections.Counter()
    for records in groups.values():
        if any(a == "9999" for a, _, _ in records):
            counts["group_has_present_record"] += 1
        elif any(y is None for _, y, _ in records):
            counts["group_unknown_end"] += 1
        elif any(link for _, _, link in records):
            counts["group_finite_end_other_relationship_requires_review"] += 1
        else:
            counts["group_finite_end_no_explicit_other_relationship_candidate_only"] += 1
            last = max(y for _, y, _ in records)
            decades[str(((last + 1)//10)*10)] += 1
    return {"row_end_class": dict(ends), "rows_missing_relationship_id": missing,
            "relationship_groups": len(groups), "group_class": dict(counts),
            "candidate_group_end_plus_one_decades_NOT_verified_closures": dict(sorted(decades.items()))}

def japan():
    import shapefile
    page, page_receipt = get(JAPAN, 3_000_000)
    url = extract_japan_link(page.decode("utf-8"))
    data, zip_receipt = get(url, 70_000_000)
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        names = [n for n in z.namelist() if pathlib.PurePosixPath(n).name.lower() == "railroadsection2.dbf"]
        if len(names) != 1:
            raise ValueError("expected one RailroadSection2.dbf")
        if z.getinfo(names[0]).file_size > 150_000_000:
            raise ValueError("DBF byte budget exceeded")
        dbf = z.read(names[0])
        # Read attributes only, no locations or names are exported.
        cpg = names[0][:-4] + ".cpg"
        encoding = "cp932"
        if cpg in z.namelist():
            label = z.read(cpg).decode("ascii").strip().upper()
            if label in {"UTF-8", "UTF8", "65001"}:
                encoding = "utf-8"
            elif label not in {"CP932", "932", "SHIFT_JIS", "SJIS", "SHIFT-JIS"}:
                raise ValueError("unrecognized documented DBF encoding")
        reader = shapefile.Reader(dbf=io.BytesIO(dbf), encoding=encoding)
        fields = [f[0] for f in reader.fields[1:]]
        required = {"N05_004", "N05_005b", "N05_005e", "N05_006", "N05_010"}
        if not required.issubset(fields):
            raise ValueError("documented history schema missing")
        stats = group_audit(r.as_dict() for r in reader.iterRecords())
        return {"status": "attribute audit only", "page": page_receipt, "archive": zip_receipt,
                "dbf_sha256": hashlib.sha256(dbf).hexdigest(), "fields": fields, "encoding": encoding,
                "stats": stats, "limits": "Noncommercial academic use only; no raw data redistribution. End year includes changes, not just closure; group counts are not route/closure counts. Freight-only lines excluded. No geometry or ecology inference."}

def ohio():
    meta, meta_receipt = get(OHIO + "?f=json", 5_000_000)
    m = json.loads(meta)
    fields = {f["name"] for f in m.get("fields", [])}
    oid = m.get("objectIdField")
    if m.get("geometryType") != "esriGeometryPolyline" or not oid or not {oid, "OPERATIONAL_STATUS", "ABANDONMENT_YEAR"}.issubset(fields):
        raise ValueError("documented Ohio geometry or attribute schema missing")
    stats = [{"statisticType": "count", "onStatisticField": oid, "outStatisticFieldName": "n"}]
    query = urllib.parse.urlencode({"where":"1=1","groupByFieldsForStatistics":"OPERATIONAL_STATUS,ABANDONMENT_YEAR",
        "outStatistics":json.dumps(stats),"returnGeometry":"false","f":"json"})
    data, receipt = get(OHIO + "/query?" + query, 5_000_000)
    obj = json.loads(data)
    if "error" in obj or obj.get("exceededTransferLimit") or not isinstance(obj.get("features"), list) or not obj["features"]:
        raise ValueError("grouped query error, truncated or empty")
    counts = collections.Counter()
    decades = collections.Counter()
    for row in obj["features"]:
        r = {k.upper():v for k,v in row["attributes"].items()}
        status = str(r.get("OPERATIONAL_STATUS") or "NULL")
        n = int(r["N"])
        if n < 0:
            raise ValueError("negative row count")
        y = year(r.get("ABANDONMENT_YEAR"))
        counts[(status, "valid_year" if y is not None else "unknown_year")] += n
        if status.lower() == "abandoned" and y is not None:
            decades[str((y//10)*10)] += n
    return {"status": "attribute audit only", "metadata": meta_receipt, "query": receipt,
        "status_year_counts": [{"operational_status": s, "year_class": y, "n_features": n}
            for (s,y),n in sorted(counts.items())],
        "explicit_abandoned_valid_year_feature_decades": dict(sorted(decades.items())),
        "limits": "Features may be split or duplicated; counts are not independent closures. No route lengths, geometry, ecology, or confirmation of final train date. Railbanking/inactive/unknown not recoded as abandoned."}

def safe_source(fn):
    try:
        return fn()
    except Exception as e:
        return {"status":"blocked", "error_type":type(e).__name__, "reason":str(e)[:300], "no_retry":True}

def self_test():
    assert extract_japan_link('<a href="../data/N05/N05-2025/N05-25_GML.zip">x</a>').endswith(FILENAME)
    assert year("9999") is None and year("999") is None and year("2020") == 2020
    rows = [
        {"N05_006":"a","N05_005e":"1999","N05_010":""},
        {"N05_006":"a","N05_005e":"9999","N05_010":""},
        {"N05_006":"b","N05_005e":"2009","N05_010":""},
        {"N05_006":"c","N05_005e":"2010","N05_010":"b"},
        {"N05_006":"d","N05_005e":"999","N05_010":""}]
    a = group_audit(rows)
    assert a["group_class"]["group_has_present_record"] == 1
    assert a["candidate_group_end_plus_one_decades_NOT_verified_closures"] == {"2010":1}
    assert a["group_class"]["group_finite_end_other_relationship_requires_review"] == 1
    assert a["group_class"]["group_unknown_end"] == 1

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", default="results/railstatus/dated_closure_intake.json")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    self_test()
    if not args.self_test:
        result = {"version":"dated-rail-closure-source-intake-v1","japan":safe_source(japan),"ohio":safe_source(ohio),
                  "decision":"Source feasibility only. No habitat or species fitting; do not treat blocked or unknown as zero."}
        out = pathlib.Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        print(json.dumps(result, ensure_ascii=False))
