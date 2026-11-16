# v1.1 — Japan official N05 rail history (1950–2025) and independently validated event design

**First successful actual N05 source audit:** [Actions #42438225215](https://github.com/zuizui0223/tetsu/actions/runs/42438225215), executed 2026-11-08 JST. [Source extraction code](../scripts/audit_n05_historical_railway.py), source [MLIT N05 2025](https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N05-2025.html).

**Important legal restriction:** the official N05 time-series railway GIS is **NONCOMMERCIAL-ONLY** because some historical information was derived from JTB publications. The GitHub repo contains only download/audit code and small **factual derived counts and company/route labels**; it does **not** redistribute the original N05 restricted shapefiles, route coordinates, station locations or full historical asset archive. A reader must check their own intended use against the official data-use terms. The present-day N02 and N05 have different license scopes.

## What N05 establishes (new empirical railway infrastructure source)

Unlike the annual N02 current-passenger rail snapshots, N05 records passenger line and station geometry together with route opening year `N05_004`, timeline start `N05_005b`, timeline end `N05_005e` and relation-group ID `N05_006`, over **1950–2025**. Freight-only lines, freight-only stations and marshalling yards are outside its scope. Thus N05 can supply *candidate* treatment chronologies for passenger rail, never measured freight movement.

| Official N05 2025 source layer | History feature rows | Unique GIS relation-group IDs |
|---|---:|---:|
| RailroadSection2.shp | **2,653** | **1,203** |
| Station2.shp | **23,864** | **17,255** |

For line relations, source-audited 2010–2025 opening-year groups = **11**; 2010–2024 non-current groups with a latest ending = **9**. Station analogues are **202 opening-year groups** and **285 terminal-year groups**, but these are NOT independent new/closed stations.

**Interpretation rules:**
- `N05_004` = initial supply opening year for a recorded segment, which can be a small new segment of an old route. It must NOT be called an entirely newly built railway system without external operator verification.
- `N05_005b` = start of the GIS historical *attribute version*, **not necessarily a railway opening**.
- `N05_005e` = the last year of this attribute version. `9999` denotes current, whereas `2022` can indicate a history period ending in 2022 and physical/legal closure in 2023. Attribute/company/name changes can also terminate versions. Never equate N05 `2022` to "2022 last train".
- `N05_006` groups are GIS relation IDs, **not guaranteed one independently operated, geometrically contiguous corridor**. One 2019 candidate ID even has unrelated legacy *Jōbu Railway* and modern Sotetsu labels; it is flagged as requiring independent verification, not accepted blindly as a new route. Several groups of East/West JR Tōkaidō Line describe new line segments on historical routes.
- `1950` is left truncation, and `999` opening unknown. Most source history feature rows have unknown/pre-1950 openings, so do not interpret the 2010–25 event counts as an exhaustive catalogue of all track lengths.
- "Route appears in city bbox plus 5km halo" is not a botanical rail-verge sampling plot.

## Three fixed metropolitan 2018–2025 photo catchments

[Full N05 event geography audit](https://github.com/zuizui0223/tetsu/actions/runs/42438444715), and [small named route candidate audit](https://github.com/zuizui0223/tetsu/actions/runs/42438558961). Fixed Tokyo/Osaka/Sendai bounding boxes are the same as v0.4.

| Pre-fixed pilot rectangle | N05 relation IDs with opening year 2010–25 | Candidate non-current terminal 2010–24 |
|---|---:|---:|
| Tokyo | 3 (2010×1; 2019×2) | 1 (2022×1) |
| Osaka | 1 (2023×1) | 0 |
| Sendai | 1 (2015×1) | 0 |

**Independent operator event verification; N05 is not sufficient on its own:**
- Keisei's *Narita Sky Access* commenced **2010-07-17** ([JRTT](https://www.jrtt.go.jp/construction/achievement/naritaskyaccess.html)), but existing sections were reused and the opened corridor is not uniformly newly built along its full length.
- Sendai Tōzai subway opened **2015-12-06** ([JRTT](https://www.jrtt.go.jp/construction/achievement/sendai-subway.html)); much is underground and cannot be modeled as 14km of exposed trackside grass subject to train slipstream.
- Sotetsu–JR through line opened **2019-11-30** ([Sotetsu](https://www.sotetsu.co.jp/life/history/)); the source 2019 line group mixed earlier unrelated route labels, so independently verify any exact new surface track sections rather than use the GIS relation ID as a treatment.
- West Japan JR Tōkaidō branch was undergrounded **2023-02-13** and the Umekita underground station opened **2023-03-18** ([JR West](https://www.westjr.co.jp/railroad/project/project2/)); this is not a new open-air vegetation strip, though it may alter development/urban land use.
- Tokyo Ueno Zoo monorail **last operated 2019-10-31**, was suspended from **2019-11-01**, and legally closed **2023-12-27** ([Tokyo Metro Transit](https://www.kotsu.metro.tokyo.jp/about/ayumi/others.html)). N05 terminal attribute year 2022 must never be used as the last *operating* year or the year trains generated gusts.
- New at-grade Utsunomiya/Haga **LRT began service 2023-08-26** ([operator](https://www.miyarail.co.jp/info/1143)), confirmed as a source N05 route segment **14.52 km** long. It is the most defensible recent open-air passenger-operation event candidate in the audited routes. But **construction began June 2018** and train trials in 2022, long before regular service ([operator/city development history](https://u-movenext.net/information/164/)); pre-2023 does NOT mean disturbance-free habitat. A properly identified comparison must account for construction separately from operations.

## What follows biologically — strict source and causal gates

Treat N05 as a source to discover **specific passenger-service treatment dates**, then construct the actual ecological analysis only after:
1. independent operator verification of physical at-grade rail and last running date, construction/maintenance timing, land-cover change and discontinuities (avoid treating subway/monorail or simple transfer as ordinary surface rail);
2. candidate route buffer and nearby road-only/old rail active controls in the same climate and urban context, geocoded to individual **plant sampling plots**, not merely photo coordinates near an N05 polyline;
3. pre/post **replicated complete vegetation surveys**, or independently repeatable satellite NDVI/green cover at fixed acquisition season with clouds/land cover controlled; species composition or phenology still needs species-specific field observations;
4. use source/photo effort, observer identity, independent site spatial blocks and placebo years. "First photographed after opening" is NOT first biological arrival.
5. for dispersal, require actual train cargo/passenger mechanism evidence and distinguish freight from passenger; N05 **explicitly excludes freight-only**. For adaptation/plasticity, require independently verified actual mowing/herbicide/wind exposure and repeated fixed-plot plants.

N05 makes a **historical opportunity contrast possible**, but does **not** by itself establish environmental adaptation, new seed colonization, flowering shift or train-induced seed transport.
