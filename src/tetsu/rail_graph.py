"""Pure railway-network topology primitives for macroecological cell-pair exposure.

Construct paths from Overture segment connectors, NOT intersecting drawn
lines. The graph describes MAPPED physical connectivity only: Overture can
include disused rail, tunnels, switches and incomplete cross-window edges.
No actual service, passenger/freight tonnage or propagule flux is observed.

Deliberately no external dependencies, so inference pipelines can audit
network logic independently of remote API/data fetching.
"""
from __future__ import annotations

from collections import defaultdict
from heapq import heappush, heappop
import math


def validate_segment(raw, allowed_classes):
    if not isinstance(raw,dict) or raw.get("subtype")!="rail":
        raise ValueError("Rail segment required")
    if raw.get("class") not in allowed_classes:
        return None
    uid=str(raw.get("id") or "").strip()
    if not uid:
        raise ValueError("Missing segment identity")
    length=raw.get("length_m")
    if isinstance(length,bool) or not isinstance(length,(int,float)) or not math.isfinite(length) or length<=0:
        raise ValueError("Rail segment length must be positive finite meters")
    raw_connectors=raw.get("connectors")
    if not isinstance(raw_connectors,list) or len(raw_connectors)<2:
        raise ValueError("At least two explicit Overture connectors required")
    parsed=[]
    for c in raw_connectors:
        if not isinstance(c,dict):
            raise ValueError("Connector must be a mapping")
        cid=str(c.get("connector_id") or "").strip()
        at=c.get("at")
        if not cid or isinstance(at,bool) or not isinstance(at,(int,float)) or not math.isfinite(at) or at<0 or at>1:
            raise ValueError("Connector id and linear reference must be valid")
        parsed.append((float(at),cid))
    parsed.sort(key=lambda item:(item[0],item[1]))
    if len(set(x[1] for x in parsed))!=len(parsed):
        raise ValueError("Duplicate connector id within one segment; ambiguous linear reference")
    if parsed[0][0]!=0. or parsed[-1][0]!=1.:
        raise ValueError("Overture segments must have endpoint connectors at 0 and 1")
    if any(b[0]<=a[0] for a,b in zip(parsed,parsed[1:])):
        raise ValueError("Two connectors may not share an 'at' position")
    return {"id":uid,"length_m":float(length),"connectors":parsed}


def build_rail_connector_graph(segments,allowed_classes=("standard_gauge","narrow_gauge","light_rail")):
    """Build weighted undirected connector graph. Use source identifiers.

    Each source segment may include interior connectors; edge weights are
    segment-length * difference in linear reference, a proxy for path
    distance. For precise geometry distances, preprocessing must measure
    each sub-line rather than use overall length fractions.
    """
    if not allowed_classes or len(set(allowed_classes))!=len(allowed_classes):
        raise ValueError("Invalid allowed rail-class set")
    adj=defaultdict(dict)
    seen=set()
    stats={"n_source_segments":0,"n_excluded_classes":0,"n_eligible_segments":0,
           "n_connector_links":0,"source_rail_length_m":0.}
    for raw in segments:
        stats["n_source_segments"]+=1
        item=validate_segment(raw,allowed_classes)
        if item is None:
            stats["n_excluded_classes"]+=1
            continue
        sid=item["id"]
        if sid in seen:
            raise ValueError("Duplicate Overture source segment ID")
        seen.add(sid)
        stats["n_eligible_segments"]+=1
        stats["source_rail_length_m"]+=item["length_m"]
        for (a,u),(b,v) in zip(item["connectors"],item["connectors"][1:]):
            length=(b-a)*item["length_m"]
            if length<=0:
                raise ValueError("Zero-length connector edge")
            if v in adj[u]:
                # parallel geometry may have different sources; retaining the
                # shortest edge reflects the shortest path only.
                length=min(length,adj[u][v])
            adj[u][v]=length
            adj[v][u]=length
            stats["n_connector_links"]+=1
    stats["n_distinct_connectors"]=len(adj)
    stats["source_rail_length_m"]=round(stats["source_rail_length_m"],3)
    stats["incomplete_at_spatial_window_boundary"]=True
    stats["freight_or_actual_train_service_unidentified"]=True
    return {k:dict(v) for k,v in adj.items()},stats


def shortest_connector_path_m(graph,start,end,max_distance_m=None):
    """Dijkstra over explicit physical connectors; None means unmapped/disconnected."""
    if start not in graph or end not in graph:
        return None
    if start==end:
        return 0.
    if max_distance_m is not None and (not math.isfinite(max_distance_m) or max_distance_m<0):
        raise ValueError("Invalid maximum path distance")
    frontier=[(0.,start)]
    best={start:0.}
    while frontier:
        dist,node=heappop(frontier)
        if dist>best[node]:
            continue
        if max_distance_m is not None and dist>max_distance_m:
            return None
        if node==end:
            return round(dist,6)
        for nxt,w in graph[node].items():
            if not math.isfinite(w) or w<=0:
                raise ValueError("Bad graph edge weight")
            proposed=dist+w
            if proposed<best.get(nxt,float("inf")):
                best[nxt]=proposed
                heappush(frontier,(proposed,nxt))
    return None


def graph_components(graph):
    """Components are in-source-window connectivity, not full national lines."""
    seen=set()
    sizes=[]
    for node in graph:
        if node in seen:
            continue
        todo=[node]
        seen.add(node)
        count=0
        while todo:
            curr=todo.pop()
            count+=1
            for nxt in graph[curr]:
                if nxt not in seen:
                    seen.add(nxt)
                    todo.append(nxt)
        sizes.append(count)
    sizes.sort(reverse=True)
    return sizes
