"""Exact undirected electrical effective resistance on physical source networks.

This is mathematically the *network pairwise* calculation implemented by
Circuitscape. It is not a claim that Julia Circuitscape itself was run.
Resistors are edge-length-equivalent (1 km track = 1 resistance unit),
not biologically calibrated propagule permeability. Parallel edges are
combined as CONDUCTANCES, unlike Dijkstra's shortest-edge representation.
Terminal-preserving leaf pruning and series reduction are exactly
resistance-preserving for all designated focal connectors.

SciPy sparse grounded-Laplacian solve. Never compute a dense full-graph
pseudoinverse. Disconnected terminal pairs are missing, not zero.
"""
from __future__ import annotations

from collections import defaultdict, deque
import math


def add_conductance(graph, u, v, conductance):
    if not isinstance(u,str) or not u or not isinstance(v,str) or not v or u==v:
        raise ValueError("Invalid physical connector IDs")
    if (isinstance(conductance,bool) or not isinstance(conductance,(int,float))
            or not math.isfinite(conductance) or conductance<=0):
        raise ValueError("Conductance must be finite positive")
    current=graph.setdefault(u,{}).get(v,0.)
    new=current+float(conductance)
    if not math.isfinite(new):
        raise ValueError("Nonfinite parallel conductance")
    graph[u][v]=new
    graph.setdefault(v,{})[u]=new


def components(graph):
    """Connected groups of source connector IDs, no geometric line crossing."""
    output=[]
    seen=set()
    for root in graph:
        if root in seen:continue
        seen.add(root)
        stack=[root];nodes=set()
        while stack:
            node=stack.pop();nodes.add(node)
            for other in graph[node]:
                if other not in seen:
                    seen.add(other);stack.append(other)
        output.append(nodes)
    return output


def simplify_resistance_graph(graph,terminals):
    """Exact two-terminal equivalent on all terminal pairs.

    Remove unprotected dangling branches, then contract degree-2
    intermediate nodes. The series resistance 1/g1+1/g2 is in parallel
    with any pre-existing link between the neighboring nodes.
    """
    terminal=set(terminals)
    if not terminal.issubset(graph):
        raise ValueError("Source terminal missing from graph")
    g={u:dict(nbrs) for u,nbrs in graph.items()}
    for u,neighbors in g.items():
        for v,c in neighbors.items():
            if (u==v or v not in g or g[v].get(u)!=c or
                not isinstance(c,(int,float)) or c<=0 or not math.isfinite(c)):
                raise ValueError("Non-symmetric or invalid source conductance")
    initial=len(g)
    leaves=deque(u for u in g if len(g[u])<=1 and u not in terminal)
    leaf_count=0
    while leaves:
        u=leaves.popleft()
        if u not in g or u in terminal or len(g[u])>1:continue
        nbrs=list(g[u])
        if nbrs:
            v=nbrs[0]
            del g[v][u]
            if len(g[v])<=1 and v not in terminal:
                leaves.append(v)
        del g[u];leaf_count+=1
    series=deque(u for u in g if len(g[u])==2 and u not in terminal)
    series_count=0
    while series:
        u=series.popleft()
        if u not in g or u in terminal or len(g[u])!=2:continue
        (a,ca),(b,cb)=list(g[u].items())
        if a==b:raise ValueError("Invalid collapsed multi-edge representation")
        equivalent=1./(1./ca+1./cb)
        del g[a][u];del g[b][u];del g[u]
        add_conductance(g,a,b,equivalent)
        series_count+=1
        for k in (a,b):
            if len(g[k])==2 and k not in terminal:series.append(k)
    if not terminal.issubset(g):
        raise AssertionError("A designated focal terminal was pruned")
    return g,{
        "n_before":initial,"n_after":len(g),
        "n_pruned_dead_end_nodes":leaf_count,
        "n_exact_series_eliminated":series_count,
        "terminal_equivalence_preserved":True,
    }


def pairwise_effective_resistance(graph,terminals,max_reduced_nodes=35000):
    """Sparse exact pairwise resistance, including only terminal-relevant components.

    Returns dict keyed by sorted (terminalA,terminalB) with values in input
    resistor units; disconnected pairs have None. Note numerical solution
    is finite-precision; residual/edge positivity is audited.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import splu
    import numpy as np

    focal=sorted(set(terminals))
    if len(focal)<2:raise ValueError("Need >=2 distinct physical focal nodes")
    if len(focal)!=len(list(terminals)):
        raise ValueError("Duplicate physical focal connectors")
    if not set(focal).issubset(graph):raise ValueError("Missing source focal connector")
    allpairs={(a,b):None for i,a in enumerate(focal) for b in focal[i+1:]}
    details={"n_anchors":len(focal),"n_usable_components":0,
             "n_exact_grounded_solves":0,"n_reduced_nodes_max":0,
             "n_computed_connected_pairs":0,"n_disconnected_pairs":len(allpairs),
             "individual_components":[]}
    for comp in components(graph):
        active=sorted(comp&set(focal))
        if len(active)<2:continue
        sub={node:{nbr:val for nbr,val in graph[node].items() if nbr in comp} for node in comp}
        reduced,info=simplify_resistance_graph(sub,active)
        n=len(reduced)
        if n>max_reduced_nodes:
            raise RuntimeError(f"Exact sparse circuit solver memory gate: reduced graph {n} > {max_reduced_nodes}")
        nodes=sorted(reduced)
        index={node:i for i,node in enumerate(nodes)}
        diag=np.zeros(n,dtype=float)
        row=[];col=[];values=[]
        edges=0
        for u,adj in reduced.items():
            i=index[u]
            for v,c in adj.items():
                if u>=v:continue
                j=index[v]
                diag[i]+=c;diag[j]+=c
                row.extend([i,j]);col.extend([j,i]);values.extend([-c,-c])
                edges+=1
        row.extend(range(n));col.extend(range(n));values.extend(diag.tolist())
        L=coo_matrix((values,(row,col)),shape=(n,n)).tocsc()
        ground=index[active[0]]
        keep=[i for i in range(n) if i!=ground]
        if not keep:raise ValueError("No Laplacian grounded source vertices")
        A=L[keep,:][:,keep].tocsc()
        factor=splu(A)
        lookup={node:j for j,node in enumerate(nodes)}
        reduced_index={old:idx for idx,old in enumerate(keep)}
        voltages={active[0]:np.zeros(n,dtype=float)}
        for focal_node in active[1:]:
            b=np.zeros(n-1,dtype=float)
            b[reduced_index[lookup[focal_node]]]=1.
            solved=factor.solve(b)
            error=np.linalg.norm(A@solved-b,ord=np.inf)
            if not math.isfinite(error) or error>1e-5:
                raise ValueError(f"Grounded Laplacian residual too large: {error}")
            v=np.zeros(n,dtype=float)
            v[keep]=solved
            voltages[focal_node]=v
            details["n_exact_grounded_solves"]+=1
        for i,a in enumerate(active):
            for b in active[i+1:]:
                ia,ib=lookup[a],lookup[b]
                va,vb=voltages[a],voltages[b]
                resistance=(va[ia]-va[ib])-(vb[ia]-vb[ib])
                if resistance< -1e-6 or not math.isfinite(resistance):
                    raise ValueError("Unphysical negative/nonfinite circuit source resistance")
                allpairs[(a,b)]=max(0.,float(resistance))
                details["n_computed_connected_pairs"]+=1
        details["n_usable_components"]+=1
        details["n_reduced_nodes_max"]=max(details["n_reduced_nodes_max"],n)
        details["individual_components"].append({
            "n_original_nodes":info["n_before"],
            "n_reduced_nodes":info["n_after"],
            "n_reduced_edges":edges,
            "n_terminals":len(active),
            "n_pruned_dead_end_nodes":info["n_pruned_dead_end_nodes"],
            "n_eliminated_series_nodes":info["n_exact_series_eliminated"],
        })
    details["n_disconnected_pairs"]=sum(x is None for x in allpairs.values())
    return allpairs,details


def resistor_edges(graph):
    """Circuitscape three-column undirected network resistor input, no mirror duplicates."""
    result=[]
    for u in sorted(graph):
        for v,c in sorted(graph[u].items()):
            if u<v:
                if not c>0:raise ValueError("Bad circuit edge conductance")
                result.append((u,v,1./c))
    return result


def write_circuitscape_network(graph,terminals,outdir,stem="rail_length_ohms"):
    """Optional local adapter for Julia Circuitscape network pairwise mode.

    Exports original network (not a 50km raster), integer connector IDs,
    pairwise focal node IDs and a minimal config. Raw source mapping files
    MUST be handled under Overture/OSM licensing and are NOT CI artifacts.
    """
    from pathlib import Path
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    nodes=sorted(graph);mapping={node:i+1 for i,node in enumerate(nodes)}
    focal=list(terminals)
    if len(set(focal))<2 or any(x not in mapping for x in focal):
        raise ValueError("At least two valid focal source connectors required")
    network=out/(stem+"_network.txt")
    focal_file=out/(stem+"_focal_nodes.txt")
    target=out/(stem+"_result")
    network.write_text("".join(
        f"{mapping[a]}\t{mapping[b]}\t{r:.12g}\n"
        for a,b,r in resistor_edges(graph)
    ),encoding="utf-8")
    focal_file.write_text("".join(f"{mapping[x]}\n" for x in dict.fromkeys(focal)),encoding="utf-8")
    cfg=out/(stem+".ini")
    cfg.write_text(
        "[Circuitscape mode]\ndata_type = network\nscenario = pairwise\n"
        "[Habitat raster or graph]\nhabitat_file = "+str(network.resolve())+
        "\nhabitat_map_is_resistances = True\n"
        "[Focal nodes]\npoint_file = "+str(focal_file.resolve())+
        "\n[Output options]\noutput_file = "+str(target.resolve())+
        "\nwrite_cur_maps = False\nwrite_volt_maps = False\n",
        encoding="utf-8",
    )
    return {"config":str(cfg),"network":str(network),"focal":str(focal_file),
            "n_original_nodes":len(nodes),"n_undirected_edges":len(resistor_edges(graph)),
            "n_focal_nodes":len(set(focal)),
            "julia_circuitscape_not_automatically_executed":True}
