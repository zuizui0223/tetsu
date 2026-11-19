#!/usr/bin/env python3
"""v29.6 fallback: read genuine 3 original V.PhyloMaker2 phylo .rda files via
pure Python rdata (NO R runtime, no apt-get Ubuntu package dependency).

This audits only tree metadata and unambiguous original phylo node structure.
Does NOT fabricate species sequences or ecological traits. The same original
Git blob hashes as R workflow are required. No original tip names exported.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections.abc import Mapping
from pathlib import Path
import numpy as np

COMMIT="7af3fb5152f691af2e4ec9d5e2e467d1b50505e9"
TREES={
    "LCVP":("0d761cb0bf36392cc13759a7ac9026db564217eb",73420,146782),
    "WP":("8d3777762cb32a80e8c6ea097a0456e2b63a94c6",72570,145083),
    "TPL":("d7bc9b69bf3cb11c3c3d090282723b6e73b1de55",74529,149000),
}
def original_git_blob_checksum(path:Path)->str:
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def unwrap_tree(tree)->dict:
    if isinstance(tree,Mapping):
        return dict(tree)
    if hasattr(tree,"_asdict"):
        return dict(tree._asdict())
    if hasattr(tree,"__dict__"):
        return dict(vars(tree))
    raise ValueError("Original author's R phylo S3 named list could not be mapped to a tree")

def inspect_tree_rdata(path:Path, name:str)->tuple[dict,dict]:
    if name not in TREES: raise ValueError("Unregistered original backbone")
    expected,n_tips,n_edges=TREES[name]
    digest=original_git_blob_checksum(path)
    if digest!=expected:raise ValueError("RData original publisher Git blob digest mismatch: "+name)
    import rdata
    try:
        content=rdata.read_rda(path,constructor_dict={"phylo":lambda obj,attrs:obj})
    except Exception as exc:
        raise ValueError("Unable to read original R phylo S3 object") from exc
    object_name="GBOTB.extended."+name
    if object_name not in content:
        raise ValueError("Original RData phylo object name missing; not an RDS or real published tree")
    obj=unwrap_tree(content[object_name])
    if not all(x in obj for x in ("edge","edge.length","tip.label","Nnode")):
        raise ValueError("Original phylo edge/nodes/branch lengths/tip.label incomplete: keys="+str(list(obj)[:20]))
    edges=np.asarray(obj["edge"])
    lengths=np.asarray(obj["edge.length"],dtype=float).reshape(-1)
    tips=np.asarray(obj["tip.label"]).reshape(-1)
    nnodes=np.asarray(obj["Nnode"]).reshape(-1)
    if edges.shape!=(n_edges,2) or len(lengths)!=n_edges or len(tips)!=n_tips:
        raise ValueError("Original V.PhyloMaker2 source tree tip/edge counts changed")
    if nnodes.size!=1 or int(nnodes[0])<1:
        raise ValueError("Invalid tree original Nnode")
    if not np.issubdtype(edges.dtype,np.integer):
        if np.any(edges!=np.floor(edges)):raise ValueError("Non-integer source phylo nodes")
        edges=edges.astype(np.int64)
    if not np.isfinite(lengths).all() or np.any(lengths<0):
        raise ValueError("Original phylogenetic branch lengths are missing or invalid")
    if edges.min()<1 or edges.max()>n_tips+int(nnodes[0]):
        raise ValueError("Original tree node indices outside R phylo 1-based indexing")
    if len(set(str(x) for x in tips))!=n_tips:
        raise ValueError("Original author tree duplicate tip labels")
    # Trust only the original DNA+super-tree era backbone already distributed;
    # do not call every tip direct molecular sequencing, do not graft missing taxa.
    out={
      "name":name,
      "source_git_commit":COMMIT,
      "original_phylo_RData_git_blob_sha":digest,
      "n_original_source_tips":n_tips,
      "n_original_source_edges":n_edges,
      "n_original_internal_nodes":int(nnodes[0]),
      "fraction_positive_original_edges":float(np.mean(lengths>0)),
      "zero_branch_length_source_edges":int(np.sum(lengths==0)),
      "original_tree_tip_sequence_names_not_exported":True,
      "no_plant_taxa_grafted":True,
      "existing_backbone_tip_does_not_imply_own_DNA_sequence":True
    }
    return obj,out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--tree-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    result={}
    for k in TREES:
        path=a.tree_dir/("GBOTB.extended."+k+".rda")
        obj,qa=inspect_tree_rdata(path,k)
        result[k]=qa
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps({"version":"v29.6","renderer":"pure Python rdata",
        "original_tree_schemas_without_tip_identities":result},indent=2),encoding="utf8")
    print(json.dumps({k:{"source_tip_count":v["n_original_source_tips"],
          "edge_count":v["n_original_source_edges"],
          "source_git_blob_SHA_verified":True} for k,v in result.items()},indent=2))
if __name__=="__main__":main()
