#!/usr/bin/env python3
"""v384b: summarize whole-continent influence from existing original v382 result."""
import argparse,json,statistics
from pathlib import Path
def summarize(obj):
    if obj["version"]!="v382-original-world-plant-richness-and-Jaccard-tree-road-falsification":
        raise ValueError("Not original v382")
    out=[]
    for p in obj["panels"]:
        r={"width_m":p["width_m"],"min_land_fraction":p["min_land_fraction"],"n_cells":p["n_matched"],"n_pairs":p["n_pairs"]}
        for label in ("richness_prediction","jaccard_prediction"):
            a=p[label]["base_road__to__both"]
            folds=a["source_independent_continent_holdouts"]
            effects=[f["delta_MSE_positive_means_rail_added_prediction_helped"] for f in folds]
            if len({f["heldout_continent"] for f in folds})!=len(folds):raise ValueError("Repeated fold")
            mean=statistics.mean(effects) if effects else None
            if mean is not None and abs(mean-a["mean_holdout_continent_delta_MSE"])>1e-12:
                raise ValueError("Original reported mean not reproduced")
            drops=[statistics.mean(effects[:i]+effects[i+1:]) for i in range(len(effects))] if len(effects)>1 else []
            base=statistics.mean([f["baseline_test_blockweighted_MSE"] for f in folds]) if folds else None
            r[label]={"n_folds":len(folds),"positive_folds":sum(v>0 for v in effects),
                      "mean_delta_MSE":mean,"relative_prediction_improvement_pct":100*mean/base if base else None,
                      "leave_one_continent_score_out_min":min(drops) if drops else None,
                      "leave_one_continent_score_out_max":max(drops) if drops else None,
                      "fold_delta":{f["heldout_continent"]:d for f,d in zip(folds,effects)}}
        out.append(r)
    return {"version":"v384b-v382-original-fold-influence-summary","panels":out,
        "limit":"Leave-one-fold-score-out sensitivity is NOT a refit without that continent; folds share training observations, not independent statistical replicates."}
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    out=summarize(json.loads(a.source.read_text()));a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False));print(json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
