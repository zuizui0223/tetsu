#!/usr/bin/env python3
"""Link validated fixed quadrat photo ROI to documented cutting dates.

No imputation, no biological claim, no p-values. Restricted plots preserved.
"""
from __future__ import annotations
import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.mowing_recovery import load_csv,summarize_mowing_roi


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--roi",required=True,help="CSV from extract_photo_roi.py")
    p.add_argument("--events",required=True,help="independently documented mowing events CSV")
    p.add_argument("--output-dir",default="results/mowing")
    args=p.parse_args()
    measurements, rejected=summarize_mowing_roi(load_csv(args.roi),load_csv(args.events))
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    if measurements:
        keys=list(dict.fromkeys(k for row in measurements for k in row))
        with (out/"plot_event_roi.csv").open("w",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=keys)
            writer.writeheader();writer.writerows(measurements)
    report={
        "description":"descriptive fixed-plot photo-ROI response aligned with independently verified mowing",
        "n_plot_cut_events":len(measurements),
        "n_valid_pre_recovery_pairs":sum(x["eligibility"]=="pass" for x in measurements),
        "n_restricted":sum(x["eligibility"]!="pass" for x in measurements),
        "photo_rejections":rejected,
        "interpretation_limits":[
            "Visible vegetation pixels / fixed frame is not plant height or absolute biomass.",
            "Flower pixels / conspecific vegetation pixels is not reproductive biomass or seed production.",
            "Before/after reflects season, observer and other events as well as cuts; not causal without matched untreated controls.",
            "Do not infer train slipstream exposure or heritable railway adaptation from ROI."
        ]
    }
    (out/"audit.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
