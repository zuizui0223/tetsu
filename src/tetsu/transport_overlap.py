"""Observation-network overlap controls; NOT transport or adaptation inference.

Rail proximity means distance to MLIT N02 passenger rail centreline.
Road proximity means distance to MLIT N13 road centreline.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict

def page_cluster_inclusion_weight(total: int, chosen_count: int, per_page: int = 200) -> float:
    """Inverse inclusion probability for ONE uniformly drawn whole ID page.

    Every observation in a page is selected with probability 1/n_pages;
    a short final page does not give each observation a larger probability.
    Chosen ID page is a correlated cluster, not n independent samples.
    """
    if not isinstance(total, int) or total < 0:
        raise ValueError("total must be a nonnegative integer")
    if not isinstance(chosen_count, int) or chosen_count < 0 or chosen_count > per_page:
        raise ValueError("invalid chosen page size")
    if not isinstance(per_page, int) or per_page <= 0:
        raise ValueError("per_page must be positive integer")
    if total == 0:
        if chosen_count:
            raise ValueError("nonempty selected page when total=0")
        return 1.0
    if chosen_count == 0:
        raise ValueError("nonempty population cannot select an empty page")
    return float(math.ceil(total / per_page))


def road_rail_status(rail_m:float, road_m:float, road_limit=100.0,
                     rail_near=250.0, rail_far=1000.0)->str:
    for val in (rail_m,road_m):
        if not math.isfinite(val) or val<0:
            raise ValueError("distances must be nonnegative and finite")
    if rail_near>=rail_far or road_limit<=0:
        raise ValueError("invalid predefined distance thresholds")
    if road_m>road_limit:
        return "road_distant"
    if rail_m<=rail_near:
        return "road_matched_rail_near"
    if rail_m>=rail_far:
        return "road_matched_rail_far"
    return "road_matched_rail_middle_excluded"


def weighted_case_control_table(rows, *, rail_near=250., road_limit=100.):
    """Compute two-by-two counts from COMPLETE target + sampled reference rows.

    Each row needs label=target/reference, rail_m, road_m, ref_weight.
    Reference weights correct only the known marginal sampling fraction:
    selection of an iNaturalist result page is NOT random inside that page.
    """
    t=Counter()
    regions=set()
    for r in rows:
        cat=road_rail_status(r["rail_m"],r["road_m"],
                             rail_near=rail_near,road_limit=road_limit)
        if cat not in ("road_matched_rail_near","road_matched_rail_far"):
            continue
        label=r["label"]
        if label not in ("target","reference"):
            raise ValueError(f"unexpected label {label}")
        weight=(float(r["ref_weight"]) if label=="reference" else 1.0)
        if not math.isfinite(weight) or weight<=0:
            raise ValueError("weights must be finite positive")
        t[(label,"near" if cat.endswith("_near") else "far")]+=weight
        if "grid2km" in r:
            regions.add(r["grid2km"])
    a=t[("target","near")]; b=t[("target","far")]
    c=t[("reference","near")]; d=t[("reference","far")]
    # Continuity correction is disclosed; not a p-value or causal estimator.
    ratio=((a+0.5)/(c+0.5))/((b+0.5)/(d+0.5))
    return {
        "target_rail_near":int(a), "target_rail_far":int(b),
        "reference_weighted_rail_near":round(c,2),
        "reference_weighted_rail_far":round(d,2),
        "approx_target_to_reference_ratio_of_ratios":round(ratio,4),
        "n_grid2km_with_valid_comparison_rows":len(regions),
        "n_target_total_valid":int(a+b),
        "caution":"target=full source query, reference=single selected API page per month, weighted only by count; not a representative random survey"
    }


def monthly_mh_ratio(rows, road_limit=100., rail_near=250.):
    """Mantel-Haenszel seasonal contrast, descriptive only."""
    cells=defaultdict(Counter)
    for row in rows:
        status=road_rail_status(row["rail_m"],row["road_m"],
                                road_limit=road_limit,rail_near=rail_near)
        if status not in ("road_matched_rail_near","road_matched_rail_far"):
            continue
        w=float(row["ref_weight"]) if row["label"]=="reference" else 1.0
        cells[int(row["month"])][(row["label"],status.endswith("_near"))]+=w
    num=0.0; den=0.0; represented=0
    for tab in cells.values():
        a=tab[("target",True)]; b=tab[("target",False)]
        c=tab[("reference",True)]; d=tab[("reference",False)]
        n=a+b+c+d
        if n<=0: continue
        if a+b>0 and c+d>0:
            represented+=1
            num+=a*d/n
            den+=b*c/n
    return {"seasonal_mantel_haenszel_ratio":round(num/den,4) if den>0 else None,
            "months_with_both_target_and_reference":represented,
            "note":"month-stratified descriptive sample contrast; not a biological population odds ratio"}
