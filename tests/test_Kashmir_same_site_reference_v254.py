"""Synthetic source-independent mathematical checks for v25.4 (no biological claim)."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import random
import unittest

MOD=Path(__file__).resolve().parents[1]/"scripts/analyze_Kashmir_same_sites_reference_year_v254.py"
spec=importlib.util.spec_from_file_location("v254",MOD)
a=importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


def synthetic():
    plots={}
    for group in a.GROUPS:
        for i in range(6):
            site=f"Between {i}" if group=="between" else f"Station {i}"
            for year in a.YEARS:
                for zone in a.ZONES:
                    seed=(i * 23 + year + (a.ZONES.index(zone)+1)*17 + (group=="between")*31)
                    rng=random.Random(seed)
                    plants=set(rng.sample(range(40),8+i%3))
                    plots[site,year,zone]=plants
    return plots

class SourceMathTests(unittest.TestCase):
    def test_matched_sites_identical_years(self):
        p=synthetic()
        s=a.matching_sites(p,"station",a.ZONES)
        self.assertEqual(len(s),6)
        p[s[0],2017,"Interior Plot 1"]=None
        self.assertEqual(len(a.matching_sites(p,"station",a.ZONES)),5)
        self.assertEqual(len(a.matching_sites(p,"station",("Track","Interior Plot 2"))),6)

    def test_track_cancels_and_bootstrap_uses_paired_sites(self):
        p=synthetic(); sites={g:a.matching_sites(p,g,a.ZONES) for g in a.GROUPS}
        dist=a.pair_distance(p,sites)
        obs=a.calc(dist)
        self.assertTrue(obs["track_cancels_in_reference_sensitivity"])
        # Change track species shared in both reference contrasts; sensitivity stays unchanged.
        for g in a.GROUPS:
            for y in a.YEARS:
                for key in dist[g][y]["Track"]:
                    dist[g][y]["Track"][key]=0.712
        altered=a.calc(dist)
        self.assertAlmostEqual(obs["reference_sensitivity_station_minus_between"],altered["reference_sensitivity_station_minus_between"],places=12)
        b=a.resample_summary(dist,sites,draws=80)
        self.assertEqual(b["bootstrap_requested"],80)
        self.assertGreater(b["by_reference"]["Interior Plot 1"]["draws"],74)

    def test_rarefaction_common_site_year_effort(self):
        p=synthetic(); sites={g:a.matching_sites(p,g,a.ZONES) for g in a.GROUPS}
        rare=a.averaged_pair_distances(p,sites,draws=4)
        self.assertTrue(a.calc(rare)["track_cancels_in_reference_sensitivity"])
        self.assertTrue(all(0<=x<=1 for g in a.GROUPS for y in a.YEARS
                            for z in a.ZONES for x in rare[g][y][z].values()))

    def test_same_site_temporal_turnover_and_site_bootstrap(self):
        p=synthetic(); sites={g:a.matching_sites(p,g,a.ZONES) for g in a.GROUPS}
        out=a.within_site_temporal_turnover(p,sites,draws=5,bootstrap=21)
        self.assertTrue(out["not_traffic_or_heritable_adaptation"])
        for mode in ("raw","rarefied"):
            v=out["results"][mode]["observation"]
            self.assertEqual(set(v["by_reference"]),set(a.REFERENCES))
            self.assertEqual(out["results"][mode]["site_bootstrap_technical"]["Interior Plot 2"]["draws"],21)

    def test_bootstrap_does_not_include_duplicate_site_self_pairs(self):
        weights={"a":4,"b":2,"c":0}
        pairs={("a","b"):0.4,("a","c"):0.9,("b","c"):0.1}
        self.assertAlmostEqual(a.mean_pair(pairs,weights),0.4)
        with self.assertRaises(ValueError):a.mean_pair(pairs,{"a":5,"b":0,"c":0})

    def test_excel_column_index(self):
        self.assertEqual([a.excel_col_index(z) for z in ("A1","Z2","AA1","DY249")],[0,25,26,128])

if __name__=="__main__":unittest.main()
