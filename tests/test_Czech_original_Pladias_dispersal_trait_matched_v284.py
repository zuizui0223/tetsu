"""v28.4 synthetic QA for exactly source-guild-preserving Czech railway
study stretch dispersal-proxy interactions; no fabricated field findings."""
from __future__ import annotations
import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Czech_35stations_same_stretch_pladias_dispersal_guild_v284 as v
from audit_Czech_2024_true_35_stations_stretch_geodistance_matched_v282 import (
    candidates
)


def synthetic_sites():
    sites=[]
    guild={}
    for g in v.GUILDS:
        guild[f"{g}_global"]=g
        for stretch in range(5):
            guild[f"{g}_stretch_{stretch}"]=g
        for station in range(7):
            guild[f"{g}_station_{station}"]=g
    for stretch in range(5):
        for i in range(7):
            taxa={t for g in v.GUILDS for t in (
                f"{g}_global",f"{g}_stretch_{stretch}",f"{g}_station_{i}")}
            sites.append({
                "site":f"{stretch}_{i}",
                "section":str(stretch),
                "lat":49+stretch*.06,"lon":14+i*.07,
                "area_ha":float(i+3),
                "known":taxa,"all_original_taxa":taxa|{f"unknown_{i}"},
            })
    return sites,guild


def synthetic_triples(sites):
    geom=[{"source_stretch":x["section"],"lat":x["lat"],"lon":x["lon"],
           "area_ha":x["area_ha"],"taxa":x["all_original_taxa"]}
           for x in sites]
    matched,q=candidates(geom,160,10,None)
    return matched


class GenuineCzechDispersalGuildTests(unittest.TestCase):
    def test_genuine_annotated_taxa_never_assigned_automatically_nonwind(self):
        sites,guild=synthetic_sites()
        m=v.guild_source_matrix(sites,guild)
        self.assertEqual(len(m),4)
        self.assertEqual(len(m["th_wind"]),35)
        for i,s in enumerate(sites):
            self.assertEqual(v.group_set_at(m,"all_544_both_traits",i),s["known"])
            self.assertEqual(
                v.group_set_at(m,"wind_strategy_proxy",i),
                m["th_wind"][i]|m["nonth_wind"][i])
            self.assertEqual(
                v.group_set_at(m,"other_strategy_proxy",i),
                m["th_other"][i]|m["nonth_other"][i])
            self.assertTrue(f"unknown_{int(s['site'].split('_')[1])}" not in
                            v.group_set_at(m,"all_544_both_traits",i))

    def test_identical_source_geographic_triples_and_wind_other_interaction(self):
        s,g=synthetic_sites()
        m=v.guild_source_matrix(s,g)
        triads=synthetic_triples(s)
        self.assertEqual(len(triads),35)
        y=v.matched_source_metrics(m,s,triads)
        self.assertTrue(y[v.MAIN]["identified"])
        self.assertEqual(y[v.MAIN]["source_usable_station_anchors"],35)
        self.assertEqual(y[v.MAIN]["source_distinct_author_stretches"],5)
        for label in v.LABELS:
            self.assertTrue(y[label]["identified"])
        self.assertAlmostEqual(y[v.MAIN]["source_stretch_equal_wind_minus_other_matched_gap"],0.)

    def test_exact_source_4guild_per_site_richness_and_taxon_occupancy_margins(self):
        s,g=synthetic_sites()
        original=v.guild_source_matrix(s,g)
        current={k:[set(x) for x in rows] for k,rows in original.items()}
        for mode in v.NULLS:
            strata=v.source_strata(s,mode)
            n=v.exact_source_margins(current,strata)
            acc=v.draw_source_matrix(current,strata,random.Random(284),125)
            self.assertGreaterEqual(acc,0)
            self.assertEqual(n,v.exact_source_margins(current,strata))
        with self.assertRaises(ValueError):
            v.source_strata(s,"selected_positive_trait_after_seeing_result")

    def test_two_technical_chains_and_frequency_null_are_reproducible(self):
        s,g=synthetic_sites()
        match=synthetic_triples(s)
        first=v.technical_null(s,g,match,v.NULLS[0],seed=4,
                               draws=31,burnin=60,thin=12)
        second=v.technical_null(s,g,match,v.NULLS[0],seed=4,
                                draws=31,burnin=60,thin=12)
        self.assertEqual(first,second)
        self.assertEqual(first["source_technical_correlated_null_draws"],62)
        self.assertEqual(len(first["source_two_independent_seeds_acceptance"]),2)
        self.assertTrue(first["predeclared_full_trait_four_guild_and_wind_other_results"][v.MAIN]["identified"])
        self.assertLessEqual(abs(first["predeclared_full_trait_four_guild_and_wind_other_results"][v.MAIN][
            "observed_minus_frequency_fixed_null"]),2.)

    def test_within_true_source_stretch_preserves_each_original_taxon_occupancy(self):
        s,g=synthetic_sites()
        match=synthetic_triples(s)
        r=v.technical_null(s,g,match,v.NULLS[1],seed=7,
                            draws=31,burnin=50,thin=10)
        self.assertEqual(r["source_technical_correlated_null_draws"],62)
        self.assertTrue(r["site_guild_richness_and_taxon_occurrence_fixed_exact"])
        self.assertEqual(r["predeclared_full_trait_four_guild_and_wind_other_results"][
            "wind_strategy_proxy"]["original_author_stretches"],5)


if __name__=="__main__":
    unittest.main()
