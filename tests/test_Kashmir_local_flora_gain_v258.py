"""v25.8 mathematical and cohort QA on synthetic source-structured flora."""
import itertools
import random
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Kashmir_2014_local_flora_predict_2017_rail_gains_v258 as m


class LocalFloraTemporalGainTests(unittest.TestCase):
    def test_temporal_gain_loss_retained_bookkeeping(self):
        sites=["Site A"]
        x={
            ("Site A",2014,"Track"):{1,2,3},
            ("Site A",2017,"Track"):{2,3,4},
        }
        out=m.time_records(x,sites,"Track")["Site A"]
        self.assertEqual(out["gained"],{4})
        self.assertEqual(out["lost"],{1})
        self.assertEqual(out["retained"],{2,3})
        self.assertEqual(out["before_n"],3)
        self.assertEqual(out["after_n"],3)

    def test_source_pool_only_from_2014_and_not_track(self):
        x={
            ("S",2014,"Track"):{999},
            ("S",2014,"Embankment"):{30,31},
            ("S",2014,"Interior Plot 1"):{1,2},
            ("S",2014,"Interior Plot 2"):{2,3},
            ("S",2017,"Interior Plot 1"):{500},
        }
        self.assertEqual(m.source_pool(x,"S","interiors_2014"),{1,2,3})
        self.assertEqual(m.source_pool(x,"S","embankment_2014"),{30,31})
        self.assertEqual(m.source_pool(x,"S","all_nontrack_2014"),{1,2,3,30,31})
        self.assertNotIn(999,m.source_pool(x,"S","all_nontrack_2014"))
        self.assertNotIn(500,m.source_pool(x,"S","all_nontrack_2014"))
        with self.assertRaises(ValueError):
            m.source_pool(x,"S","illegal")

    def test_exact_leave_same_site_out_expectation(self):
        records={
            "A":{"gained":frozenset({11}),"before_n":1,"after_n":2},
            "B":{"gained":frozenset({12}),"before_n":1,"after_n":2},
            "C":{"gained":frozenset({13}),"before_n":1,"after_n":2},
        }
        pools={
            "A":frozenset({11}),
            "B":frozenset({12}),
            "C":frozenset({13}),
        }
        r=m.observed_and_analytic_expected(records,pools,"gained")
        self.assertEqual(r["observed_site_equal_hit_share"],1.)
        self.assertEqual(r["source_other_site_exact_expected_site_equal"],0.)
        self.assertEqual(r["local_excess_exact"],1.)
        self.assertEqual(r["n_target_taxon_events"],3)
        # All possible 3-site derangements yield 0 when all gain taxa
        # are strictly site-private in 2014 interior samples.
        from analyze_rail_local_site_coupling_v257 import derangement
        for i in range(100):
            mapping=derangement(list(records),random.Random(i))
            self.assertEqual(m.permuted_score(records,pools,"gained",mapping),0.)

    def test_common_regional_species_not_falsely_enriched(self):
        rec={x:{"gained":frozenset({1}),"before_n":1,"after_n":2} for x in ("A","B","C","D")}
        donors={x:frozenset({1}) for x in rec}
        out=m.observed_and_analytic_expected(rec,donors,"gained")
        self.assertEqual(out["observed_site_equal_hit_share"],1.)
        self.assertEqual(out["source_other_site_exact_expected_site_equal"],1.)
        self.assertEqual(out["local_excess_exact"],0.)

    def test_absent_gains_are_unknown_not_zero(self):
        r=m.observed_and_analytic_expected(
            {s:{"gained":frozenset()} for s in ("A","B","C")},
            {s:frozenset({1}) for s in ("A","B","C")},"gained")
        self.assertFalse(r["identified"])
        self.assertNotIn("local_excess_exact",r)

    def test_source_original_site_strata_not_pseudoreplicated(self):
        # 12 synthetic sites, unique local 2014 species; all four zones
        # exist in both years. Original site is the one matching unit.
        plots={}
        sites=[f"Station {i}" for i in range(12)]
        for i,s in enumerate(sites):
            local={100+i}
            for year in m.YEARS:
                for zone in m.ZONES:
                    if zone=="Track":
                        plot={1} if year==2014 else {1,100+i}
                    elif zone=="Embankment":
                        plot={2} if year==2014 else {2,100+i}
                    elif zone=="Interior Plot 1":
                        plot={100+i,200+i}
                    else:
                        plot={100+i,300+i}
                    plots[(s,year,zone)]=plot
        r=m.group_analysis(plots,sites,199,99)
        self.assertEqual(r["original_source_sites"],12)
        a=r["scenario_metrics"]["Track__interiors_2014"]["gain"]
        b=r["scenario_metrics"]["Embankment__interiors_2014"]["gain"]
        self.assertAlmostEqual(a["local_excess_exact"],1.)
        self.assertAlmostEqual(b["local_excess_exact"],1.)
        self.assertAlmostEqual(r["track_vs_embankment_using_identical_interior_2014_source_pool"]["track_minus_embankment_local_excess_above_exact_own_null"],0.)
        self.assertEqual(a["n_target_taxon_events"],12)

    def test_technical_rank_and_math(self):
        r=m.technical_diagnostic(.7,[.1,.2,.4]*50)
        self.assertAlmostEqual(r["upper_tail_rank_plusone_mapping_only_NOT_field_p"],1/151)
        self.assertLess(r["null_mean"],.4)


if __name__=="__main__":
    unittest.main()
