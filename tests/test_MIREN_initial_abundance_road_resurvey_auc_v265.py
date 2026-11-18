"""v26.5 mathematical/source grain tests (all synthetic, not field evidence)."""
from __future__ import annotations
import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_MIREN_original_road_resurvey_initial_abundance_auc_v265 as m


class MIRENRoadAUCTests(unittest.TestCase):
    def test_numeric_and_empty_source(self):
        self.assertEqual(m.numerical("3"),3.)
        self.assertIsNone(m.numerical(""))
        for s in ("0","-1","not-a-number"):
            with self.assertRaises(ValueError):m.numerical(s)

    def test_within_early_plot_auc_with_ties(self):
        early={"Plant A":(3.,2.),"Plant B":(2.,3.),"Plant C":(1.,3.)}
        late={"Plant A":(1.,1.),"Plant X":(3.,1.)}
        x=m.score_auc(early,late,"Abundance")
        y=m.score_auc(early,late,"Cover")
        self.assertEqual(x["persistent_taxa"],1)
        self.assertEqual(x["not_re_recorded_taxa"],2)
        self.assertAlmostEqual(x["auc"],1.)
        self.assertAlmostEqual(y["auc"],0.)
        self.assertTrue(x["identified"])

    def test_undefined_if_abundance_missing_or_no_both_fates(self):
        a={"Plant A":(1.,None),"Plant B":(2.,3.)}
        b={"Plant A":(1.,1.)}
        self.assertFalse(m.score_auc(a,b,"Cover")["identified"])
        self.assertTrue(m.score_auc(a,b,"Abundance")["identified"])
        self.assertFalse(m.score_auc({"Plant A":(2.,2.)},b,"Abundance")["identified"])
        self.assertFalse(m.score_auc(a,{"Plant A":(2.,1.),"Plant B":(1.,2.)},"Abundance")["identified"])

    def test_original_site_road_region_three_level_weighting(self):
        rows=[]
        # Each original road 4 or 8 sites, region 2 roads; therefore
        # this is NOT a plot/taxon weighted 16+ observation mean.
        for reg,eff in [("A",.5),("B",-.5)]:
            for road,n in [("one",4),("two",8)]:
                for site in range(n):
                    edge=.5+eff/2
                    inner=.5-eff/2
                    row={"region":reg,"road":road,"site":str(site),
                         "initial_year":2017,"final_year":2022}
                    for metric in m.METRICS:
                        row[("1",metric)]={"identified":True,"auc":edge}
                        row[("3",metric)]={"identified":True,"auc":inner}
                    rows.append(row)
        summary=m.aggregate(rows,"Abundance","post2012")
        self.assertEqual(summary["source_roads_in_supported_equal_region_summary"],4)
        self.assertEqual(summary["source_sites_in_supported_equal_region_summary"],24)
        self.assertAlmostEqual(summary["mean_original_region_equal_roadside_minus_interior_AUC"],0.)
        self.assertEqual(summary["regions_roadside_AUC_greater"],1)
        self.assertEqual(summary["regions_roadside_AUC_lower"],1)
        same=m.matched_metric_comparison(rows)
        self.assertEqual(same["same_source_sites_for_Abundance_and_Cover"],24)
        self.assertEqual(same["ordinal_cover_on_identical_source_sites"]["source_sites_in_supported_equal_region_summary"],24)

    def test_source_year_cohorts_not_homogenized(self):
        original=[]
        for i,year in enumerate([2007,2012]):
            row={"region":"R","road":"A","site":str(i),"initial_year":year}
            for metric in m.METRICS:
                row[("1",metric)]={"identified":True,"auc":.7}
                row[("3",metric)]={"identified":True,"auc":.6}
            original.append(row)
        self.assertEqual(m.aggregate(original,"Abundance","post2012")[
            "source_original_sites_with_two_plot_AUC_identified_all_road_sizes"],1)
        self.assertEqual(m.aggregate(original,"Abundance","early2007_08")[
            "source_original_sites_with_two_plot_AUC_identified_all_road_sizes"],1)
        self.assertEqual(m.aggregate(original,"Abundance","all_first_last")[
            "source_original_sites_with_two_plot_AUC_identified_all_road_sizes"],2)


if __name__=="__main__":unittest.main()
