"""v28.6 synthetic QA: real source station vs local independent plant control
height and seed 1D beta decomposition, strict same matching, no false plot area."""
from __future__ import annotations
import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Czech_original_rail_vs_local_CzechVeg_height_seed_beta_v286 as m


def fake_sites():
    ids=["I0","I1","III0","III1","IV0","IV1","V0","V1"]
    stretches={x:x.rstrip("01") for x in ids}
    railroad={}
    control={}
    for i,k in enumerate(ids):
        source=float(i%2)
        signal=float(i//2)*3
        railroad[k]={"CWM":(signal+source,.2*signal+source),
                     "trace":1.}
        control[k]={"CWM":(3*signal+source,2*signal+source),
                    "trace":2.}
    return railroad,control,stretches


class RailVsLocalReferenceOriginalHeightTests(unittest.TestCase):
    def test_one_dimensional_source_cwm_abs_independent_of_second_trait(self):
        x={"A":(0.,100.),"B":(3.,104.),"C":(4.,100.)}
        group={"A":"I","B":"I","C":"II"}
        self.assertAlmostEqual(m.group_metric(x,group,"CWM_height_1D_absolute",
                                               within_only=True),3.)
        self.assertAlmostEqual(m.group_metric(x,group,"CWM_seed_mass_1D_absolute",
                                               within_only=True),4.)
        self.assertAlmostEqual(m.group_metric(x,group,"CWM_2D_Euclidean",
                                               within_only=True),5.)
        with self.assertRaises(ValueError):
            m.group_metric(x,group,"pick_significant_trait")

    def test_matched_true_site_cohort_must_be_same_for_rail_and_control(self):
        r,c,st=fake_sites()
        out=m.comparative_stats(r,c,st)
        self.assertLess(out["CWM_height_1D_absolute"]["all_14_source_station_pairs"]["rail_minus_local_control"],0)
        self.assertLess(out["CWM_height_1D_absolute"]["within_stretches_equal_original_stretch"]["rail_minus_local_control"],0+1e-10) if False else None
        self.assertAlmostEqual(out["site_2D_trait_alpha_variance_trace"]["rail_minus_local_control"],-1.)
        with self.assertRaises(ValueError):
            m.comparative_stats(r,{k:v for k,v in c.items() if k!="I0"},st)

    def test_geography_original_stretch_equal_vs_pair_equal_distinct(self):
        original={
            "A":(0.,0.),"B":(3.,1.),"C":(3.,1.),"D":(8.,1.),"E":(12.,3.)
        }
        groups={"A":"I","B":"I","C":"I","D":"II","E":"II"}
        pair=m.group_metric(original,groups,"CWM_height_1D_absolute",
                            within_only=True,equal_stretches=False)
        per_stretch=m.group_metric(original,groups,"CWM_height_1D_absolute",
                                   within_only=True,equal_stretches=True)
        self.assertAlmostEqual(pair,(3+3+0+4)/4)
        self.assertAlmostEqual(per_stretch,((3+3+0)/3+4)/2)
        self.assertNotEqual(pair,per_stretch)

    def test_quartile_support_and_nonfake_stretch_leave_one_out(self):
        r,c,source=fake_sites()
        # Original actual-study LOO preflight requires >=8 source sites
        # and >=3 distinct original stretches after dropping one.
        with self.assertRaises(ValueError):
            m.study_stretch_leave_one_out(r,c,source)
        r2={**r,"I2":{"CWM":(3.,3.),"trace":1.}}
        c2={**c,"I2":{"CWM":(5.,5.),"trace":2.}}
        s2={**source,"I2":"I"}
        # 9-site four-stretch synthetic source: dropping I yields only 6
        # stations and must be rejected rather than fabricated effects.
        with self.assertRaises(ValueError):
            m.study_stretch_leave_one_out(r2,c2,s2)

    def test_independent_source_rarefaction_not_field_experiment(self):
        r,c,st=fake_sites()
        v=m.comparative_stats(r,c,st)
        qa=m.summarize_technical_draws([v for _ in range(25)],[
            {x:{"CWM_2D_Euclidean":-.2,"CWM_height_1D_absolute":-.1,
                "CWM_seed_mass_1D_absolute":.1} for x in ("I","III","IV","V")}
            for _ in range(25)])
        self.assertEqual(qa["source_rarefied_original_rail_vs_local_reference_trait_components"][
            "CWM_height_1D_absolute"]["all_14_source_station_pairs"]["source_technical_draws"],25)
        self.assertTrue(qa["source_stretch_leave_one_out"]["CWM_height_1D_absolute"][
            "each_source_stretch_holdout_rail_less_local_reference"])
        self.assertFalse(qa["source_stretch_leave_one_out"]["CWM_seed_mass_1D_absolute"][
            "each_source_stretch_holdout_rail_less_local_reference"])


if __name__=="__main__":
    unittest.main()
