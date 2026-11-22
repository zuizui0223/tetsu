"""Behavioral tests for fixed matched-country donor ranking."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
try:
    import analyze_matched_foreign_rail_road_v314 as m
except ModuleNotFoundError as e:
    if e.name != "analyze_matched_foreign_rail_road_v314": raise
    m=None


def example(n=30):
    names=[f"Taxon species{i}" for i in range(n)]
    f={t:(i+1)/(n+1) for i,t in enumerate(names)}
    r={t:1-f[t] for t in names}
    y=f.copy()
    return names,f,r,y


class MatchedDonorRanking(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(m,"v314 implementation not yet present")

    def test_habitat_exchange_flips_only_signed_effect(self):
        taxa,f,r,y=example()
        a=m.evaluate(f,r,y,taxa)
        b=m.evaluate(r,f,y,taxa)
        self.assertTrue(a["estimable"])
        self.assertAlmostEqual(a["delta_rho"],2.)
        self.assertAlmostEqual(b["delta_rho"],-2.)
        self.assertEqual(a["n_evaluated"],b["n_evaluated"])

    def test_recipient_outcomes_cannot_change_fixed_donors_or_pool(self):
        rail={"CZ":{"Alpha alba":.5,"Beta beta":.2},"IT":{"Alpha alba":.1,"Gamma gamma":.6},
              "IN":{"Alpha alba":1.}}
        road={"CZ":{"Alpha alba":.1,"Gamma gamma":.2},"IT":{"Beta beta":.3}}
        a=m.donor_scores(rail,road,{})
        rail["IN"]={"Unknown taxon":.99}
        b=m.donor_scores(rail,road,{})
        self.assertEqual(a,b)
        self.assertAlmostEqual(a["rail"]["Alpha alba"],.3)
        self.assertAlmostEqual(a["road"]["Alpha alba"],.05)
        self.assertEqual(a["shared"],["Alpha alba","Beta beta","Gamma gamma"])
        reverse=m.donor_scores({"CZ":road["CZ"],"IT":road["IT"]},
                               {"CZ":rail["CZ"],"IT":rail["IT"]},{})
        self.assertEqual(a["shared"],reverse["shared"])

    def test_unknown_frequency_is_excluded_but_source_nonrecord_is_zero(self):
        taxa,f,r,y=example()
        del y[taxa[0]]
        y[taxa[1]]=None
        a=m.evaluate(f,r,y,taxa,unavailable={taxa[2]})
        self.assertEqual(a["n_evaluated"],28)
        self.assertEqual(a["n_recorded"],27)
        self.assertEqual(a["excluded_unavailable"],2)

    def test_low_support_and_constant_scores_are_not_zero_effects(self):
        taxa,f,r,y=example()
        a=m.evaluate(f,r,{t:0. for t in taxa},taxa)
        self.assertFalse(a["estimable"])
        self.assertIsNone(a["delta_rho"])
        b=m.evaluate({t:.5 for t in taxa},r,y,taxa)
        self.assertFalse(b["estimable"])
        self.assertIn("constant_rail_score",b["support_failures"])

    def test_rational_frequency_ties_survive_floating_point_arithmetic(self):
        # The same road-weighted rational can have different floating encodings.
        expected=m.rho([.3,.3,.5,.8],[.1,.2,.3,.4])
        self.assertAlmostEqual(m.rho([.1+.2,.3,.5,.8],[.1,.2,.3,.4]),expected,places=14)
        taxa,f,r,y=example()
        tied={t:(.1+.2 if i%2 else .3) for i,t in enumerate(taxa)}
        result=m.evaluate(tied,r,y,taxa)
        self.assertFalse(result["estimable"])
        self.assertIn("constant_rail_score",result["support_failures"])
        y={t:(.1+.2 if i%2 else .3) for i,t in enumerate(taxa)}
        result=m.evaluate(f,r,y,taxa)
        self.assertEqual(result["positive_frequency_levels"],1)

    def test_recorded_subset_is_separate_from_full_nondetection_estimand(self):
        taxa,f,r,y=example()
        for t in taxa[:15]: y[t]=0.
        a=m.evaluate(f,r,y,taxa)
        b=m.evaluate(f,r,y,taxa,recorded_only=True)
        self.assertEqual(a["n_evaluated"],30)
        self.assertEqual(b["n_evaluated"],15)
        self.assertTrue(a["estimable"])
        self.assertTrue(b["estimable"])
        self.assertAlmostEqual(b["delta_rho"],2.)

    def test_studies_are_equal_weight_and_unsupported_study_is_omitted(self):
        out=m.study_summary([
            {"study":"A","estimable":True,"delta_rho":.4,"n_evaluated":1000},
            {"study":"B","estimable":True,"delta_rho":-.2,"n_evaluated":20},
            {"study":"C","estimable":False,"delta_rho":None,"n_evaluated":50}])
        self.assertAlmostEqual(out["mean_delta_rho"],.1)
        self.assertEqual(out["included_studies"],["A","B"])
        self.assertEqual(out["leave_one_study_mean_range"],[-.2,.4])

    def test_habitat_contrast_uses_same_taxa_after_union_of_uncertain_masks(self):
        taxa,f,r,y=example()
        out=m.habitat_contrast(f,r,y,r,taxa,rail_unavailable={taxa[0]},
                              background_unavailable={taxa[1]})
        self.assertEqual(out["rail"]["n_evaluated"],28)
        self.assertEqual(out["background"]["n_evaluated"],28)
        self.assertAlmostEqual(out["difference_in_advantage"],4.)

    def test_ambiguous_taxon_labels_do_not_become_named_species(self):
        for raw in ["Taraxacum sect. Ruderalia","Rubus fruticosus agg.",
                    "Poa cf. annua","Iris × germanica L.","Festuca rubra s.l.",
                    "Poa sp.","Quercus robur group"]:
            self.assertIsNone(m.normalize_name(raw),raw)
        self.assertEqual(m.normalize_name("Poa annua L."),"Poa annua")
        self.assertEqual(m.normalize_name("Festuca rubra subsp. rubra"),"Festuca rubra")

if __name__=="__main__":
    unittest.main()

