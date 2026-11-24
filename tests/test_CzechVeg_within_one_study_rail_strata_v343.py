"""v343 synthetic within one source study rail yard vs embankment check."""
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_CzechVeg_within_one_study_rail_strata_v343 as m

def source_pair(study,typ,i):
    return {"rail":{"id":f"{study}r{i}","reference":study,"lat":49.0+.001*i,
                    "lon":15.0,"altitude":250},
            "reference":{"id":f"{study}c{i}","reference":study,
                         "lat":49.0+.001*i,"lon":15.0001,"altitude":250},
            "rail_habitat":typ}

class WithinSourceRailCategories(unittest.TestCase):
    def test_only_study_with_both_habitats_selected(self):
        rows=[source_pair("A",m.HABITATS[0],i) for i in range(18)]
        rows += [source_pair("A",m.HABITATS[1],i+18) for i in range(9)]
        rows += [source_pair("B",m.HABITATS[0],28)]
        got,token=m.original_study_with_both(rows)
        self.assertEqual(len(got),27)
        self.assertEqual(len(token),10)
    def test_no_two_railway_habitat_types_does_not_prove_selection(self):
        with self.assertRaises(ValueError):
            m.original_study_with_both([source_pair("A",m.HABITATS[0],1)])
    def test_same_study_unmatched_stratum_effects_additive(self):
        samples=[
          {"habitat":m.HABITATS[0],"effect":.5,"cover_effect":.3,
           "source_rail_species":8,"source_ref_species":10},
          {"habitat":m.HABITATS[1],"effect":-.3,"cover_effect":-.2,
           "source_rail_species":10,"source_ref_species":8}]
        out,delta=m.type_diffs(samples)
        self.assertAlmostEqual(delta,.8)
        self.assertEqual(out[m.HABITATS[0]]["positive_site_pairs"],1)
        self.assertEqual(out[m.HABITATS[1]]["negative_site_pairs"],1)
    def test_mixed_plot_types_nearby_only_one_to_one(self):
        yr=[{"rail":{"id":"r1","lat":49.,"lon":15.,"altitude":250},
             "effect":.5,"habitat":m.HABITATS[0]},
            {"rail":{"id":"r2","lat":49.001,"lon":15.,"altitude":250},
             "effect":.3,"habitat":m.HABITATS[0]},
            {"rail":{"id":"r3","lat":49.0012,"lon":15.,"altitude":251},
             "effect":-.2,"habitat":m.HABITATS[1]}]
        d=m.nearby_rail_type_pairs(yr)
        self.assertEqual(d["n_source_rail_types_nearby"],1)
        self.assertGreater(d["mean_rail_yard_minus_embankment_thermal_effect_for_nearby_source_sites"],0)
    def test_no_nearby_habitat_pair_not_imputed_zero(self):
        yr=[{"rail":{"id":"r1","lat":49.,"lon":15.,"altitude":250},
             "effect":.5,"habitat":m.HABITATS[0]},
            {"rail":{"id":"r3","lat":50.,"lon":16.,"altitude":500},
             "effect":-.2,"habitat":m.HABITATS[1]}]
        d=m.nearby_rail_type_pairs(yr)
        self.assertEqual(d["n_source_rail_types_nearby"],0)
        self.assertIsNone(d["mean_rail_yard_minus_embankment_thermal_effect_for_nearby_source_sites"])
    def test_source_body_not_fake_independent_publication(self):
        self.assertEqual(m.GATE,"G2_same_source_extended")
        self.assertEqual(m.MAIN,(3,5,.5))
    def test_original_study_empty_sensitivity_is_unestimable(self):
        d,delta=m.type_diffs([])
        self.assertIsNone(delta)
        self.assertEqual(d[m.HABITATS[1]]["n_site_vs_nearby_same_study_pairs"],0)

if __name__=="__main__":unittest.main()
