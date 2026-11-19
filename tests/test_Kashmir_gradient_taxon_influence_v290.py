import importlib.util
import pathlib
import unittest

SCRIPT=pathlib.Path(__file__).resolve().parents[1]/"scripts"/"analyze_Kashmir_gradient_taxon_influence_v290.py"
spec=importlib.util.spec_from_file_location("kashmir_v290",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class KashmirGradientTaxonSensitivityTests(unittest.TestCase):
    def setUp(self):
        self.names=[
           "Cynodon dactylon (L.) Pers.",
           "Breea arvensis (L.) Less.",
           "Taraxacum officinale Weber",
           "Lactuca serriola L.",
           "Senecio chrysanthemoides DC.",
           "Lotus corniculatus L.",
           "Poa annua L.",
           "Plantago major L.",
        ]
        self.traits={
           0:"perennial", 1:"perennial", 6:"annual",7:"perennial",
        }
        self.plots={}
        for site in ("Station A","Station B","Between A","Between B"):
            for yr in (2014,2017):
                self.plots[(site,yr,"Track")]={0,1,2,6}
                self.plots[(site,yr,"Embankment")]={1,3,6,7}
                self.plots[(site,yr,"Interior Plot 1")]={3,4,6,7}
                self.plots[(site,yr,"Interior Plot 2")]={3,5,6,7}
        self.plots[("Station B",2017,"Track")]={0,2,6,7}

    def test_main_table_contain_all_physical_reference_pairs(self):
        r=mod.produce(self.plots,self.traits,self.names)
        self.assertEqual(r["source_site_count"],4)
        self.assertEqual(len(r["all_zone_paired_references"]),5)
        self.assertEqual(r["all_zone_paired_references"]["Track__minus__Interior Plot 1"]["all"]["two_year_mean"]["original_sites"],4)

    def test_site_leaveout_uses_one_spatial_site_with_two_years(self):
        r=mod.leave_one_original_site(self.plots,self.traits,len(self.names))
        a=r["all"]["Track__minus__Interior Plot 1"]
        self.assertEqual(a["LOSO_total_folds"],4)
        self.assertEqual(r["between"]["Track__minus__Interior Plot 1"]["LOSO_total_folds"],2)
        self.assertEqual(r["station"]["Track__minus__Interior Plot 1"]["LOSO_total_folds"],2)
        self.assertLessEqual(a["LOSO_negative_upper_folds"],4)

    def test_taxon_drop_reduces_species_denominator_not_just_uncertainty(self):
        x=mod.exclude_taxon_everywhere(self.plots,self.traits,self.names,
                                      "Cynodon dactylon (L.) Pers.")
        self.assertEqual(x["remaining_complete_source_sites"],4)
        self.assertTrue(x["dropped_taxon_was_life_span_known"])
        self.assertEqual(x["two_year_mean"]["Track__minus__Interior Plot 1"]["all"]["original_sites"],4)

    def test_taxon_drop_name_must_be_original(self):
        with self.assertRaisesRegex(ValueError,"not an original"):
            mod.exclude_taxon_everywhere(self.plots,self.traits,self.names,
                                         "invented plant taxon")

    def test_label_ablations_are_not_species_omission(self):
        r=mod.label_uncertainty_ablation(self.plots,self.traits,self.names)
        self.assertEqual(r["curated_baseline_upper"],mod.analyze(
            self.plots,self.traits,len(self.names))["comparisons"][
             "Track__minus__Interior Plot 1"]["all"]["two_year_mean"]["sharp_upper"])
        self.assertTrue(r["Cynodon dactylon (L.) Pers."]["drop_a_confirmed_label_not_drop_the_species"])

    def test_missing_plot_remains_gated(self):
        self.plots[("Station A",2017,"Interior Plot 2")]=None
        r=mod.produce(self.plots,self.traits,self.names)
        self.assertEqual(r["source_site_count"],3)
        self.assertEqual(r["all_zone_paired_references"]["Track__minus__Embankment"]["all"]["two_year_mean"]["original_sites"],3)

    def test_different_reference_zone_not_silent_averaged(self):
        r=mod.produce(self.plots,self.traits,self.names)
        contrasts=r["all_zone_paired_references"]
        self.assertIn("Track__minus__Interior Plot 2",contrasts)
        self.assertIn("Track__minus__Embankment",contrasts)
        self.assertNotIn("Track__minus__reference_averaged",contrasts)

if __name__=="__main__":
    unittest.main()
