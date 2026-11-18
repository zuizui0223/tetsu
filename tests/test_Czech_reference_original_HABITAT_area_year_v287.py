"""v28.7 synthetic source audit: HABITAT text classification independent of
species and syntaxon, preselected true source survey year / physical plot
area / altitude availability, matching eligibility without data leakage."""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Czech_original_14station_reference_habitat_era_area_altitude_v287 as m


class CzechOriginalHabitatReferenceMetadataQA(unittest.TestCase):
    def test_open_and_gravel_habitat_categories_are_not_all_bare_grass(self):
        self.assertTrue(bool(m.OPEN_WORDS.search("suchá louka a travnatý porost")))
        self.assertTrue(bool(m.OPEN_WORDS.search("ruderal railway-free open gravel grassland")))
        self.assertTrue(bool(m.GRAVEL_WORDS.search("písčitý štěrkový násep")))
        self.assertFalse(bool(m.GRAVEL_WORDS.search("meadow pasture")))
        self.assertTrue(bool(m.NEGATIVE_WORDS.search("forest woodland peatland")))
        # No botanical SYNTAXON or species output is consulted anywhere
        self.assertFalse(bool(m.OPEN_WORDS.search("TH original publisher vegetation syntaxon")))

    def test_metadata_conditions_all_fail_closed_for_missing_year_area_altitude(self):
        record={"open":True,"gravel":True,"year":None,
                "plot_area_m2":None,"altitude_m":None}
        self.assertTrue(m.belongs("habitat_text_open",record))
        self.assertTrue(m.belongs("explicit_rock_sand_gravel",record))
        self.assertFalse(m.belongs("habitat_open_recent",record))
        self.assertFalse(m.belongs("habitat_open_area",record))
        self.assertFalse(m.belongs("open_recent_area",record))
        self.assertFalse(m.belongs("broad_open_altitude",record))
        record["year"]=2022.
        record["plot_area_m2"]=100.
        record["altitude_m"]=435.
        self.assertTrue(m.belongs("habitat_open_recent",record))
        self.assertTrue(m.belongs("open_recent_area",record))
        self.assertTrue(m.belongs("broad_open_altitude",record))
        with self.assertRaises(ValueError):m.belongs("known_best_biology_positive",record)

    def test_same_real_16_rail_station_cohort_and_publication_requirement(self):
        stations=[]
        metadata={}
        for rail in range(16):
            candidates=[]
            for j in range(6):
                key=f"plot{rail}_{j}"
                candidates.append({"id":key,"study":f"independent_paper{j%3}"})
                metadata[key]={"open":j<5,"gravel":j<2,
                    "year":2019 if j<2 else 2021,
                    "plot_area_m2":25.,"altitude_m":450.}
            stations.append({"id":f"rail{rail}","stretch":str(rail%4),
                             "control_options":candidates})
        same,qs=m.panel_stations(stations,metadata,"baseline_any_nontransport_context")
        open_,q=m.panel_stations(stations,metadata,"habitat_text_open")
        gravel,qg=m.panel_stations(stations,metadata,"explicit_rock_sand_gravel")
        self.assertEqual(qs["eligible_source_reference_plot_opportunities"],96)
        self.assertEqual(q["eligible_source_reference_plot_opportunities"],80)
        self.assertEqual(q["original_rail_stations_with_both_min5_and_min2papers"],16)
        self.assertEqual(qg["original_rail_stations_with_both_min5_and_min2papers"],0)
        self.assertEqual(len(same),16)
        self.assertEqual(len(open_),16)
        self.assertEqual(len(gravel),16)
        # Main railway identities unchanged; alternative source controls
        # are chosen without inspecting original botanical taxa.
        self.assertEqual([x["id"] for x in stations],[x["id"] for x in gravel])

    def test_metadata_reporting_does_not_serialize_exact_original_plot_ids(self):
        candidates=[{"id":"SECRETOriginalPlot17","study":"paper1"},
                    {"id":"SECRETOriginalPlot18","study":"paper2"}]
        mdata={x["id"]:{"open":True,"gravel":False,"year":2022.,
                       "plot_area_m2":25.,"altitude_m":420.,
                       "slope_deg":2.,
                       "source_original_habitat_present":True,
                       "source_geology_text_recorded":False,
                       "source_soil_text_recorded":False}
               for x in candidates}
        result=m.chosen_control_metadata({0:candidates},
                      [{"id":"station0","stretch":"I","control_options":candidates}],mdata)
        self.assertEqual(result["source_reference_original_releve_count"],2)
        self.assertEqual(result["original_railway_field_stretches"],1)
        self.assertEqual(result["source_control_original_plot_area_m2_median"],25)
        self.assertNotIn("SECRETOriginalPlot17",str(result))
        self.assertNotIn("paper1",str(result))

    def test_source_1D_and_2D_metrics_not_interchangeable(self):
        self.assertEqual(len(m.MEASURES),3)
        self.assertIn("CWM_height_1D_absolute",m.MEASURES)
        self.assertIn("CWM_seed_mass_1D_absolute",m.MEASURES)
        self.assertNotEqual(m.MEASURES[1],m.MEASURES[2])


if __name__=="__main__":
    unittest.main()
