"""v17.4 all small numerical fixtures are SYNTHETIC SOURCE MATH TESTS
and NOT any actual taxonomic CzechVeg railway relevé, Preslia railway station
survey or evidence of original ecology/trait adaptation.
"""
from __future__ import annotations
import csv,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_original_CzechVeg_independent_railroad_plot_baseline_39rail_stations_v174 import (
    standardized_taxon_panel,rank_independent_rail_context,
    source_vegetation_scan,station_source,retain_blocks,
    SOURCE_CZECHVEG_PLANT_ROWS,MATCHING_S2
)

def toy_plot():
    block={
      ("paper1","region1",2020,"gt16_to_64"):{"rail":["r1","r2"],"road":["c1","c2"]},
      ("paper2","region1",2020,"gt16_to_64"):{"rail":["r3","r4","r5","r6"],"road":["c3","c4"]}
    }
    names={
       "r1":{"Poa annua","Echium vulgare"},
       "r2":{"Poa annua"},
       "c1":{"Echium vulgare"},
       "c2":set(),
       "r3":{"Echium vulgare"},
       "r4":set(),"r5":set(),"r6":set(),
       "c3":{"Poa annua"},"c4":set()
    }
    stations={"Poa annua":10,"Echium vulgare":5}
    background={"Poa annua":35,"Echium vulgare":8}
    return block,names,stations,background

class TrueOriginalCzechVegStudyPlotGrainSafety(unittest.TestCase):
    def test_actual_source_count_contracts_not_two_point_five_million_railway_areas(self):
        self.assertEqual(SOURCE_CZECHVEG_PLANT_ROWS,2510340)
        self.assertEqual(MATCHING_S2[0],"S2_study_region_decade_area4_min2")

    def test_source_block_same_study_region_decade_area_standardized_species_Q(self):
        blocks,names,stations,background=toy_plot()
        obj,support=standardized_taxon_panel(blocks,names,stations,background,100)
        # Poa: first block rail2/2-road0/2 =+1; second
        # rail0/4-road1/2=-.5 => w=2 both => +.25.
        self.assertAlmostEqual(
          obj["Poa annua"]["CzechVeg_study_region_decade_area_STANDARDIZED_rail_minus_road_plot_prevalence"],.25)
        # Echium: block1 rail1/2-road1/2=0; block2 rail1/4-road0/2=+.25
        self.assertAlmostEqual(
          obj["Echium vulgare"]["CzechVeg_study_region_decade_area_STANDARDIZED_rail_minus_road_plot_prevalence"],.125)
        self.assertEqual(support["n_authentic_same_study_plot_balance_blocks"],2)
        self.assertEqual(support["n_original_rail_text_plots"],6)
        self.assertEqual(support["n_original_road_text_plots"],4)
        self.assertEqual(support["n_original_source_field_studies"],2)
        self.assertEqual(obj["Poa annua"]["CzechVeg_same_study_source_N_distinct_publications_with_taxon"],2)
        self.assertTrue(obj["Poa annua"]["station_occurrence_ratio_NOT_a_genuine_rail_affinity_effect"])
        self.assertEqual(obj["Poa annua"]["CzechVeg_all_original_source_AUTHOR_STUDY_scored_loo_full_source_studies"],2)
        self.assertEqual(obj["Poa annua"]["CzechVeg_original_AUTHOR_STUDY_scored_loo_N_with_positive_rail_minus_road"],1)
        self.assertAlmostEqual(obj["Poa annua"]["CzechVeg_original_AUTHOR_STUDY_scored_loo_min_delta"],-.5)
        self.assertAlmostEqual(obj["Poa annua"]["CzechVeg_original_AUTHOR_STUDY_scored_loo_max_delta"],1.)
        self.assertFalse(obj["Poa annua"]["CzechVeg_original_AUTHOR_STUDY_scored_loo_POSITIVE_for_ALL_eligible_omissions"])
        self.assertTrue(obj["Poa annua"]["CzechVeg_source_loo_deletes_entire_original_publication_relevés_not_station_cohort"])

    def test_source_candidate_does_not_depend_on_photo_or_station_prevalence_ratio(self):
        blocks,names,stations,background=toy_plot()
        obj,_=standardized_taxon_panel(blocks,names,stations,background,100)
        candidates,support=rank_independent_rail_context(obj)
        # Poa source plots 3, no taxon has >=5 verified study-matched source
        # observed plot taxon records and >=2 independent studies simultaneously
        self.assertEqual(candidates,[])
        self.assertEqual(support["n_eligible_station_taxa_with_positive_independent_CzechVeg_railroad_proxy_and_source_support"],0)
        obj["Poa annua"]["CzechVeg_same_study_region_decade_area_source_N_RAIL_text_plot_with_this_taxon"]=5
        obj["Poa annua"]["CzechVeg_same_study_region_decade_area_source_N_ROAD_text_plot_with_this_taxon"]=1
        a,_=rank_independent_rail_context(obj)
        self.assertEqual([x["scientific_name_original_exact_normalized"] for x in a],["Poa annua"])
        self.assertTrue(a[0]["CzechVeg_observational_proxy_NOT_independently_confirmed_active_railway_affinity"])

    def test_one_rail_species_plot_discovery_is_not_robust_twenty_taxon_ranking(self):
        obj={
          "A":{
            "CzechVeg_study_region_decade_area_STANDARDIZED_rail_minus_road_plot_prevalence":.20,
            "station_original_2020_22_n_of_39":39,
            "CzechVeg_same_study_region_decade_area_source_N_RAIL_text_plot_with_this_taxon":1,
            "CzechVeg_same_study_region_decade_area_source_N_ROAD_text_plot_with_this_taxon":0,
            "CzechVeg_same_study_source_N_distinct_publications_with_taxon":1,
          }
        }
        rows,meta=rank_independent_rail_context(obj)
        self.assertEqual(len(rows),0)
        self.assertEqual(meta["source_selected_top20_max_not_forced_to_fill"],0)

    def test_full_2510340_author_species_source_records_required_not_fake_small_CSV(self):
        with tempfile.TemporaryDirectory() as temp:
            f=Path(temp)/"toy.csv"
            with f.open("w",encoding="utf-8",newline="") as out:
                writer=csv.DictWriter(out,fieldnames=("RELEVE_NR","KAPLAN_PLUS_NAME_SIMPLE","PLANT_GROUP"))
                writer.writeheader()
                writer.writerow({"RELEVE_NR":"r1","KAPLAN_PLUS_NAME_SIMPLE":"Poa annua","PLANT_GROUP":"vascular"})
            with self.assertRaisesRegex(ValueError,"2.51034 million"):
                source_vegetation_scan(f,{"r1":{"area":2}},{"r1"},{"Poa annua"})

    def test_wrong_Preslia_39_railway_stations_cannot_be_used_as_actual_source(self):
        mocked={"AA":{"Poa annua"},"BB":{"Echium vulgare"}}
        with patch(
          "analyze_original_CzechVeg_independent_railroad_plot_baseline_39rail_stations_v174.read_original_czech_publisher",
          return_value=(mocked,{},{"original_source_Czech_railway_Preslia_DOI":"10.23855/preslia.2024.247"})
        ):
            with self.assertRaisesRegex(ValueError,"39 railway stations"):
                station_source("fakeS1","fakeS2")

if __name__=="__main__":unittest.main()
