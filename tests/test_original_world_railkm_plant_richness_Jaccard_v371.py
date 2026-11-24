"""v371 original-global source schema/blocked modelling synthetic validation."""
import unittest,sys,math,random
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_original_world_railkm_plant_richness_Jaccard_v371 as m

class WorldPlantRailModel(unittest.TestCase):
    def test_original_rail_and_plant_same_global_grid_requires_source_integrity(self):
        fakep={"source_summary":{"original_source_iNat_obs_S3_ETag":m.OBS_ETAG,
                      "original_source_iNat_taxa_S3_ETag":m.TAXA_ETAG,
                      "source_global_photo_observations_QA":{
                        "original_global_iNat_photo_linked_rows_read":60000000}},
               "source_grids":{w:{
                    "source_full_world_100km_or_200km_gridded_plant_effort_and_taxonomic_counts":[
                        {"source_year_2025_photo_linked_plants_observation_count":100}]
                    } for w in ("100000","200000")}}
        faker={"source_original_fullworld_rail_128_files_SHA_verified":True,
               "source_Oveture_release":m.RAIL_RELEASE,
               "source_fullworld_candidate_normal_rail_km":m.RAIL_SOURCE_KM,
               "source_global_gridded_rail_exposure_km_by_cell_width":{
                  w:[{"source_real_original_Overture_normal_rail_km_clipped_into_cell":m.RAIL_SOURCE_KM}]
                  for w in ("100000","200000")}}
        out=m.verify_original_sources(fakep,faker)
        self.assertEqual(out["original_2025_full_world_plant_detections_in_grid"],100)
        fakep["source_summary"]["source_global_photo_observations_QA"][
           "original_global_iNat_photo_linked_rows_read"]=6000
        with self.assertRaises(ValueError):
            m.verify_original_sources(fakep,faker)
    def test_no_observer_effort_treated_as_no_plants(self):
        self.assertFalse(m.valid_source_count(50,4,30))
        self.assertFalse(m.valid_source_count(100,5,19))
        self.assertTrue(m.valid_source_count(50,5,20))
    def test_grid_transform_source_origin_shared_across_100km_200km(self):
        self.assertEqual(m.cell_center(0,0,100000),(-17350000,-7350000))
        self.assertEqual(m.cell_center(0,0,200000),(-17300000,-7300000))
        self.assertEqual(m.block_id(1,2,100000),m.block_id(0,1,200000))
    def test_original_rail_observation_cell_does_not_impute_plant_absences(self):
        plants=[{"source_equal_area_grid_ix":10,"source_equal_area_grid_iy":20,
          "source_year_2025_photo_linked_plants_observation_count":100,
          "source_2025_distinct_observers_with_photo_linked_plant_records":10,
          "source_2025_plant_species_richness":40}]
        climate={(10,20):{"continent":"Europe","block":(0,0),
                          "bio1":10,"abs_lat":48}}
        out,qa,rail=m.original_cell_rows(plants,[],climate)
        self.assertEqual(len(out),1)
        self.assertEqual(out[0]["rail_km"],0)
        self.assertEqual(out[0]["species"],40)
    def test_original_beta_record_requires_same_flora_source_cell_counts(self):
        plantrows=[]
        for ix in (10,11):
            plantrows.append({"source_equal_area_grid_ix":ix,
                  "source_equal_area_grid_iy":20,
                  "source_year_2025_photo_linked_plants_observation_count":100,
                  "source_2025_distinct_observers_with_photo_linked_plant_records":10,
                  "source_2025_plant_species_richness":40})
        geo={(ix,20):{"continent":"Europe","block":(0,0),
                        "bio1":10,"abs_lat":48} for ix in (10,11)}
        x,qa,rail=m.original_cell_rows(plantrows,[],geo)
        source=[{"source_cell_A_original_grid_ix_iy":[10,20],
                 "source_cell_B_original_grid_ix_iy":[11,20],
                 "source_biological_Jaccard_similarity_between_ENTIRE_landscape_cells":.25,
                 "source_grid_A_original_plant_species_n":40,
                 "source_grid_B_original_plant_species_n":40}]
        rows,qq=m.pair_rows(source,{a["key"]:a for a in x},rail)
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]["jaccard"],.25)
        source[0]["source_grid_B_original_plant_species_n"]=30
        with self.assertRaises(ValueError):m.pair_rows(source,{a["key"]:a for a in x},rail)
    def test_original_block_equal_weighted_not_oversample_dense_block(self):
        rows=[{"continent":"Europe","block":(0,0)}]*10 + [
             {"continent":"Europe","block":(1,0)}]
        w=m.block_weights(rows)
        # True frozen ridge objective is normalized block weighted MSE.
        self.assertAlmostEqual(float(sum(w)),1)
        self.assertAlmostEqual(float(sum(w[:10])),.5)
        self.assertAlmostEqual(float(w[-1]),.5)
    def test_ridge_holdout_correctly_improves_if_global_true_rail_signal(self):
        rows=[]
        for ci,cont in enumerate(m.CONTINENTS):
            for i in range(140):
                rail=float((i*7+ci)%29)
                obs=100+i%30
                visitors=10+i%8
                temp=10+ci%4+i%5
                lat=abs(15+ci*8+i%3)
                richness=20+.1*obs+.05*visitors+.4*rail
                features={"continent":cont,
                   "block":(i//10,ci),"log_species":math.log1p(richness),
                   "log_observations":math.log1p(obs),
                   "log_observers":math.log1p(visitors),
                   "bio1":temp,"bio1_sq":temp**2,
                   "abs_lat":lat,"abs_lat_sq":lat**2,
                   "log_rail_km":math.log1p(rail)}
                # Ensure out-of-continent holds actual rail effect rather than
                # a coincidental within-one-climate source pattern.
                features["log_species"]=(
                    2+.12*features["log_observations"]+
                    .06*features["log_observers"]+
                    .03*features["bio1"]+
                    .20*features["log_rail_km"])
                rows.append(features)
        stats=m.continent_holdout(rows,m.RICH_BASE,m.RICH_RAIL,"log_species",
            False,m.MIN_CONTINENT_RICH_CELLS,300,"richness")
        self.assertTrue(stats["independent_cross_continent_test_supported"])
        self.assertEqual(stats["n_original_continents_eligible_as_heldout_folds"],6)
        self.assertEqual(stats["n_holdout_continents_with_positive_delta_MSE"],6)
        self.assertGreater(stats["mean_holdout_continent_delta_MSE"],0)
    def test_inadequate_original_continents_fail_closed(self):
        rows=[{"continent":"Europe","block":(1,1),"log_species":3,
            **{f:1. for f in m.RICH_RAIL}} for _ in range(200)]
        stats=m.continent_holdout(rows,m.RICH_BASE,m.RICH_RAIL,
            "log_species",False,100,300,"richness")
        self.assertFalse(stats["independent_cross_continent_test_supported"])
        self.assertIsNone(stats["mean_holdout_continent_delta_MSE"])
    def test_EIVE_plant_temperature_not_conflated_with_site_WorldClim(self):
        self.assertEqual(m.SOURCE_YEAR,2025)
        self.assertEqual(m.RICH_RAIL[-1],"log_rail_km")
        self.assertNotIn("species_thermal_affinity",m.RICH_RAIL)

if __name__=="__main__":unittest.main()
