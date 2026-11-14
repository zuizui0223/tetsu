"""v12.9 real author 2013 tracked CAR seed displacement experiment Table 2 source."""
from __future__ import annotations
import sys,tempfile,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_road_car_airflow_original_seed_morphology_v129 import (
    read_original_experiment_table,summarize_original_source,run
)

CSV=Path(__file__).resolve().parents[1]/"data"/"road_car_airflow_vonderlippe2013_original_table2_numeric_v129.csv"


class OriginalCarSlipstreamPhysicalSeedDistance(unittest.TestCase):
    def test_exact_published_4_species_5_pass_levels_and_group_medians(self):
        rows=read_original_experiment_table(CSV)
        self.assertEqual(len(rows),20)
        r=run(CSV)
        self.assertEqual(r["n_original_botanical_species"],4)
        self.assertEqual(r["n_original_experimental_replicates_per_species"],10)
        self.assertEqual(r["source_speed_kmh"],48)
        self.assertEqual(r["source_max_longitudinal_measurement_m"],45)
        p=r["original_table2_direct_ROAD_car_airflow_seed_movement_by_passes"]
        self.assertAlmostEqual(p["80"]["wind_type_vs_other_source_summary"][
            "wind_morphology"]["mean_of_published_species_median_displacement_m"],8)
        self.assertAlmostEqual(p["80"]["wind_type_vs_other_source_summary"][
            "no_wind_morphology"]["mean_of_published_species_median_displacement_m"],1)
        self.assertAlmostEqual(p["80"]["ratio_of_group_mean_source_medians"],8)
        self.assertAlmostEqual(p["80"]["ratio_of_group_means_original_mean_displacement"],
                               10.645/1.65)
        self.assertFalse(r["direct_railway_train_gust_experiment_or_cargo_seed_measurement"])

    def test_original_four_seed_species_and_cumulative_pass_observations_not_20_independent_biological_replicates(self):
        rows=read_original_experiment_table(CSV)
        species=set(x["species"] for x in rows)
        self.assertEqual(len(species),4)
        self.assertEqual(set(x["passes"] for x in rows),{1,10,20,40,80})
        self.assertTrue(all(sum(x["species"]==k for x in rows)==5 for k in species))

    def test_source_value_corruption_not_accepted_even_if_other_source_cells_unchanged(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"table2.csv"
            raw=CSV.read_text(encoding="utf-8")
            p.write_text(raw.replace("Ailanthus altissima,wind_morphology,1,5.14",
                                     "Ailanthus altissima,wind_morphology,1,5.15"),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"endpoint summaries"):
                read_original_experiment_table(p)

    def test_source_morphology_relabeling_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"table2.csv"
            raw=CSV.read_text(encoding="utf-8")
            p.write_text(raw.replace("Brassica napus,no_wind_morphology,",
                                     "Brassica napus,wind_morphology,"),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"morphology label"):
                read_original_experiment_table(p)

    def test_empty_original_source_rows_not_invented(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"table2.csv"
            rows=CSV.read_text(encoding="utf-8").splitlines()
            p.write_text("\n".join(rows[:-1])+"\n",encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"Missing genuine 20"):
                read_original_experiment_table(p)


if __name__=="__main__":unittest.main()
