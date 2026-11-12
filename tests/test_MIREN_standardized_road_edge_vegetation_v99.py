"""v9.9 tests: same plot-size/road/site/YEAR original MIREN
alien richness pairing, duplicate taxon and native-status conflicts,
first site/year robustness; no railway causal effects invented.
"""
import hashlib,json,sys,tempfile,unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_MIREN_standardized_road_edge_vegetation_v99 as model
from analyze_MIREN_standardized_road_edge_vegetation_v99 import (
    SITE_YEAR,SITE,STRICT_PLOTS,native_status_safe,
    plot_species_richness,matched_siteyears,region_summary,
    region_signs,repeated_site_temporal,compute,source_integrity)

def fixture():
    """Real source grain: 1 species per plot/year. Synthetic only."""
    rows=[]
    def add(region,road,site,plot,year,status,sp):
        rows.append({
            "Plot_id":f"{region}.{road}.{site:02d}.{plot}",
            "Region":region,"Road":road,"Site":site,
            "Plot":plot,"Year":year,"Elevation":550.,
            "Latitude":42.,"Longitude":-90.,
            "Accepted.Name.MIREN":sp,"Abundance":1.,
            "Cover":2.,"Status":status,"Protocol":"v.2017"})
    # two regions, two separate ROAD-sites each and two true revisits.
    for region in ("M1","M2"):
        for site in (1,2):
            for year in (2017,2022):
                for plot in ("1","2","3"):
                    for i in range(4 if plot=="1" else (2 if plot=="2" else 1)):
                        add(region,"R",site,plot,year,"Alien",f"Alien{i}")
                    for i in range(3):
                        add(region,"R",site,plot,year,"Native",f"Native{i}")
    # These exact duplicate lines are not independent taxa.
    rows.append(rows[0].copy())
    # This species has conflict of native status in region M2:
    # must NOT be decided by higher richness on any plot.
    add("M2","R",1,"1",2017,"Alien","ConflictX")
    add("M2","R",1,"3",2017,"Native","ConflictX")
    # replacement plot should not be treated as exact original site.
    add("M2","R",1,"1A",2017,"Alien","NewSiteTaxon")
    return pd.DataFrame(rows)

class MIRENIndependentRoadPlantSourceTests(unittest.TestCase):
    def test_same_year_and_same_road_site_are_necessary_no_fake_no_plant_plots(self):
        data=fixture()
        clean,q=native_status_safe(data)
        plot=plot_species_richness(clean)
        grid,paired,tri=matched_siteyears(plot)
        self.assertEqual(len(paired),8)
        self.assertEqual(len(tri),8)
        self.assertEqual(len(clean[clean.Plot=="1A"]),0)
        self.assertTrue(q["n_removed_repeated_species_same_plot_year"]>=1)
        self.assertEqual(q["n_region_species_with_ambiguous_native_alien_status_excluded_from_both"],1)
        self.assertEqual(paired["Alien_roadside_minus_interior"].min(),3)
    def test_observed_ROAD_pair_1_2_3_is_distance_gradient_not_train(self):
        clean,q=native_status_safe(fixture())
        _,pairs,triplets=matched_siteyears(plot_species_richness(clean))
        self.assertEqual(triplets["Alien_1"].median(),4)
        self.assertEqual(triplets["Alien_2"].median(),2)
        self.assertEqual(triplets["Alien_3"].median(),1)
        self.assertTrue((pairs["Alien_roadside_minus_interior"]>0).all())
    def test_first_site_visit_is_not_pseudoreplicated_with_2022_visit(self):
        clean,q=native_status_safe(fixture())
        _,paired,_=matched_siteyears(plot_species_richness(clean))
        first=paired.sort_values("Year").drop_duplicates(list(SITE))
        self.assertEqual(len(paired),8)
        self.assertEqual(len(first),4)
        reg=region_summary(first)
        self.assertEqual(len(reg),2)
        self.assertEqual(region_signs(reg,2)["n_codes_roadside_alien_richness_higher"],2)
    def test_temporal_first_last_when_no_change_is_not_false_adaptation(self):
        clean,q=native_status_safe(fixture())
        _,paired,_=matched_siteyears(plot_species_richness(clean))
        t=repeated_site_temporal(paired)
        self.assertEqual(t["n_original_road_sites_with_two_or_more_matched_survey_years"],4)
        self.assertEqual(t["median_original_road_site_change_latest_minus_first"],0)
        self.assertTrue(t["NO_road_management_or_railway_development_causal_trend"])
    def test_source_exact_MD5_and_row_provenance_required(self):
        df=fixture()
        with tempfile.TemporaryDirectory() as t:
            csv=Path(t)/"source.csv"
            audit=Path(t)/"audit.json"
            df.to_csv(csv,index=False)
            md5=hashlib.md5(csv.read_bytes()).hexdigest()
            oldmd5,oldrows=model.SOURCE_MD5,model.SOURCE_ROWS
            try:
                model.SOURCE_MD5,model.SOURCE_ROWS=md5,len(df)
                audit.write_text(json.dumps({
                    "complete":True,
                    "source_original_binary_metadata":{
                        "source_original_zenodo_md5":md5,
                        "source_original_Zenodo_DOI":model.DOI},
                    "actual_CSV_source_schema":{"n_original_csv_rows":len(df)}
                }))
                actual,actualmd5=source_integrity(csv,audit)
                self.assertEqual(actualmd5,md5)
                self.assertEqual(len(actual),len(df))
                audit.write_text('{"complete":false}')
                with self.assertRaises(KeyError):
                    source_integrity(csv,audit)
            finally:
                model.SOURCE_MD5,model.SOURCE_ROWS=oldmd5,oldrows
    def test_compute_reports_region_code_and_source_not_true_measured_rail(self):
        out=compute(fixture())
        self.assertTrue(out["complete"])
        self.assertEqual(out["geography"]["n_region_codes_with_original_road_plot1_and_interior_plot3_same_site_year"],2)
        self.assertEqual(out["geography"]["n_matched_original_plot1_plot3_same_site_year_pairs"],8)
        self.assertEqual(out["geography"]["n_original_STRICT_complete_three_plot_site_year_sets"],8)
        self.assertTrue(out["NO_train_seed_dispersal_or_mowing_wind_heritable_adaptation_conclusions"])

if __name__=="__main__":
    unittest.main()
