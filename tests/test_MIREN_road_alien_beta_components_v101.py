"""v10.1 genuine field-survey source-method tests for roadside vs interior
Alien + Native plant Jaccard beta = turnover + nestedness, with strict
source-plot matching and no fabricated railway disturbance inference.
Synthetic fixtures only; full original 167804 rows run in GitHub Action.
"""
import math,sys,unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_MIREN_road_alien_beta_components_v101 import (
    jaccard_components,original_spatial_incidence,group_summary,
    region_summary,SOURCE_GROUPS,REGION_BOTH_POSITIVE_MIN)


def fixture():
    rows=[]
    for region in ("R1","R2"):
        for site in range(1,13):
            for year in (2017,2022):
                for plot in ("1","3"):
                    for name in (
                        ("a","b","c") if plot=="1" else ("b","d")):
                        rows.append({
                            "Region":region,"Road":"ROAD",
                            "Site":site,"Year":year,"Plot":plot,
                            "Accepted.Name.MIREN":name,
                            "Status":"Alien"})
                    for name in (
                        ("n1","n2","n3","n4") if plot=="1" else
                        ("n1","n2","n3","n5")):
                        rows.append({
                            "Region":region,"Road":"ROAD",
                            "Site":site,"Year":year,"Plot":plot,
                            "Accepted.Name.MIREN":name,
                            "Status":"Native"})
    rows.append(rows[0].copy())
    rows.append({**rows[0],"Plot":"1A","Accepted.Name.MIREN":"EXCLUDE"})
    return pd.DataFrame(rows)


class MIRENJaccardSourceComponents(unittest.TestCase):
    def test_predefined_alien_native_not_invented_train_mode(self):
        self.assertEqual(SOURCE_GROUPS,("Alien","Native"))
        self.assertEqual(REGION_BOTH_POSITIVE_MIN,(10,20))
    def test_asymmetric_shared_and_unique_sets_conserve_jaccard(self):
        z=jaccard_components({"a","b","c"},{"b","d"})
        self.assertEqual((z["shared_species"],z["roadside_only_species"],
                          z["interior_only_species"]),(1,2,1))
        self.assertAlmostEqual(z["jaccard_total"],.75)
        self.assertAlmostEqual(z["jaccard_turnover"],2/3)
        self.assertAlmostEqual(z["jaccard_nestedness_resultant"],1/12)
        self.assertTrue(z["turnover_higher_than_nestedness"])
    def test_nestedness_without_turnover_and_disjoint_replacement(self):
        nested=jaccard_components({"a","b","c"},{"a","b"})
        self.assertEqual(nested["jaccard_turnover"],0)
        self.assertAlmostEqual(nested["jaccard_nestedness_resultant"],1/3)
        disjoint=jaccard_components({"a","b"},{"c","d"})
        self.assertEqual(disjoint["jaccard_total"],1)
        self.assertEqual(disjoint["jaccard_turnover"],1)
        self.assertEqual(disjoint["jaccard_nestedness_resultant"],0)
    def test_unknown_alien_group_in_one_plot_not_fabricated_zero_turnover(self):
        z=jaccard_components({"a","b"},set())
        self.assertFalse(z["has_both_recorded_species_pools"])
        self.assertIsNone(z["jaccard_total"])
        self.assertIsNone(z["jaccard_turnover"])
        self.assertIsNone(z["jaccard_nestedness_resultant"])
        self.assertEqual(z["roadside_minus_interior_richness"],2)
        self.assertFalse(jaccard_components(set(),set())[
            "has_both_recorded_species_pools"])
    def test_repeated_siteyears_not_pseudoreplicated_in_primary(self):
        clean,qc,allpairs,first=original_spatial_incidence(fixture())
        self.assertEqual(len(allpairs),48)
        self.assertEqual(len(first),24)
        self.assertEqual(len(clean[clean.Plot=="1A"]),0)
        self.assertEqual(qc["n_removed_repeated_species_same_plot_year"],1)
        self.assertTrue((allpairs.Alien_roadside_only_species==2).all())
        self.assertTrue((allpairs.Alien_interior_only_species==1).all())
        self.assertTrue(allpairs.Alien_jaccard_total.notna().all())
    def test_region_requires_min_unique_first_source_sites_and_native_2017_possible(self):
        clean,qc,allpairs,first=original_spatial_incidence(fixture())
        report=group_summary(first,"Alien")
        self.assertEqual(report["n_same_original_site_year_plot1_vs_plot3_comparisons"],24)
        self.assertEqual(report["n_both_original_plot_group_taxa_nonzero"],24)
        self.assertEqual(report["n_source_siteyears_without_defined_Jaccard_component"],0)
        self.assertAlmostEqual(report["group_both_positive_mean_Jaccard_turnover"],2/3,6)
        self.assertEqual(region_summary(first,"Alien",10)["n_region_codes_turnover_dominant"],2)
        self.assertEqual(region_summary(first,"Alien",20)["n_region_codes_with_min_support"],0)
        self.assertEqual(region_summary(first[first.Year>=2017],"Native",10)[
            "n_region_codes_with_min_support"],2)
    def test_missing_complete_pair_site_not_zero_species(self):
        data=fixture()
        data=data[~((data.Region=="R1")&(data.Site==1)&
                    (data.Year==2017)&(data.Plot=="3"))]
        _,qc,allpairs,first=original_spatial_incidence(data)
        self.assertEqual(len(allpairs),47)
        self.assertEqual(len(first),24)
        self.assertFalse("unknown_species_absences" in qc)
    def test_jaccard_site_plot_groups_are_species_sets_not_cover_ordinal(self):
        d=fixture()
        d["Cover"]="3"
        _,_,pairs,first=original_spatial_incidence(d)
        a=group_summary(first,"Alien")
        self.assertAlmostEqual(a["mean_roadside_only_recorded_group_taxa_per_pair"],2)
        self.assertAlmostEqual(a["mean_interior_only_recorded_group_taxa_per_pair"],1)

if __name__=="__main__":unittest.main()
