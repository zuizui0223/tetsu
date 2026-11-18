"""v26.9 synthetic mathematical QA for source-authentic same-site RAIL/ROAD beta."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_original_RAIL_vs_global_ROAD_taxonomic_microhabitat_beta_v269 as v
from analyze_rail_microhabitat_replacement_nestedness_v256 import components as prior


class RailRoadBetaTests(unittest.TestCase):
    def test_pure_nestedness_not_species_replacement(self):
        x=v.beta({1,2},{1,2,3,4})
        self.assertEqual(x["beta_jtu_replacement"],0.)
        self.assertEqual(x["beta_jne_nestedness"],.5)
        self.assertEqual(x["corridor_only_taxa"],0.)
        self.assertEqual(x["corridor_strict_subset"],1.)

    def test_only_replacement_and_exact_source_v256_equivalence(self):
        x=v.beta({1,2},{3,4})
        self.assertEqual(x["beta_jac"],1.)
        self.assertEqual(x["beta_jtu_replacement"],1.)
        self.assertEqual(x["beta_jne_nestedness"],0.)
        q=prior({1,2},{3,4})
        self.assertAlmostEqual(q["beta_jac"],x["beta_jac"])
        self.assertAlmostEqual(q["beta_replacement_jtu"],x["beta_jtu_replacement"])

    def test_mixed_taxonomic_beta_partition(self):
        x=v.beta({1,2,3},{2,3,4,5})
        self.assertAlmostEqual(x["beta_jac"],.6)
        self.assertAlmostEqual(x["beta_jtu_replacement"],.5)
        self.assertAlmostEqual(x["beta_jne_nestedness"],.1)
        self.assertAlmostEqual(x["replacement_fraction"],5/6)
        self.assertAlmostEqual(x["corridor_unshared_fraction"],1/3)

    def test_identical_and_invalid_empty(self):
        self.assertIsNone(v.beta({1},{1})["replacement_fraction"])
        with self.assertRaises(ValueError):v.beta(set(),{1})

    def test_source_site_rarefaction_exactly_equalizes_richness(self):
        a={1,2,3,4,5};b={2,4,5,6,7,8,9}
        q=v.rarefied_site(a,b,721,79)
        self.assertEqual(q["corridor_richness"],5.)
        self.assertEqual(q["reference_richness"],5.)
        self.assertAlmostEqual(q["beta_jac"],q["beta_jtu_replacement"]+q["beta_jne_nestedness"])
        self.assertEqual(a,{1,2,3,4,5})
        self.assertEqual(b,{2,4,5,6,7,8,9})
        self.assertEqual(q,v.rarefied_site(a,b,721,79))

    def test_original_first_year_filter_frozen_not_future_imputation(self):
        x=[{"year":2007,"region":"R"},{"year":2012,"region":"R"},{"year":2017,"region":"R"}]
        self.assertEqual(len(v.select_miren_sites(x,"first_anyyear")),3)
        self.assertEqual(len(v.select_miren_sites(x,"first_2012plus")),2)
        self.assertEqual(len(v.select_miren_sites(x,"first_2017plus")),1)
        with self.assertRaises(ValueError):v.select_miren_sites(x,"future_only_not_source")

    def test_road_region_two_stage_equal_weight_and_within_site(self):
        rows=[]
        for region,dif in (("R1",0),("R2",1)):
            for road in ("A","B"):
                for site in range(4 if road=="A" else 7):
                    rows.append({"region":region,"road":road,"site":str(site),"year":2017,
                        "corridor":({1,2} if dif==0 else {1,4}),
                        "reference":{1,2,3}})
        out=v.road_macro(rows,"first_anyyear",7)
        self.assertEqual(out["original_regions_meeting_two_road_support"],2)
        self.assertEqual(out["original_roads_in_supported_regions"],4)
        self.assertEqual(out["original_sites_in_supported_regions"],22)
        self.assertAlmostEqual(
            out["observed_mean_equal_road_equal_region"]["raw"]["beta_jac"],
            (1/3+3/4)/2)
        self.assertEqual(out["observed_mean_equal_road_equal_region"]["raw"]["n_region_codes"],2)


if __name__=="__main__":unittest.main()
