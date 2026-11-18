"""v27.1 synthetic MIREN source first-year-block fixed margin and chain test."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_MIREN_year_stratified_frequency_null_v271 as v


def fake_original_sites():
    sites=[]
    for region in ("R1","R2"):
        for road in ("a","b"):
            for i in range(8):
                sites.append({"region":region,"road":road,"site":str(i),
                              "year":2012 if i<4 else 2017,
                              "corridor":{1,2,i+10},
                              "reference":{1,i+10,i+60}})
    return sites


class YearSpecificRoadNullTests(unittest.TestCase):
    def test_source_original_road_by_first_year_and_support(self):
        blocks,q=v.blocks_by_source_year(fake_original_sites())
        self.assertEqual(len(blocks),2)
        self.assertEqual(q["n_original_source_roads"],4)
        self.assertEqual(q["n_source_regions_with_min_two_original_roads"],2)
        self.assertEqual(q["n_road_year_strata_with_four_original_sites"],8)
        self.assertEqual(q["n_original_spatial_sites"],32)
        for reg,roads in blocks.items():
            for road,strata in roads.items():
                self.assertEqual(len(strata),2)
                for sites in strata:
                    self.assertEqual(len(sites),4)
                    self.assertEqual(len({row["year"] for row in sites}),1)

    def test_source_2007_alien_focused_sites_excluded(self):
        source=fake_original_sites()
        for reg in ("R1","R2"):
            source.append({"region":reg,"road":"a","site":"new2007",
                           "year":2007,"corridor":{1},"reference":{1,2}})
        _,qa=v.blocks_by_source_year(source)
        self.assertEqual(qa["n_source_original_first_2012plus_site_rows"],32)
        self.assertEqual(qa["n_original_spatial_sites"],32)

    def test_original_site_nesting_road_region_effect_means(self):
        source=fake_original_sites()
        blocks,_=v.blocks_by_source_year(source)
        out=v.summarize_groups(blocks,1001,draws=19,burnin=100,spacing=20)
        self.assertTrue(out["identified"])
        self.assertEqual(out["n_region_codes"],2)
        self.assertEqual(out["n_roads"],4)
        self.assertEqual(out["n_sites"],32)
        self.assertEqual(out["n_source_road_year_strata"],8)
        self.assertTrue(out["source_fixed_site_and_species_frequency_results"]["beta_jac"]["observed"]<=1.)
        self.assertEqual(out["source_fixed_site_and_species_frequency_results"]["beta_jac"]["frequency_null_mean"],
                         v.summarize_groups(blocks,1001,draws=19,burnin=100,spacing=20)["source_fixed_site_and_species_frequency_results"]["beta_jac"]["frequency_null_mean"])

    def test_first_2012plus_road_cohort_cannot_include_future_other_group(self):
        data=fake_original_sites()
        # Legacy source one site's earlier original year must never
        # enter the modern first-year selected panel.
        data[0]["year"]=2007
        blocks,q=v.blocks_by_source_year(data)
        self.assertEqual(q["n_source_original_first_2012plus_site_rows"],31)
        self.assertEqual(q["n_original_spatial_sites"],28)
        # One 2012 group in one road now has only three sites,
        # its original 2017 block remains but the site is excluded.
        self.assertEqual(q["n_road_year_strata_with_four_original_sites"],7)


if __name__=="__main__":unittest.main()
