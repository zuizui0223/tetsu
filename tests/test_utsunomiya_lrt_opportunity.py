"""Synthetic corridor exposure and source/photo audit regression tests."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_utsunomiya_lrt_photo_opportunity import (
    exposure_bin,organism_group,summarize_photo_year
)


def photo(species="Setaria viridis",accuracy=10):
    return {
        "id":1,
        "taxon":{"name":species},
        "quality_grade":"research",
        "photos":[{"id":1}],
        "observed_on":"2024-06-15",
        "obscured":False,"geoprivacy":None,"captive":False,
        "positional_accuracy":accuracy,
        "geojson":{"coordinates":[140.0,36.55]},
        "user":{"id":75}
    }


class UtsunomiyaAuditTests(unittest.TestCase):
    def test_fixed_rail_bins(self):
        self.assertEqual(exposure_bin(250),"near_0_250m")
        self.assertEqual(exposure_bin(251),"excluded_mid_or_outer")
        self.assertEqual(exposure_bin(1000),"control_1_3km")
        self.assertEqual(exposure_bin(3000),"control_1_3km")
        self.assertEqual(exposure_bin(3001),"excluded_mid_or_outer")

    def test_setaria_is_genus_not_single_taxon(self):
        self.assertEqual(organism_group(photo()),"Setaria")
        self.assertEqual(organism_group(photo("Solidago altissima")),"Solidago altissima")
        self.assertIsNone(organism_group(photo("Solidago gigantea")))
        self.assertEqual(organism_group(photo("Equisetum arvense")),"Equisetum arvense")

    def test_identical_photos_dont_create_new_independent_grid_cells(self):
        samples=[photo(),photo("Equisetum arvense")]
        result=summarize_photo_year(samples,lambda lon,lat:(120,(500,600)))
        self.assertEqual(result["photo_counts_by_rail_proximity"]["near_0_250m"],2)
        self.assertEqual(result["distinct_1km_cells_by_bin"]["near_0_250m"],1)
        self.assertEqual(result["distinct_photographers_by_bin"]["near_0_250m"],1)

    def test_bad_accuracy_cannot_be_rail_site(self):
        result=summarize_photo_year([photo(accuracy=200)],lambda lon,lat:(10,(1,2)))
        self.assertEqual(result["photo_counts_by_rail_proximity"].get("near_0_250m",0),0)
        self.assertEqual(result["api_geo_quality"]["accuracy_not_within_100m"],1)


if __name__=="__main__":
    unittest.main()
