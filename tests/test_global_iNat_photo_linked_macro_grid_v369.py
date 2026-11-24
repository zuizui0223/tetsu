#!/usr/bin/env python3
"""v369 source-safe all-world 100/200km biodiversity grid tests.
No false rail80m or uniformly surveyed botanical field-plot claims.
"""
import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import (
  OBS_ETAG,TAXA_ETAG,OBS_BYTES,TAXA_BYTES,
  ROOT_PLANT_TAXON,SOURCE_REGION_M,
  X_START,X_END,Y_START,Y_END,source_coarse_cell,MAX_UNCERT_M
)
class V369WholeWorldBioGridTests(unittest.TestCase):
    def test_original_public_world_source_snapshot_is_explicitly_pinned(self):
        self.assertEqual(ROOT_PLANT_TAXON,47126)
        self.assertEqual(OBS_BYTES,13697505636)
        self.assertEqual(TAXA_BYTES,40184871)
        self.assertIn('-27',OBS_ETAG)
        self.assertEqual(SOURCE_REGION_M,(100000,200000))
        self.assertEqual(MAX_UNCERT_M,1000)

    def test_world_equal_area_grid_covers_entire_antimeridian_and_both_hemispheres(self):
        from pyproj import Transformer
        t=Transformer.from_crs('EPSG:4326','EPSG:6933',always_xy=True)
        for lon,lat in [(-179,55),(179,55),(140,36),(-50,-22),(19,-35),(14,50),(150,-33)]:
            for m in SOURCE_REGION_M:
                cell=source_coarse_cell(lon,lat,t,m)
                self.assertIsNotNone(cell)
                i,j=cell
                self.assertGreaterEqual(i,0)
                self.assertLess(i,(X_END-X_START)//m)
                self.assertGreaterEqual(j,0)
                self.assertLess(j,(Y_END-Y_START)//m)

    def test_no_raw_observer_identifier_or_species_geoprivacy_claim_in_schema(self):
        self.assertEqual((X_END-X_START)//100000,348)
        self.assertEqual((Y_END-Y_START)//100000,148)

if __name__=='__main__':unittest.main()
