"""v17.1 synthetic test: changing photo-rights cohort BEFORE matching original
Echium CC0/BY public images, with no observed botanical flower/height labels.
"""
import sys,unittest,copy
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tests"))
from test_Echium_original_licensed_pixel_review_v170 import (
   source_fixtures,fake_source_photos,SourceRailStub
)
from reselect_Echium_CC0BY_before_photo_pairing_v171 import (
   licensed_first_matching,run,V171
)

class RightsFirstRematchingOriginalSourceTests(unittest.TestCase):
    def test_exact_source_original29_CC_all_first_then_rights_first_29(self):
        a,b=source_fixtures()
        orig=fake_source_photos()
        records,qc=licensed_first_matching(
            orig,SourceRailStub(),
            b["source_species_country_counts"][0],
            a["source_restricted_original_4_taxa_matched_opportunity"][0])
        self.assertEqual(qc["n_source_original_v168_allCC_pairs_reproduced"],29)
        self.assertEqual(qc["n_source_original_near30_CC0_BY_observations"],29)
        self.assertEqual(qc["n_source_original_rail100to1000_CC0_BY_observations"],29)
        self.assertEqual(len(records),29)
        self.assertEqual(qc["n_source_rights_first_CC0_BY_same_yearmonth_10km_different_observer_pairs"],29)

    def test_original_pairing_29_still_invariant_if_some_pictures_are_NC_and_not_redistributable_in_artifact(self):
        a,b=source_fixtures()
        orig=fake_source_photos()
        orig[0]["photos"][0]["license_code"]="cc-by-nc"
        pairs,qc=licensed_first_matching(
            orig,SourceRailStub(),
            b["source_species_country_counts"][0],
            a["source_restricted_original_4_taxa_matched_opportunity"][0])
        self.assertEqual(qc["n_source_original_v168_allCC_pairs_reproduced"],29)
        self.assertEqual(qc["n_source_original_near30_CC0_BY_observations"],28)
        self.assertEqual(len(pairs),28)
        self.assertGreater(qc["source_in_band_but_no_CC0BY_image_or_attribution"],0)

    def test_bad_original_source_observation_total_not_counted_as_no_pair(self):
        a,b=source_fixtures()
        with self.assertRaisesRegex(ValueError,"4,235"):
            licensed_first_matching(
               fake_source_photos()[:-1],SourceRailStub(),
               b["source_species_country_counts"][0],
               a["source_restricted_original_4_taxa_matched_opportunity"][0])

if __name__=="__main__":unittest.main()
