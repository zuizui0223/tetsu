#!/usr/bin/env python3
"""v359b: species frequency-preserving source null verifies exact margins, no Jaccard result synthesized."""
import random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Slovenia_same_line_geodistance_matched_frequency_null_v281 import (
  source_groups,exact_source_margins,swap_within_source_groups
)
from analyze_Slovenia_networkmatched_species_frequency_null_v359b import (
    EXPECTED_GEO_NETWORK_ANCHORS,EXPECTED_GEO_NETWORK_AUTHOR_SECTIONS,
    EXPECTED_MATCHED_NETWORK_TAXONOMIC_GAP,PANEL_NULLS,DRAWS,BURN,THIN
)

class V359bTrueSourceFrequencyNullTests(unittest.TestCase):
    def setUp(self):
        self.sites=[{"year":2022 if i<6 else 2023,"layout":"SINGLE" if i%2 else "DOUBLE",
             "section":str(i%4)} for i in range(12)]
        self.original=[{f"species{k}" for k in range(10) if ((i+3*k)%7)<3}
                       for i in range(12)]

    def test_same_original_n27_line_n7_frozen_and_two_null_families(self):
        self.assertEqual(EXPECTED_GEO_NETWORK_ANCHORS,27)
        self.assertEqual(EXPECTED_GEO_NETWORK_AUTHOR_SECTIONS,7)
        self.assertAlmostEqual(EXPECTED_MATCHED_NETWORK_TAXONOMIC_GAP,0.027033311999075067)
        self.assertEqual(len(PANEL_NULLS),2)
        self.assertEqual((DRAWS,BURN,THIN),(99,3000,350))

    def test_original_source_species_frequency_by_year_layout_held_exact(self):
        for family in PANEL_NULLS:
            groups=source_groups(self.sites,family)
            before=exact_source_margins(self.original,groups)
            after=[set(x) for x in self.original]
            swap_within_source_groups(after,groups,random.Random(20261124),2000)
            self.assertEqual(exact_source_margins(after,groups),before)

    def test_source_within_year_layout_does_not_move_taxa_across_original_study_year(self):
        groups=source_groups(self.sites,"within_original_survey_year_and_track_layout")
        for group in groups:
            self.assertEqual(len({(self.sites[i]["year"],self.sites[i]["layout"]) for i in group}),1)
        self.assertEqual(sum(len(group) for group in groups),12)

if __name__=="__main__":unittest.main()
