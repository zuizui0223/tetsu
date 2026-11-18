"""v28.0 synthetic source nearest reference controls; never empirical plant findings."""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Slovenia_source_same_line_vs_other_geodistance_matched_v280 as m


def make_source():
    # 5 synthetic rail lines x 3 locations: many genuinely near
    # alternative-source line options for geodesic matching.
    sites=[]
    for line in range(5):
        for site in range(3):
            sites.append({
                "source_id":f"line{line}_site{site}","section":str(line),
                "class":"main_G" if line%2 else "regional_R",
                "lat":46.+.008*line,"lon":14.+.033*site,
                "routekm":site*3.,"year":2022 if line%2 else 2023,
                "layout":"SINGLE" if line<3 else "DOUBLE",
                "BS":{1,line+10,site+30},
                "BT":{1,line+10,site+40},
            })
    return sites


class DistanceMatchedLineTests(unittest.TestCase):
    def setUp(self):
        self.old=(m.N_SITES,m.N_LINES,m.N_SAME_UNORDERED,m.MIN_ANCHORS,m.MIN_LINES)
        m.N_SITES=15;m.N_LINES=5;m.N_SAME_UNORDERED=15
        m.MIN_ANCHORS=5;m.MIN_LINES=2

    def tearDown(self):
        m.N_SITES,m.N_LINES,m.N_SAME_UNORDERED,m.MIN_ANCHORS,m.MIN_LINES=self.old

    def test_field_line_is_distinct_from_geography_and_identical_site(self):
        source=make_source()
        rows,qa=m.match_candidates(source,20.,2.5,"geo_only")
        self.assertEqual(qa["source_pairs_same_line_undirected"],15)
        self.assertGreater(qa["source_anchor_sites_with_both_matched_reference_options"],0)
        for r in rows:
            anchor=source[r["anchor_index"]]
            same=source[r["within_index"]]
            other=source[r["other_index"]]
            self.assertNotEqual(r["anchor_index"],r["within_index"])
            self.assertNotEqual(r["anchor_index"],r["other_index"])
            self.assertEqual(anchor["section"],same["section"])
            self.assertNotEqual(anchor["section"],other["section"])
            self.assertLessEqual(r["source_pair_physical_distance_gap_km"],2.5)
            for typ in m.HABITATS:
                self.assertAlmostEqual(
                    r[typ]["matched_same_line_minus_other_line"],
                    r[typ]["same_line"]-r[typ]["different_line"])

    def test_strict_year_and_track_can_be_unidentified_without_chasing_support(self):
        source=make_source()
        rows,qa=m.match_candidates(source,20,2.5,"same_year_layout")
        self.assertLessEqual(len(rows),15)
        self.assertTrue(all(source[r["anchor_index"]]["year"]==
                            source[r["within_index"]]["year"]==
                            source[r["other_index"]]["year"] for r in rows))
        self.assertTrue(all(source[r["anchor_index"]]["layout"]==
                            source[r["within_index"]]["layout"]==
                            source[r["other_index"]]["layout"] for r in rows))

    def test_source_richness_matching_secondary_only(self):
        source=make_source()
        source[4]["BS"].update(range(100,200))
        source[4]["BT"].update(range(100,200))
        baseline,_=m.match_candidates(source,20,5,"geo_only")
        richer,_=m.match_candidates(source,20,5,"source_richness")
        self.assertLessEqual(len(richer),len(baseline))
        for r in richer:
            same=source[r["within_index"]]
            other=source[r["other_index"]]
            self.assertLessEqual(abs(len(m.source_community(same,"BS_and_BT"))-
                                     len(m.source_community(other,"BS_and_BT"))),5)

    def test_no_pair_support_does_not_invent_positive_effect(self):
        source=make_source()
        x,qa=m.match_candidates(source,0.00001,.00001,"geo_only")
        val=m.line_equal_aggregate(x,qa)
        self.assertFalse(val["identified"])
        self.assertEqual(val["source_anchor_sites"],0)
        self.assertNotIn("habitat_comparisons",val)

    def test_source_anchored_line_aggregates_and_section_loo(self):
        source=make_source()
        trials,qa=m.match_candidates(source,20,5,"geo_only")
        res=m.line_equal_aggregate(trials,qa)
        self.assertTrue(res["identified"])
        self.assertEqual(res["source_anchor_sites"],len(trials))
        self.assertEqual(res["source_lines_with_comparable_same_and_other_reference"],5)
        r=res["habitat_comparisons"]["BS_and_BT"]
        self.assertLessEqual(abs(r["equal_source_line_same_minus_different_line_similarity"]),1.)
        self.assertEqual(len(r["original_line_leave_one_out_range"]),2)
        # A Jaccard identity baseline should not become an assumed
        # taxonomic "same line" effect just because line labels exist.
        self.assertTrue(r["NOT_field_randomized_treatment_or_independent_anchor_pairs"])


if __name__=="__main__":
    unittest.main()
