"""v27.9 synthetic QA: original identical source railway sites evaluated under
route-chainage vs straight-line 5km source taxon-frequency matrix nulls."""
from __future__ import annotations
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Slovenia_rail_route_vs_geo_distance_source_null_v279 as m
from audit_Slovenia_rail_nearby_route_distance_fixed_frequency_null_v278 import eligible_edges


def make_source(nlines=4):
    panel={}
    for line in ("10","20","34","80")[:nlines]:
        sites=[]
        for i,km in enumerate((0.,1.,3.,6.,9.)):
            sites.append({
                "source_id":f"{line}_{i}","section":line,
                "routekm":km,"lat":46.,"lon":14.+km/200.,
                "year":2022 if i<4 else 2023,
                "BS":{1,10+i},"BT":{1,10+i,20+i}
            })
        panel[line]=sites
    return panel


class RouteVsGeodesicReferenceTests(unittest.TestCase):
    def test_same_original_site_set_and_geometry_candidate_edges(self):
        source=make_source(1)["10"]
        geo=m.geodesic_edges(source,5,"all")
        route,_=eligible_edges(source,5,"all")
        self.assertGreater(len(geo),len(route))
        self.assertEqual(set(m.endpoints(geo)),set(range(5)))
        self.assertEqual(set(m.endpoints(route)),set(range(5)))
        self.assertTrue(all(source[i]["section"]==source[j]["section"] for i,j in geo))

    def test_exact_row_and_species_occ_margins_both_reference_systems(self):
        orig=make_source(1)["10"]
        q=m.section_pair(orig,5,"all",draws=29,burnin=120,spacing=12,seed=99)
        self.assertTrue(q["identified"])
        self.assertEqual(q["support"]["same_source_sites_with_BOTH_real_distance_neighborhoods"],5)
        self.assertAlmostEqual(
            q["source_route_limited_reference"]["original_source_observed_jaccard_beta"],
            q["source_geodesic_limited_reference"]["original_source_observed_jaccard_beta"])
        self.assertTrue(q["source_route_limited_reference"]["original_source_site_and_taxon_frequencies_fixed_exact"])
        self.assertTrue(q["source_geodesic_limited_reference"]["original_source_site_and_taxon_frequencies_fixed_exact"])

    def test_source_neighborhood_does_not_fake_unmatched_site(self):
        source=make_source(1)["10"]
        source[-1]["routekm"]=100.
        source[-1]["lon"]=15.
        output=m.section_pair(source,5,"all",draws=29,burnin=120,spacing=12,seed=99)
        self.assertEqual(output["support"]["same_source_sites_with_BOTH_real_distance_neighborhoods"],4)
        self.assertEqual(output["support"]["original_sites_without_both_type_neighborhoods"],1)
        self.assertEqual(output["source_route_limited_reference"]["source_test_site_count"],4)
        self.assertEqual(output["source_geodesic_limited_reference"]["source_test_site_count"],4)

    def test_full_original_section_equal_hierarchical_support(self):
        out=m.aggregate_panel(make_source(),5,"all",draws=29,burnin=100,spacing=10,seed=19)
        self.assertTrue(out["identified"])
        self.assertEqual(out["source_original_railway_sections_with_both_nulls_informative"],4)
        self.assertEqual(out["source_original_actual_matched_neighborhood_sites"],20)
        self.assertEqual(out["original_sections_positive_route_excess"],4)
        self.assertEqual(out["original_sections_positive_geodesic_excess"],4)
        self.assertIsNotNone(out["route_minus_geo_leave_original_section_out_range"])

    def test_low_actual_source_support_fails_inference_gate(self):
        out=m.aggregate_panel(make_source(2),5,"all",draws=29,burnin=100,spacing=10,seed=19)
        self.assertFalse(out["identified"])
        self.assertEqual(out["source_original_railway_sections_with_both_nulls_informative"],2)

    def test_no_artificial_region_hopping_and_survey_year_gate(self):
        source=make_source(1)["10"]
        edge=m.geodesic_edges(source,5,"same_survey_year")
        self.assertTrue(all(source[i]["year"]==source[j]["year"] for i,j in edge))
        with self.assertRaises(ValueError):m.geodesic_edges(source,-1,"all")


if __name__=="__main__":
    unittest.main()
