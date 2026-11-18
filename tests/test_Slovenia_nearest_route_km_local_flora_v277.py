"""v27.7 synthetic checks: source route km nearest OTHER rail line location,
matching original survey year/track and preserving section aggregation."""
from __future__ import annotations
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Slovenia_rail_route_nearest_plant_microhabitat_v277 as m


def locations(line, entries):
    out=[]
    for i,(km,year,layout) in enumerate(entries):
        out.append({"source_id":f"{line}_{i}","section":line,
             "class":"main_G" if line in ("10","20") else "regional_R",
             "year":year,"layout":layout,
             "routekm":float(km),
             "lat":46.0,"lon":14.+km/200.,
             "BS":{1,i+10},"BT":{1,i+10,i+100}})
    return out


class RouteNearestFloraTests(unittest.TestCase):
    def test_route_km_nearest_vs_farther_other_location(self):
        q=locations("10",[(0,2022,"SINGLE"),(2,2022,"SINGLE"),
                          (100,2023,"DOUBLE")])
        # At km=0, closest 2km neighbor must be selected regardless of
        # plant species overlap; source route km is an observation-independent criterion.
        q[0]["BS"]={1,3}
        q[0]["BT"]={1,3}
        q[1]["BT"]={1,4}
        q[2]["BT"]={1,3}
        pairs=m.source_neighbors(q,math.inf,"all")
        source0=next(x for x in pairs if x["source_id"]=="10_0")
        self.assertEqual(source0["route_gap_km"],2.)
        self.assertAlmostEqual(source0["same_similarity"],1.)
        self.assertAlmostEqual(source0["neighbor_similarity"],1/3)
        self.assertAlmostEqual(source0["same_minus_nearest_other"],2/3)
        # If taxon similarity had been used to pick a different location,
        # the far-away km100 plant community would have won, which is forbidden.

    def test_source_radius_and_original_survey_year_layout_gates(self):
        q=locations("10",[(0,2022,"SINGLE"),(3,2023,"DOUBLE"),
                          (5,2022,"SINGLE")])
        near=m.source_neighbors(q,4,"all")
        self.assertEqual(len(near),2)  # site at 5km lacks another within 4km? 5-3=2, yes 3 sites
        # Strict same-year/site-layout must choose actual 5km site, never
        # use source-year 2023 observation at 3km for the 2022 site at 0.
        within=m.source_neighbors(q,5,"same_year_and_layout")
        self.assertEqual({x["source_id"] for x in within},{"10_0","10_2"})
        self.assertTrue(all(x["route_gap_km"]==5 for x in within))
        with self.assertRaises(ValueError):
            m.source_neighbors(q,5,"not-an-original-source-condition")

    def test_line_equal_section_mean_not_plant_dyad_pseudoreplicates(self):
        blocks={
            "10":locations("10",[(0,2022,"SINGLE"),(2,2022,"SINGLE"),(7,2022,"SINGLE")]),
            "20":locations("20",[(0,2022,"SINGLE"),(1,2022,"SINGLE"),(3,2022,"SINGLE")]),
            "34":locations("34",[(0,2022,"SINGLE"),(4,2022,"SINGLE"),(12,2022,"SINGLE")]),
            "80":locations("80",[(0,2023,"SINGLE"),(4,2023,"SINGLE"),(20,2023,"SINGLE")]),
        }
        out=m.summarise(blocks,"all",math.inf)
        self.assertTrue(out["identified"])
        self.assertEqual(out["source_real_rail_sections"],4)
        self.assertEqual(out["source_eligible_original_spatial_sites"],12)
        self.assertEqual(out["source_original_G_main_sections"],2)
        self.assertEqual(out["source_original_R_regional_sections"],2)
        self.assertEqual(len(out["source_section_results_aggregate_only"]),4)
        self.assertEqual(out["source_original_line_sections_with_positive_same_nearest_gap"],4)
        within2=m.summarise(blocks,"all",2)
        self.assertFalse(within2["identified"])
        self.assertLess(within2["source_eligible_original_spatial_sites"],12)

    def test_route_distance_bin_and_actual_original_pair_counts(self):
        q=locations("80",[(0,2022,"SINGLE"),(3,2022,"SINGLE"),
                          (9,2022,"SINGLE"),(45,2023,"SINGLE")])
        x=m.section_distance_profile(q)
        self.assertEqual(x["n_different_original_location_directed_source_pairs"],12)
        self.assertEqual(sum(v["n_dependent_directed_pairs_NOT_independent"] for v in x["bins"].values()),12)
        self.assertEqual(m.bin_for_distance(1),"0-5")
        self.assertEqual(m.bin_for_distance(12),"5-20")
        self.assertEqual(m.bin_for_distance(30),"20-50")
        self.assertEqual(m.bin_for_distance(60),"50+")
        with self.assertRaises(ValueError):m.bin_for_distance(float("nan"))

    def test_no_nearest_site_cannot_create_fake_0km_match(self):
        q=locations("80",[(0,2022,"SINGLE"),(30,2023,"DOUBLE")])
        self.assertEqual(m.source_neighbors(q,5,"same_year_and_layout"),[])
        self.assertEqual(m.source_neighbors(q,20,"all"),[])
        out=m.summarise({"80":q},"all",math.inf)
        self.assertFalse(out["identified"])
        self.assertEqual(out["source_real_rail_sections"],1)
        self.assertEqual(out["source_eligible_original_spatial_sites"],2)

if __name__=="__main__":
    unittest.main()
