"""v27.8 synthetic scientific QA: route-windowed original railway BS/BT
checkerboard swaps conserve site richness and section-specific taxon frequency."""
import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Slovenia_rail_nearby_route_distance_fixed_frequency_null_v278 as m
from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import source_freq


def rail_source(n_lines=4,n_site=5):
    panel={}
    for line in ("10","20","34","80")[:n_lines]:
        group=[]
        for i in range(n_site):
            km=(0,1,3,6,9,12,15)[i]
            group.append({
                "source_id":f"{line}_{i}","section":line,"routekm":float(km),
                "year":2022 if i<4 else 2023,
                "lat":46.0,"lon":14.+km/200,
                "BS":{1,i+10},"BT":{1,i+10,i+20},
            })
        panel[line]=group
    return panel


class LocalWindowFrequencyNullTests(unittest.TestCase):
    def test_edges_never_cross_original_line_or_route_radius(self):
        original=rail_source(1)["10"]
        edges,qa=m.eligible_edges(original,5,"all")
        self.assertGreater(len(edges),0)
        self.assertLess(len(edges),10)
        self.assertTrue(all(abs(original[i]["routekm"]-original[j]["routekm"])<=5 for i,j in edges))
        same,qa2=m.eligible_edges(original,5,"same_survey_year")
        self.assertLessEqual(len(same),len(edges))
        self.assertGreater(qa2["source_year_mismatch_edges_excluded"],0)
        with self.assertRaises(ValueError):
            m.eligible_edges(original,5,"other_region")

    def test_exact_original_row_and_species_frequency_margins(self):
        original=rail_source(1)["10"]
        rows=[set(x["BT"]) for x in original]
        sizes=[len(x) for x in rows]
        occ=source_freq(rows)
        edges,_=m.eligible_edges(original,5,"all")
        n=m.propose(rows,edges,random.Random(9),300)
        self.assertGreater(n,0)
        self.assertEqual(sizes,[len(x) for x in rows])
        self.assertEqual(occ,source_freq(rows))

    def test_dedicated_unexchangeable_reference_is_not_zero_effect(self):
        source=rail_source(1)["10"]
        for x in source:x["BT"]={1,2}
        out=m.original_rail_section_null(source,5,"all",draws=29,burnin=50,spacing=10)
        self.assertFalse(out["informative"])
        self.assertEqual(out["reason"],"SOURCE_INCIDENCE_MATRIX_NOT_EXCHANGEABLE_AT_ROUTE_RADIUS")
        self.assertEqual(out["observed"]["beta_jac"],.75)
        s2=rail_source(1)["10"]
        for i,x in enumerate(s2):x["routekm"]=float(i*30)
        out=m.original_rail_section_null(s2,5,"all",draws=29,burnin=50,spacing=10)
        self.assertFalse(out["informative"])
        self.assertEqual(out["reason"],"NO_DISTANCE_ELIGIBLE_NEARBY_REAL_SITE_PAIR")

    def test_common_same_source_section_panel_across_all_radii(self):
        base=rail_source()
        panels={"all":{}}
        for radius in m.RADII:
            k="ALL" if math.isinf(radius) else f"{radius:g}km"
            aggregate,internal=m.output_panel(base,mode="all",radius=radius,
                                               draws=29,burnin=100,spacing=15,seed=100)
            panels["all"][k]={"source_aggregate":aggregate,"source_internal":internal}
            self.assertTrue(aggregate["identified"])
            self.assertEqual(aggregate["n_informative_original_sections"],4)
            self.assertEqual(aggregate["source_original_sites_in_informative_sections"],20)
        common=m.intersection_stats(panels,"all",draws=29)
        self.assertTrue(common["identified"])
        self.assertEqual(common["original_sections_in_all_radii"],4)
        self.assertEqual(len(common["results"]),5)
        observed=[x["original_site_source_beta_jac"] for x in common["results"].values()]
        self.assertTrue(all(abs(v-observed[0])<1e-12 for v in observed))
        self.assertTrue(all(q["original_section_count_local_coupling_positive"]<=4 for q in common["results"].values()))

    def test_group_minimal_sample_eligibility_not_fabricated(self):
        s=rail_source(2)
        aggregate,internal=m.output_panel(s,mode="all",radius=5,
                                          draws=29,burnin=50,spacing=10)
        self.assertFalse(aggregate["identified"])
        self.assertEqual(aggregate["n_informative_original_sections"],2)
        self.assertEqual(aggregate["source_original_sites_in_informative_sections"],10)
        self.assertTrue(all(q["informative"] for q in internal.values()))


if __name__=="__main__":
    unittest.main()
