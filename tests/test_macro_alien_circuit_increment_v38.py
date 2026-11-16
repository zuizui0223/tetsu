"""Only synthetic checks for optional circuit incremental ecological association."""
import sys
from pathlib import Path
from itertools import combinations
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_macro_alien_circuit_increment_v38 import run,analyze_city,STRATEGIES


def fake_sources():
    flora={"complete":True,"regions":{}}
    rail={"complete":True,"regions":{}}
    road={"complete":True,"regions":{}}
    circuit={"complete":True,"regions":{}}
    for block,country in (("us_midwest","US"),("fr_burgundy","FR")):
        ids=[f"{block}_x{i%4}_y{i//4}" for i in range(16)]
        cells={}
        for i,cid in enumerate(ids):
            species={str(100+j):20 for j in range(3)}
            species[str(200+i)]=20
            species[str(300+i%4)]=20
            row={"full_facet_source_pagination_completed":True,
                 "alien_observed_minimum_5_species_and_20_records":True,
                 "observed_alien_specieskey_counts":species,
                 "n_indexed_GRIIS_introduced_plant_records":100,
                 "indexed_all_plantae_source_records":200+i*4,
                 "prior_2026_GBIF_indexed_source_record_count":200+i*4}
            cells[cid]={
                "source_pre_outcome_both_period_recording_eligible":True,
                "periods":{"2010_2017":dict(row),"2018_2025":dict(row)}
            }
        flora["regions"][block]={"country":country,"cells":cells}
        railrows=[];roadrows=[];circuitrows=[]
        for i,(a,b) in enumerate(combinations(sorted(ids),2)):
            dist=50.*(1+i%5)
            railkm=dist*(1.35+.04*(i%13))
            roadkm=dist*(1.10+.025*(i%11))
            ratio=.5+.025*(i%12)
            railrows.append({
                "a":a,"b":b,
                "rail_to_geodesic_ratio":railkm/dist,
                "cell_centroid_snap_augmented_rail_path_km":railkm,
                "source_cell_centroid_geodesic_km":dist,
            })
            roadrows.append({
                "a":a,"b":b,"road_to_geodesic_ratio":roadkm/dist,
                "centroid_snap_augmented_road_path_km":roadkm,
            })
            circuitrows.append({
                "a":a,"b":b,
                "centroid_augmented_R_eff_divided_by_Dijkstra":ratio,
                "independent_path_redundancy_1minusRoverD":1-ratio,
                "physical_circuit_reachable":True,
            })
        rail["regions"][block]={"routes":{s:{"pair_routes":[dict(x) for x in railrows]} for s in STRATEGIES}}
        road["regions"][block]={"routes":{s:{"pair_paths":[dict(x) for x in roadrows]} for s in STRATEGIES}}
        circuit["regions"][block]={
            "actual_Julia_Circuitscape_executed":False,
            "strategies":{s:{"pairwise_source_diagnostics":[dict(x) for x in circuitrows]}
                          for s in STRATEGIES}
        }
    return flora,rail,road,circuit


class CircuitEcologicalSourceTests(unittest.TestCase):
    def test_rail_vs_road_and_circuit_increment_two_periods_two_anchors(self):
        data=fake_sources()
        out=run(*data)
        self.assertTrue(out["actual_Julia_Circuitscape_was_not_run"])
        self.assertTrue(out["only_two_independent_country_macro_blocks"])
        self.assertEqual(set(out["regions"]),{"us_midwest","fr_burgundy"})
        for block in out["regions"]:
            for period in ("2010_2017","2018_2025"):
                for strategy in STRATEGIES:
                    row=out["regions"][block]["periods"][period][strategy]
                    self.assertEqual(row["n_dependent_cell_pair_rows"],120)
                    self.assertEqual(len(row["Q_alternative_path_fraction_quartiles"]),3)
                    self.assertNotIn("significance_p_value",row)
                    self.assertIn("partial_rank_extra_path_redundancy_Q_at_fixed_shortest_rail",
                                  row["responses"]["rarefied20_source_similarity"])

    def test_prior_plant_source_quality_is_required(self):
        a,b,c,d=fake_sources()
        a["regions"]["us_midwest"]["cells"]["us_midwest_x0_y0"]["source_pre_outcome_both_period_recording_eligible"]=False
        out=run(a,b,c,d)
        self.assertEqual(out["regions"]["us_midwest"]["n_two_period_eligible_GRIIS_50km_cells"],15)
        self.assertEqual(out["regions"]["us_midwest"]["periods"]["2018_2025"]["nearest_in_cell"]["n_dependent_cell_pair_rows"],105)

    def test_missing_circuit_pair_is_error_not_imputed_zero(self):
        a,b,c,d=fake_sources()
        d["regions"]["us_midwest"]["strategies"]["nearest_in_cell"]["pairwise_source_diagnostics"].pop()
        with self.assertRaisesRegex(ValueError,"Incomplete pairwise"):
            run(a,b,c,d)

    def test_non_source_circuitscape_claim_blocked(self):
        a,b,c,d=fake_sources()
        d["regions"]["fr_burgundy"]["actual_Julia_Circuitscape_executed"]=True
        with self.assertRaisesRegex(ValueError,"Julia Circuitscape"):
            run(a,b,c,d)

    def test_empty_real_prerequisites_never_become_ecological_significance(self):
        a,b,c,d=fake_sources();d["complete"]=False
        with self.assertRaises(ValueError):
            run(a,b,c,d)
        a,b,c,d=fake_sources();c["regions"].pop("fr_burgundy")
        with self.assertRaises(ValueError):
            run(a,b,c,d)


if __name__=="__main__":
    unittest.main()
