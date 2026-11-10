"""Synthetic country-block integrity tests for comparative RAIL/major ROAD circuit ecology."""
import sys
from pathlib import Path
import unittest
from itertools import combinations
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_rail_vs_major_road_circuit_plants_v42 import run,source_pairs,correlate_block,STRATEGIES,PERIODS


def case():
    flora={"complete":True,"regions":{}}
    rail={"complete":True,"regions":{}}
    road={"complete":True,"regions":{}}
    circuit={"complete":True,"regions":{}}
    us={"complete":True,"regions":{}}
    fr={"complete":True,"regions":{}}
    for block,country in [("us_midwest","US"),("fr_burgundy","FR")]:
        cells={}
        for i in range(16):
            name=f"{block}_x{i%4}_y{i//4}"
            counts={str(100+k):20 for k in range(3)}
            counts[str(200+i)]=20
            counts[str(300+i%4)]=20
            row={
                "full_facet_source_pagination_completed":True,
                "alien_observed_minimum_5_species_and_20_records":True,
                "observed_alien_specieskey_counts":counts,
                "n_indexed_GRIIS_introduced_plant_records":100,
                "indexed_all_plantae_source_records":250+i*4,
                "prior_2026_GBIF_indexed_source_record_count":250+i*4
            }
            cells[name]={"source_pre_outcome_both_period_recording_eligible":True,
                          "periods":{period:dict(row) for period in PERIODS}}
        flora["regions"][block]={"country":country,"cells":cells}
        rr=[];rs=[];cr=[];dr=[]
        for i,(a,b) in enumerate(combinations(sorted(cells),2)):
            geographic=50*(1+i%4)
            road_len=geographic*(1.15+.04*(i%10))
            rail_len=geographic*(1.45+.02*(i%11))
            rr.append({"a":a,"b":b,
                       "cell_centroid_snap_augmented_rail_path_km":rail_len,
                       "source_cell_centroid_geodesic_km":geographic})
            rs.append({"a":a,"b":b,"centroid_snap_augmented_road_path_km":road_len})
            cr.append({"a":a,"b":b,
                       "independent_path_redundancy_1minusRoverD":.3+.01*(i%13)})
            dr.append({"a":a,"b":b,
                       "source_backbone_Q_one_minus_Reff_over_shortest":.4+.009*(i%17)})
        rail["regions"][block]={"routes":{s:{"pair_routes":[dict(x) for x in rr]} for s in STRATEGIES}}
        road["regions"][block]={"routes":{s:{"pair_paths":[dict(x) for x in rs]} for s in STRATEGIES}}
        circuit["regions"][block]={"actual_Julia_Circuitscape_executed":False,
                                   "strategies":{s:{"pairwise_source_diagnostics":[dict(x) for x in cr]}
                                                 for s in STRATEGIES}}
        dest=us if country=="US" else fr
        dest["regions"][block]={"country":country,"road_backbone_not_full_network":True,
                               "strategies":{s:{"pair_source_diagnostics":[dict(x) for x in dr]}
                                             for s in STRATEGIES}}
    return flora,rail,road,circuit,us,fr


class RoadRailCircuitPlantTests(unittest.TestCase):
    def test_all_two_macro_blocks_and_both_periods(self):
        sources=case()
        result=run(*sources)
        self.assertEqual(len(result["regions"]),2)
        self.assertTrue(result["rail_source_class_vs_road_backbone_class_not_equivalent"])
        self.assertTrue(result["only_two_independent_geographic_blocks"])
        self.assertTrue(result["not_actual_train_seed_transport"])
        for block in result["regions"]:
            for period in PERIODS:
                for strategy in STRATEGIES:
                    x=result["regions"][block]["periods"][period][strategy]
                    self.assertEqual(x["n_dependent_source_pairs_with_both_circuits"],120)
                    self.assertEqual(len(x["Qroad_major_backbone_quartiles"]),3)
                    self.assertTrue(x["NO_pairwise_p_values"])
                    self.assertIn("rail_Q_given_road_backbone_Q_as_extra_covariate",
                                  x["partial_rank_association"])

    def test_unreported_road_source_not_invented_zero(self):
        data=list(case())
        del data[-1]["regions"]["fr_burgundy"]
        with self.assertRaisesRegex(ValueError,"Independent US and French"):
            run(*data)

    def test_mismatch_duplicate_pair_fails_without_ecological_claim(self):
        data=list(case())
        x=data[-2]["regions"]["us_midwest"]["strategies"]["nearest_in_cell"]["pair_source_diagnostics"]
        x.pop()
        with self.assertRaisesRegex(ValueError,"Incomplete 16-cell ROAD"):
            run(*data)

    def test_circuits_source_must_declare_major_backbone_only(self):
        data=list(case())
        data[-2]["regions"]["us_midwest"]["road_backbone_not_full_network"]=False
        with self.assertRaisesRegex(ValueError,"class/fidelity"):
            run(*data)

    def test_source_frame_incomplete_gated(self):
        data=list(case())
        data[0]["complete"]=False
        with self.assertRaisesRegex(ValueError,"indispensable"):
            run(*data)


if __name__=="__main__":
    unittest.main()
