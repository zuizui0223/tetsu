"""v8.0 source-only 2022 FAF6 state-level directed freight OD algebra tests."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_FAF6_crop_OD_rail_vs_truck_network_v80 import (
    GROUPS,THOUSAND_SHORT_TON_THRESHOLDS,build_networks,mode_network,
    paired_origin_contrasts,audit_source_invariants)

def fake_rows():
    # 51 synthetic state category labels, 510+ original rows.
    for state in range(1,52):
        o=f"{state:02d}"
        next_state=f"{(state%51)+1:02d}"
        for mode in ("1","2"):
            for code in ("2","3","10","11","12"):
                yield {"dms_origst":o,"dms_destst":next_state,
                       "sctg2":code,"dms_mode":mode,"trade_type":"1",
                       "tons_2022":"5" if mode=="2" else "10"}
                # Real state OD matrix often also has within-state flows.
                yield {"dms_origst":o,"dms_destst":o,
                       "sctg2":code,"dms_mode":mode,"trade_type":"1",
                       "tons_2022":"2" if mode=="2" else "8"}

class StateODNetworkSource(unittest.TestCase):
    def test_original_state_code_scope_and_crop_resource_controls(self):
        self.assertEqual(set(GROUPS),{
            "grain_sctg02","other_agriculture_sctg03",
            "crop_sctg02_03","resource_sctg10_12"})
        self.assertEqual(THOUSAND_SHORT_TON_THRESHOLDS,(0,10,100))
    def test_directional_OD_not_physical_track_network(self):
        flows,states,trace=build_networks(fake_rows())
        self.assertEqual(len(states),51)
        self.assertEqual(trace["original_source_rows"],1020)
        rail=mode_network(flows,"crop_sctg02_03","RAIL",states)
        truck=mode_network(flows,"crop_sctg02_03","TRUCK",states)
        self.assertEqual(rail["n_positive_interstate_directed_OD_state_pairs"],51)
        self.assertEqual(rail["n_undirected_state_pairs_with_any_estimated_interstate_flow"],51)
        self.assertEqual(rail["n_origin_states_with_positive_interstate_flow"],51)
        self.assertAlmostEqual(rail["inter_state_2022_estimated_thousand_short_tons"],510)
        self.assertAlmostEqual(truck["intra_state_2022_estimated_thousand_short_tons"],816)
        self.assertGreater(rail["share_interstate_in_mode"],truck["share_interstate_in_mode"])
        self.assertTrue(rail["ORIGINAL_ECONOMIC_OD_TONNAGE_NOT_TRACK_ROUTING"])
    def test_source_reciprocity_ignores_same_state_loop(self):
        f={("crop_sctg02_03","RAIL","01","01"):80,
           ("crop_sctg02_03","RAIL","01","02"):5,
           ("crop_sctg02_03","RAIL","02","01"):3}
        m=mode_network(f,"crop_sctg02_03","RAIL",["01","02"])
        self.assertEqual(m["n_positive_interstate_directed_OD_state_pairs"],2)
        self.assertEqual(m["n_state_pairs_with_both_OD_directions_positive"],1)
        self.assertAlmostEqual(m["interstate_directional_weight_imbalance_0reciprocal_to1oneway"],.25)
        self.assertEqual(m["n_undirected_state_pairs_with_any_estimated_interstate_flow"],1)
    def test_matched_state_origin_partner_diversity_not_species_diversity(self):
        flows,states,_=build_networks(fake_rows())
        rail=mode_network(flows,"crop_sctg02_03","RAIL",states)
        truck=mode_network(flows,"crop_sctg02_03","TRUCK",states)
        r=paired_origin_contrasts(rail,truck)
        self.assertEqual(r["n_origin_states_with_both_RAIL_and_TRUCK_interstate_flow"],51)
        self.assertEqual(r["comparisons"]["effective_destination_count_shannon"][
            "n_state_origin_codes_with_both_modes_positive_interstate_outgoing"],51)
        self.assertTrue(r["NOT_independent_species_or_route_based_replication"])
    def test_v79_source_marginals_must_conserve_within_and_interstate(self):
        flows,states,trace=build_networks(fake_rows())
        output={"source_trace":trace,"mode_networks":{}}
        for group in ("crop_sctg02_03","resource_sctg10_12"):
            output["mode_networks"][group]={
                mode:mode_network(flows,group,mode,states) for mode in ("RAIL","TRUCK")}
        prev={"complete":True,"FAF_release":"FAF6.0 / 2022 BASE YEAR / STATE OD CSV",
              "data":{"n_original_official_FAF6_state_OD_commodity_mode_rows":1020,
              "n_distinct_state_code_categories_in_source":51,
              "national_trade_type_x_commodity_group_2022":{"domestic":{
                  "crop_codes_02_03":{"rail_thousand_short_tons":714,
                                      "truck_thousand_short_tons":1836},
                  "resource_codes_10_11_12":{"rail_thousand_short_tons":1071,
                                            "truck_thousand_short_tons":2754}}}}}
        # State rows: 51 x 2 crop categories x (5 interstate + 2 within) rail.
        audit_source_invariants(output,prev)
        self.assertTrue(output["source_trace"]["previous_verified_v79_total_and_modes_reconciled"])
        prev["data"]["national_trade_type_x_commodity_group_2022"]["domestic"][
            "crop_codes_02_03"]["rail_thousand_short_tons"]+=1
        with self.assertRaisesRegex(ValueError,"not conserved"):
            audit_source_invariants(output,prev)
    def test_source_sample_scope_cannot_be_inferred_from_tiny_toy(self):
        with self.assertRaisesRegex(ValueError,"population"):
            build_networks([{"dms_origst":"01","dms_destst":"02",
                             "trade_type":"1","sctg2":"2","dms_mode":"2",
                             "tons_2022":"1.0"}])
if __name__=="__main__":unittest.main()
