"""Independent six-country circuit-vs-plants source QA without fake pairwise p-values."""
import sys
from pathlib import Path
from unittest.mock import patch
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_6country_circuit_macro_plants_v46 as study
from audit_gbif_equalarea_macroblocks_v29 import BLOCKS
from analyze_6country_50km_network_plants_v44 import BLOCKS_SIX,STRATEGIES,PERIODS


def network_source(road=False):
    return {"complete":True,"regions":{
        k:{"source":{"country":BLOCKS[k][0]},"routes":{s:{} for s in STRATEGIES}}
        if road else {"country":BLOCKS[k][0],"routes":{s:{} for s in STRATEGIES}}
        for k in BLOCKS_SIX
    }}


def circuits_by_one():
    out=[]
    for block in BLOCKS_SIX:
        out.append({"complete":True,"regions":{
            block:{"country":BLOCKS[block][0],
                   "actual_Julia_Circuitscape_executed":False,
                   "strategies":{s:{} for s in STRATEGIES}}
        }})
    return out


def fake_result(coeff=.2,regions=12,pairs=30):
    return {
        "n_two_period_eligible_GRIIS_50km_cells":regions,
        "periods":{
            period:{
                strategy:{
                    "n_dependent_cell_pair_rows":pairs,
                    "responses":{
                        response:{
                            study.Q_FIELD:coeff,
                            study.SHORTEST_FIELD:-.1,
                            study.EFF_FIELD:-.09,
                        } for response in study.RESPONSES
                    }
                } for strategy in STRATEGIES
            } for period in PERIODS
        }
    }


class SixCountryCircuitTests(unittest.TestCase):
    def test_source_frame_has_six_distinct_countries(self):
        self.assertEqual(len(BLOCKS_SIX),6)
        self.assertEqual(len({BLOCKS[b][0] for b in BLOCKS_SIX}),6)
        combined=study.merge_circuit_source_chunks(circuits_by_one())
        self.assertEqual(set(combined["regions"]),set(BLOCKS_SIX))

    def test_old_JP_source_present_but_not_an_independent_ecological_replicate(self):
        chunks=circuits_by_one()
        chunks[0]["regions"]["jp_tohoku"]={
            "country":"JP","actual_Julia_Circuitscape_executed":False,
            "strategies":{s:{} for s in STRATEGIES}
        }
        result=study.merge_circuit_source_chunks(chunks)
        self.assertEqual(result["predeclared_known_upstream_extra_JP_excluded"],["jp_tohoku"])
        self.assertNotIn("jp_tohoku",result["regions"])

    def test_missing_or_forged_Julia_metadata_rejected(self):
        partial=circuits_by_one()
        with self.assertRaisesRegex(ValueError,"incomplete"):
            study.merge_circuit_source_chunks(partial[:-1])
        bad=circuits_by_one()
        bad[0]["regions"][BLOCKS_SIX[0]]["actual_Julia_Circuitscape_executed"]=True
        with self.assertRaisesRegex(ValueError,"unsupported real-Julia"):
            study.merge_circuit_source_chunks(bad)

    def test_country_blocks_weighted_equally_not_as_120_independent_pairs(self):
        vals=[.2,.1,-.1,.05,.3,.25]
        src={k:fake_result(r,pairs=120) for k,r in zip(BLOCKS_SIX,vals)}
        q=study.blocks_summary(src,"2018_2025","nearest_in_cell",study.RESPONSES[0],study.Q_FIELD)
        self.assertEqual(q["n_independent_country_blocks_with_positive_network_and_flora_source"],6)
        self.assertEqual(q["n_blocks_in_railway_redundancy_hypothesis_positive_direction"],5)
        self.assertTrue(q["at_least_5_independent_country_blocks_source_support"])
        self.assertEqual(len(q["leave_one_country_out_medians"]),6)
        self.assertTrue(q["NO_pairwise_p_value_or_propagule_causality"])

    def test_low_source_support_stays_missing_not_null_effect(self):
        src={k:fake_result(.1) for k in BLOCKS_SIX}
        src["gb_yorkshire"]=fake_result(None,regions=7,pairs=18)
        agg=study.blocks_summary(src,"2018_2025","nearest_in_cell",study.RESPONSES[1],study.Q_FIELD)
        self.assertEqual(agg["n_independent_country_blocks_with_positive_network_and_flora_source"],5)
        self.assertFalse(agg["source_quality_and_identification"]["gb_yorkshire"]["block_support"])
        self.assertNotIn("gb_yorkshire",agg["regional_partial_rhos"])

    def test_complete_analysis_retains_original_two_period_and_anchor_sensitivity(self):
        flora={"complete":True,"regions":{k:{"country":BLOCKS[k][0]} for k in BLOCKS_SIX}}
        rail=network_source()
        road=network_source(road=True)
        circuits=study.merge_circuit_source_chunks(circuits_by_one())
        with patch.object(study,"analyze_city",side_effect=lambda block,*args,**kw:fake_result(.12)) as fake:
            out=study.run(flora,rail,road,circuits)
        self.assertEqual(fake.call_count,6)
        self.assertTrue(out["NO_railway_seed_dispersal_identified"])
        self.assertTrue(out["actual_Julia_Circuitscape_was_run_on_SYNTHETIC_TESTS_only"])
        self.assertEqual(len(out["country_equal_weight_results"]),2)
        self.assertEqual(out["country_equal_weight_results"]["2018_2025"]["nearest_in_cell"]["rarefied20_source_similarity"][study.Q_FIELD]["n_independent_country_blocks_with_positive_network_and_flora_source"],6)
        road["regions"].pop("ca_ontario")
        with self.assertRaisesRegex(ValueError,"geography"):
            study.run(flora,rail,road,circuits)


if __name__=="__main__":
    unittest.main()
