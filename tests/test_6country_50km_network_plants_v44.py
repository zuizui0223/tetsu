"""Offline independent-country macro source aggregation and missingness gates."""
import sys
from pathlib import Path
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_6country_50km_network_plants_v44 as study


def network_chunks(road=False):
    chunks=[]
    for block in study.BLOCKS_SIX:
        row={"routes":{"nearest_in_cell":{},"nearest_largest_component_in_cell":{}}}
        country=study.BLOCKS[block][0]
        if road:
            row["source"]={"country":country}
        else:
            row["country"]=country
        chunks.append({"complete":True,"regions":{block:row},"errors":{}})
    return chunks


def synthetic_region(rho,n=30,cells=10):
    return {
        "n_source_qualified_two_period_plant_cells":cells,
        "periods":{
            period:{
                anchor:{
                    "n_plant_pair_routes_with_both_network_lengths":n,
                    "responses":{
                        response:{"descriptive_partial_rank_correlation":rho}
                        for response in study.RESPONSES
                    }
                } for anchor in study.STRATEGIES
            } for period in study.PERIODS
        }
    }


class SixCountryMacroTests(unittest.TestCase):
    def test_country_registry_frozen_and_geographic_independent(self):
        self.assertEqual(len(study.BLOCKS_SIX),6)
        self.assertEqual(len(set(study.BLOCKS_SIX)),6)
        self.assertEqual(len({study.BLOCKS[x][0] for x in study.BLOCKS_SIX}),6)
        self.assertEqual(set(study.INDEPENDENT_ADDITION),
                         {"ca_ontario","gb_yorkshire","de_saxony","au_southeast"})

    def test_merge_source_true_source_chunks_without_replicates(self):
        rail=study.combine_source_regions(network_chunks())
        road=study.combine_source_regions(network_chunks(road=True))
        self.assertEqual(rail["n_geographically_independent_blocks"],6)
        self.assertEqual(road["n_true_source_region_chunks"],6)
        self.assertEqual(set(rail["regions"]),set(study.BLOCKS_SIX))

    def test_missing_or_duplicate_survey_source_fails_closed(self):
        chunks=network_chunks()
        with self.assertRaisesRegex(ValueError,"incomplete"):
            study.combine_source_regions(chunks[:-1])
        with self.assertRaisesRegex(ValueError,"twice"):
            study.combine_source_regions(chunks+[chunks[0]])
        bad=network_chunks()
        bad[0]["regions"][study.BLOCKS_SIX[0]]["country"]="KR"
        with self.assertRaisesRegex(ValueError,"country"):
            study.combine_source_regions(bad)

    def test_country_equal_weight_direction_and_leave_one_out(self):
        values=[-.3,-.2,-.1,.05,-.05,-.02]
        regions={block:synthetic_region(v) for block,v in zip(study.BLOCKS_SIX,values)}
        d=study.aggregate(regions,"2018_2025","nearest_in_cell",study.RESPONSES[0])
        self.assertEqual(d["n_independent_regions_with_qualified_network_and_flora"],6)
        self.assertEqual(d["n_blocks_in_hypothesis_compatible_negative_direction"],5)
        self.assertEqual(d["n_blocks_in_opposite_positive_direction"],1)
        self.assertTrue(d["minimum_five_independent_regions_for_descriptive_generalization"])
        self.assertEqual(len(d["leave_one_independent_country_out_medians"]),6)
        self.assertNotIn("p_value",d)
        self.assertTrue(d["has_no_hypothesis_significance_p_value"])

    def test_partial_source_coverage_does_not_get_zero_effect(self):
        vals={block:synthetic_region(-.1) for block in study.BLOCKS_SIX}
        vals["ca_ontario"]=synthetic_region(None,n=12,cells=7)
        q=study.aggregate(vals,"2010_2017","nearest_largest_component_in_cell",study.RESPONSES[1])
        self.assertEqual(q["n_independent_regions_with_qualified_network_and_flora"],5)
        self.assertFalse(q["pre_outcome_source_support_audit"]["ca_ontario"]["regional_source_support"])
        self.assertNotIn("ca_ontario",q["regional_rhos"])

    def test_run_checks_country_footprints_and_no_pairwise_significance(self):
        flora={"complete":True,"source_errors":{},
               "regions":{block:{"country":study.BLOCKS[block][0]}
                          for block in study.BLOCKS_SIX}}
        rail=study.combine_source_regions(network_chunks())
        road=study.combine_source_regions(network_chunks(road=True))
        with patch.object(study,"analyze_one_block",side_effect=lambda name,*args:synthetic_region(-.1)) as mocked:
            x=study.analyze(flora,rail,road)
        self.assertEqual(mocked.call_count,6)
        self.assertTrue(x["observational_and_dependent_pairwise_no_p_values"])
        self.assertEqual(x["country_equal_weight_aggregate"]["2018_2025"]["nearest_in_cell"][study.RESPONSES[0]]["n_independent_regions_with_qualified_network_and_flora"],6)
        flora["regions"]["ca_ontario"]["country"]="FR"
        with patch.object(study,"analyze_one_block",side_effect=lambda name,*args:synthetic_region(-.1)):
            with self.assertRaisesRegex(ValueError,"country mismatch"):
                study.analyze(flora,rail,road)


if __name__=="__main__":
    unittest.main()
