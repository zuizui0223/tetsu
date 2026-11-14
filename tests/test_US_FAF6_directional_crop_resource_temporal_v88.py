"""v8.8: synthetic math and data-integrity gates for 2022 directional
CROP vs RESOURCE RAIL cargo against pre/post GBIF recorded species.
No actual biological colonization or seed contamination inference.
"""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_US_FAF6_directional_crop_resource_temporal_v88 import (
    COMMODITIES,OUTCOMES,PREDICTORS,SOURCE_COHORTS,
    augment_directed_rows,source_cohort,cohort_summary,
    conditional_partial_rank,node_deletion)

def test_row(o="01",d="02"):
    return {
      "o":o,"d":d,
      "x_rail":math.log1p(11),"x_truck":math.log1p(12),
      "rail_2022_estimated_thousand_short_tons":11,
      "truck_2022_estimated_thousand_short_tons":12,
      "responses":{
         OUTCOMES[0]:.06,OUTCOMES[1]:.3,
      },
      "covars":{"log1p_truck_OD":math.log1p(12)}
    }

def sources():
    return {
        "crop_sctg02_03":{
           "RAIL":{("01","02"):11},
           "TRUCK":{("01","02"):12}},
        "resource_sctg10_12":{
           "RAIL":{("01","02"):9,("02","01"):2},
           "TRUCK":{("01","02"):7}},
    }

class CargoSpecificityTemporalQA(unittest.TestCase):
    def test_commodity_labels_do_not_claim_sterile_negative_control(self):
        self.assertEqual(COMMODITIES,("crop_sctg02_03","resource_sctg10_12"))
        self.assertEqual(PREDICTORS,("crop_rail","resource_rail"))
        self.assertIn("raw_source_linked_first_record_fraction",OUTCOMES)
    def test_origin_destination_2022_rail_direction_is_not_symmetrized(self):
        data=augment_directed_rows([test_row()],sources())
        self.assertEqual(len(data),1)
        r=data[0]
        self.assertEqual(r["crop_rail"],11)
        self.assertEqual(r["resource_rail"],9)
        self.assertEqual(r["resource_truck"],7)
        self.assertAlmostEqual(r["covars"]["log1p_resource_origin_rail_out_strength"],math.log1p(9))
        self.assertAlmostEqual(r["covars"]["log1p_resource_destination_rail_in_strength"],math.log1p(9))
    def test_official_crop_OD_edited_without_provenance_is_rejected(self):
        row=test_row()
        row["rail_2022_estimated_thousand_short_tons"]=20
        with self.assertRaisesRegex(ValueError,"drift"):
            augment_directed_rows([row],sources())
    def test_cargo_positive_cohorts_fixed_without_plant_response_sign(self):
        a=augment_directed_rows([test_row()],sources())
        b=dict(a[0],o="02",d="01",crop_rail=0,crop_truck=0,resource_rail=2)
        x=[*a,b]
        self.assertEqual(len(source_cohort(x,SOURCE_COHORTS[0])),2)
        self.assertEqual(len(source_cohort(x,SOURCE_COHORTS[1])),1)
        self.assertEqual(len(source_cohort(x,SOURCE_COHORTS[2])),1)
        self.assertEqual(len(source_cohort(x,SOURCE_COHORTS[3])),1)
        with self.assertRaisesRegex(ValueError,"independent"):
            source_cohort(x,"only_positive_botanical_response")
    def test_insufficient_independent_source_nodes_produces_null_rho(self):
        rows=augment_directed_rows([test_row()],sources())
        a=cohort_summary(rows)
        self.assertIsNone(a["botanical_record_source_results"][OUTCOMES[0]][
            "crop_RAIL_full_adjusted_rho"])
        self.assertFalse(a["source_quality_gate_at_least_20_origin_and_dest_states"])
        self.assertIsNone(conditional_partial_rank(rows,OUTCOMES[0],"crop_rail"))
        jack=node_deletion(rows,"crop_rail")
        self.assertEqual(jack["origin"]["n_estimable"],0)
        self.assertEqual(jack["destination"]["n_estimable"],0)
    def test_unsupported_response_or_cargo_cannot_create_estimated_effect(self):
        rows=augment_directed_rows([test_row()],sources())
        with self.assertRaisesRegex(ValueError,"Unpredeclared"):
            conditional_partial_rank(rows,"species_colonization","crop_rail")
        with self.assertRaisesRegex(ValueError,"Unpredeclared"):
            conditional_partial_rank(rows,OUTCOMES[0],"seed_laden_trains")
if __name__=="__main__":unittest.main()
