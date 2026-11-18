"""v27.5 synthetic QA for real 9 source physical traffic intervals and original line-specified freight permutations."""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Slovenia_historical_train_frequency_local_species_coupling_v275 as m
from validate_SI_actual_2020_freight_vegetation_beta_holdout_v233 import exact_within_line_labels


def nine_original_physical_segments():
    # Original interval counts by line are synthetic but retain exactly
    # the v23.3 supported 4!×2!=48 within-national-line arrangements.
    line_group=[("10",2),("20",1),("34",1),("50",4),("80",1)]
    out=[]
    for lin,n in line_group:
        for j in range(n):
            freight=120*(j+1)+int(lin)*3
            passenger=450+27*j+int(lin)//2+(17 if j%2 else 0)
            out.append({
                "original_line":lin,"class":m.official_class(lin),
                "goods_2020":freight,"passenger_2020":passenger,
                "mean_local_excess":.07+.1*math.log1p(freight)
                                  +.02*(j%3),
                "source_plant_site_count":2
            })
    return out


class ActualHistoricalTrafficSourceTests(unittest.TestCase):
    def test_nine_segments_five_lines_exact_48_assignments(self):
        data=nine_original_physical_segments()
        self.assertEqual(len(data),9)
        meta={i:{"national_line":x["original_line"],"freight_2020":x["goods_2020"]}
              for i,x in enumerate(data)}
        perms=list(exact_within_line_labels(meta))
        self.assertEqual(perms[0],(48,None))
        self.assertEqual(len(perms)-1,48)
        self.assertTrue(all(set(assignment)==set(range(9)) for _,assignment in perms[1:]))

    def test_within_line_positive_freight_coupling_and_passenger_control(self):
        data=nine_original_physical_segments()
        out=m.centred_slope(data)
        self.assertTrue(out["identified"])
        self.assertGreater(out["within_line_goods_log_count_local_excess_slope"],0)
        adj=m.centred_slope(data,passenger_adjusted=True)
        self.assertTrue(adj["identified"])
        self.assertGreater(adj["within_line_goods_log_count_local_excess_slope"],0)

    def test_line_constant_traffic_fails_closed_not_false_signal(self):
        data=nine_original_physical_segments()
        for x in data:x["goods_2020"]=900
        self.assertFalse(m.centred_slope(data)["identified"])
        self.assertFalse(m.centred_slope(data,passenger_adjusted=True)["identified"])

    def test_true_segment_is_unit_not_source_site(self):
        traffic=[]
        plant={}
        for seg,source in enumerate(nine_original_physical_segments()):
            for ix in range(3):
                sid=f"plant{seg}_{ix}"
                plant[sid]={"section":source["original_line"],
                           "local_excess_beta":source["mean_local_excess"],
                           "source_same_site_jaccard_beta":.5}
                traffic.append({
                    "sid":sid,"section":source["original_line"],
                    "source_segment":f"segment{seg}",
                    "freight":source["goods_2020"],
                    "passenger":source["passenger_2020"],
                    "year":2022,
                })
        segment,qc=m.join_and_group(traffic,plant)
        self.assertEqual(qc["original_27_2020_mapped_sites"],27)
        self.assertEqual(qc["original_mapped_sites_with_valid_frequency_fixed_null"],27)
        self.assertEqual(qc["matched_original_physical_segments"],9)
        self.assertEqual(qc["matched_national_lines"],5)
        self.assertTrue(all(z["source_plant_site_count"]==3 for z in segment))
        answer=m.traffic_diagnostic(segment)
        self.assertEqual(answer["source_9_true_segments"],9)
        self.assertTrue(answer["within_line_goods_slope_with_source_exact_48_label_redistributions"]["identified"])
        self.assertEqual(answer["within_line_goods_slope_with_source_exact_48_label_redistributions"]["source_exact_2020_freight_permutations_within_original_national_lines"],48)

    def test_missing_source_site_no_zero_imputation_and_source_line_consistency(self):
        traffic=[{"sid":"A","section":"10","source_segment":"X",
                  "freight":150,"passenger":500,"year":2022}]
        # Unsupported zero-observed traffic sites must FAIL CLOSED, not be
        # manufactured as a fake 9-segment zero-exposure ecological panel.
        with self.assertRaises(ValueError):
            m.join_and_group(traffic,{})
        with self.assertRaises(ValueError):
            m.join_and_group(traffic,{"A":{
                "section":"20","local_excess_beta":.1,
                "source_same_site_jaccard_beta":.8}})


if __name__=="__main__":
    unittest.main()
