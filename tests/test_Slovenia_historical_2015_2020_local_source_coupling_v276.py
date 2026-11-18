"""v27.6 synthetic time-resolved actual train predictors with fixed plant response QA."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Slovenia_2015_2020_real_train_local_coupling_temporal_v276 as m


def source_years():
    line_count=[("10",2),("20",1),("34",1),("50",4),("80",1)]
    a=[];b=[]
    idx=0
    for line,n in line_count:
        for j in range(n):
            z={"original_line":line,"source_plant_site_count":j+1,
               "mean_local_excess":.06+.011*idx,
               "mean_observed_Jaccard_beta":.5+.01*idx,
               "goods_2020":150+idx*100,
               "passenger_2020":250+idx*15}
            y={**z,"goods_2020":200+idx*110,"passenger_2020":300+idx*20}
            a.append(z);b.append(y);idx+=1
    return a,b


class GenuineHistoricalIntervalTests(unittest.TestCase):
    def test_identical_plant_source_same_nine_segments_two_predictor_years(self):
        x,y=source_years()
        s=m.paired_historical(x,y)
        self.assertEqual(s["original_source_physical_segments"],9)
        self.assertTrue(s["original_identical_physical_segment_plant_response_in_two_years"])
        self.assertIsNotNone(s["2015_2020_same_physical_segments_goods_spearman"])
        self.assertIsNotNone(s["2015_2020_same_physical_segments_passenger_spearman"])

    def test_changed_plants_fail_year_sensitivity_not_misreported_train_trend(self):
        x,y=source_years()
        y[0]["mean_local_excess"]+=.2
        with self.assertRaises(ValueError):
            m.paired_historical(x,y)
        x,y=source_years()
        y=y[:-1]
        with self.assertRaises(ValueError):
            m.paired_historical(x,y)

    def test_ambiguous_source_line_segment_identity_fails_closed(self):
        x,y=source_years()
        y[1]["mean_local_excess"]=y[0]["mean_local_excess"]
        y[1]["mean_observed_Jaccard_beta"]=y[0]["mean_observed_Jaccard_beta"]
        y[1]["source_plant_site_count"]=y[0]["source_plant_site_count"]
        with self.assertRaises(ValueError):
            m.paired_historical(x,y)

    def test_missing_trains_prohibit_pretending_old_year_equals_new(self):
        # No 2015 source segment match in authentic key-space is not imputed
        with self.assertRaises(ValueError):
            m.historical_site_rematch(
                [{"source_segment":"a+b"} for _ in range(27)],
                [],
                []
            )


if __name__=="__main__":
    unittest.main()
