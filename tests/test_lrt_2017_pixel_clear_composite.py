"""Synthetic tests prevent mixing S2 baseline dates across near/far sites."""
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from measure_lrt_2017_pixel_clear_composite import baseline_pair_first_clear


def value(v):
    return {"provider_flag_branch":v}


class StrictSameDateBaselineTests(unittest.TestCase):
    def test_near_far_must_come_from_same_source_date(self):
        pairs=[("n1","f1"),("n2","f2"),("n3","f3")]
        by_date={
            "2017-07-22":{
                "n1":value(.3),"f1":value(.4),
                "n2":value(.2),  # far f2 is clouded on 7/22
                "f3":value(.5)
            },
            "2017-07-07":{
                "n2":value(.25),"f2":value(.35),
                "n3":value(.1),  # near n3 clear 7/07, far f3 was 7/22
            }
        }
        chosen,counts,pair_dates=baseline_pair_first_clear(pairs,by_date)
        self.assertEqual(counts,{"2017-07-22":1,"2017-07-07":1})
        self.assertEqual(len(chosen),4)
        self.assertNotIn("n3",chosen)
        self.assertNotIn("f3",chosen)
        self.assertEqual(pair_dates[("n1","f1")],"2017-07-22")
        self.assertEqual(pair_dates[("n2","f2")],"2017-07-07")

    def test_selection_is_not_by_positive_ndvi_effect(self):
        pairs=[("n","f")]
        by_date={
            "2017-07-22":{"n":value(.05),"f":value(.9)},
            "2017-07-07":{"n":value(.9),"f":value(.05)}
        }
        chosen,counts,_=baseline_pair_first_clear(pairs,by_date)
        self.assertEqual(counts,{"2017-07-22":1})
        self.assertAlmostEqual(chosen["n"]["harmonized"],.05)

    def test_no_valid_pair_does_not_invent_absence(self):
        x,counts,_=baseline_pair_first_clear(
            [("n","f")],
            {"2017-07-22":{"n":value(.2)}}
        )
        self.assertEqual(x,{})
        self.assertEqual(counts,{})


if __name__=="__main__":
    unittest.main()
