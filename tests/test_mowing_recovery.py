"""Synthetic fixed plot photos verify mowing alignment without causal claims."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from tetsu.mowing_recovery import summarize_mowing_roi


def photo(i,date,site="railA",plot="P1",taxon="Fictitia exemplaris",
          cover=0.5,flower=0.2,valid=True):
    return {
        "image_id":str(i),"photo_date":date,"site_id":site,
        "plot_id":plot,"taxon":taxon,
        "quality_status":"patch_pass" if valid else "restricted",
        "valid_for_fixed_plot_flower_fraction":str(valid),
        "plant_frame_fraction":cover,
        "flower_area_over_plant_area":flower,
    }


class MowingROI(unittest.TestCase):
    def setUp(self):
        self.events=[{
            "site_id":"railA","event_date":"2026-06-10",
            "event_type":"cutting","evidence":"documented_notice"
        }]
        self.photos=[
            photo(1,"2026-05-20",cover=0.65,flower=0.12),
            photo(2,"2026-06-08",cover=0.60,flower=0.16),
            photo(3,"2026-06-24",cover=0.45,flower=0.01),
            photo(4,"2026-07-15",cover=0.55,flower=0.08)
        ]

    def test_pre_after_descriptive_changes(self):
        result,rejected=summarize_mowing_roi(self.photos,self.events)
        self.assertEqual(len(result),1)
        r=result[0]
        self.assertEqual(r["eligibility"],"pass")
        self.assertEqual(r["n_pre"],2)
        self.assertEqual(r["n_recovery"],2)
        self.assertAlmostEqual(r["delta_vegetation_frame_occupancy"],-0.125)
        self.assertAlmostEqual(r["delta_floral_area_per_vegetation"],-0.095)
        self.assertEqual(rejected,{})

    def test_event_overlap_restricts(self):
        events=self.events+[{
            "site_id":"railA","event_date":"2026-07-01",
            "event_type":"cutting","evidence":"documented_log"
        }]
        result,_=summarize_mowing_roi(self.photos,events)
        self.assertTrue(all(x["eligibility"]=="restricted" for x in result))
        self.assertTrue(all(x["restriction_reason"]=="overlapping_cutting_events" for x in result))

    def test_not_assuming_no_photo_is_zero(self):
        result,_=summarize_mowing_roi(self.photos[:3],self.events)
        self.assertEqual(result[0]["eligibility"],"restricted")
        self.assertEqual(result[0]["restriction_reason"],"insufficient_repeated_fixed_plot_photos")

    def test_unreviewed_patch_excluded(self):
        result,rejected=summarize_mowing_roi(
            self.photos+[photo(5,"2026-06-25",valid=False)],self.events
        )
        self.assertEqual(result[0]["n_recovery"],2)
        self.assertEqual(rejected["not_verified_fixed_conspecific_patch"],1)

    def test_missing_evidence_refused(self):
        with self.assertRaises(ValueError):
            summarize_mowing_roi(self.photos,[{
                "site_id":"railA","event_date":"2026-06-10",
                "event_type":"cutting","evidence":""
            }])

    def test_duplicate_photos_refused(self):
        with self.assertRaises(ValueError):
            summarize_mowing_roi(self.photos+self.photos[:1],self.events)


if __name__=="__main__":
    unittest.main()
