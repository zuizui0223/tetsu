"""v28.2 independent Czech author botanical STATION and 2026 registry centroid
geography synthetic math checks; NO field data or fitted outcomes in tests."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Czech_2024_true_35_stations_stretch_geodistance_matched_v282 as m


def fake_czech_stations():
    stations=[]
    for section in range(5):
        for j in range(7):
            stations.append({
                "source_site":f"{section}-{j}","source_stretch":str(section),
                "lat":49+section*.06,"lon":14+j*.07,
                "area_ha":1.0+j*.1,
                "taxa":{1,10+section,100+j}
            })
    return stations


class CzechSourceStationGeodistanceTests(unittest.TestCase):
    def test_different_source_stretches_and_same_source_not_self(self):
        stations=fake_czech_stations()
        pairs,qa=m.candidates(stations,160,10,None)
        self.assertEqual(qa["source_n_actual_geolocated_author_stations"],35)
        self.assertEqual(qa["source_n_original_stretch_labels"],5)
        self.assertGreater(len(pairs),20)
        for p in pairs:
            i,j,k=p["source_station_i"],p["source_station_same"],p["source_station_other"]
            self.assertTrue(len({i,j,k})==3)
            self.assertEqual(stations[i]["source_stretch"],stations[j]["source_stretch"])
            self.assertNotEqual(stations[i]["source_stretch"],stations[k]["source_stretch"])
            self.assertLessEqual(p["geographic_matching_gap_km"],10)
            self.assertAlmostEqual(
                p["source_all_taxa_same_minus_other"],
                p["original_source_jaccard_same"]-p["original_source_jaccard_other"])

    def test_paired_source_match_not_fitted_using_plant_species(self):
        s=fake_czech_stations()
        a,_=m.candidates(s,160,10,None)
        # Modify species identities while preserving station geography.
        for r in s:r["taxa"]={r["source_site"],"fake"}
        b,_=m.candidates(s,160,10,None)
        refs=lambda pairs:[(p["source_station_i"],p["source_station_same"],p["source_station_other"]) for p in pairs]
        self.assertEqual(refs(a),refs(b))

    def test_source_area_sensitivity_has_no_additional_matches(self):
        s=fake_czech_stations()
        baseline,_=m.candidates(s,160,10,None)
        strict,_=m.candidates(s,160,10,.5)
        self.assertLessEqual(len(strict),len(baseline))
        for p in strict:
            j,k=p["source_station_same"],p["source_station_other"]
            import math
            self.assertLessEqual(abs(math.log(s[j]["area_ha"])-math.log(s[k]["area_ha"])),.5)

    def test_no_source_overlap_is_unidentified_not_negative(self):
        s=fake_czech_stations()
        pairs,qa=m.candidates(s,.00001,.00001,None)
        out=m.aggregate(pairs,qa)
        self.assertFalse(out["identified"])
        self.assertNotIn("real_stretch_equal_same_minus_other_jaccard",out)

    def test_true_stretch_not_station_pair_n_inference(self):
        s=fake_czech_stations()
        pairs,qa=m.candidates(s,160,10,None)
        r=m.aggregate(pairs,qa)
        self.assertTrue(r["identified"])
        self.assertEqual(r["source_original_stretches_with_real_geographic_match"],5)
        self.assertEqual(r["source_matched_original_station_anchors"],35)
        self.assertEqual(len(r["original_stretch_leave_one_out_minmax"]),2)


if __name__=="__main__":
    unittest.main()
