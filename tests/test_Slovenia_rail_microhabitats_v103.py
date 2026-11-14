"""v10.3 synthetic matched SAME railway location and side ballast
microhabitat zone tests. Do not claim actual seed/carriage, mowing/gust,
Native status from blank codes, or equal physical zone plot area.
"""
import sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Slovenia_rail_microhabitats_v103 import (
    SOURCE_SHA256, SOURCE_BYTES,SHEETS,SHEET_ROWS,AREAS,
    trait_shares,float_field,pair_measure,locations_independent,
    location_summary,section_summary,validate_original_keys)

def fake_traits():
    return {
      "AN":{"Life.form":"Th","Lifespan.categories":"annual",
            "Nativeness.status":"invasive alien species"},
      "MX":{"Life.form":"Th/He","Lifespan.categories":"annual/biennial",
            "Nativeness.status":None},
      "PE":{"Life.form":"He","Lifespan.categories":"perennial",
            "Nativeness.status":None},
    }
def source_pair(lid,side,coverBS,coverBT,strictAnnual=.5,strictThero=.4):
    return {
     "location_id":lid,"track_side":side,"section_number":10,
     "survey_year":2022,"track_type":"SINGLE",
     "source_cover_value_above_100_on_either_side":max(coverBS,coverBT)>100,
     "BS_cover_source_numeric_not_guaranteed_percentage":coverBS,
     "BT_cover_source_numeric_not_guaranteed_percentage":coverBT,
     "BS":{
        "strict_therophyte_share_of_known_lifeforms":strictThero,
        "strict_annual_share_of_known_lifespans":strictAnnual,
        "species_richness_per_original_zone_not_area_standardized":5,
        "broad_therophyte_share_including_mixed_Th_lifeforms":strictThero+.1,
        "broad_annual_share_including_mixed_lifespans":strictAnnual+.1,
        "explicitly_labelled_alien_share_of_all_recorded_taxa":.1},
     "BT":{
        "strict_therophyte_share_of_known_lifeforms":.2,
        "strict_annual_share_of_known_lifespans":.3,
        "species_richness_per_original_zone_not_area_standardized":10,
        "broad_therophyte_share_including_mixed_Th_lifeforms":.3,
        "broad_annual_share_including_mixed_lifespans":.4,
        "explicitly_labelled_alien_share_of_all_recorded_taxa":.2},
    }


class SlovenianActualRailwaySourceAlgebra(unittest.TestCase):
    def test_frozen_official_field_survey_source_not_MIREN_road(self):
        self.assertEqual(len(SOURCE_SHA256),64)
        self.assertEqual(SOURCE_BYTES,401279)
        self.assertEqual(SHEETS,("README","locations","zones","species","traits","species_codebook"))
        self.assertEqual(SHEET_ROWS,{"locations":92,"zones":482,"species":6331,"traits":312,"species_codebook":312})
        self.assertEqual(AREAS,("BS","BT"))
    def test_strict_annual_and_therophyte_are_not_mixed_life_categories(self):
        x=trait_shares({"AN","MX","PE"},fake_traits())
        self.assertAlmostEqual(x["strict_annual_share_of_known_lifespans"],1/3)
        self.assertAlmostEqual(x["broad_annual_share_including_mixed_lifespans"],2/3)
        self.assertAlmostEqual(x["strict_therophyte_share_of_known_lifeforms"],1/3)
        self.assertAlmostEqual(x["broad_therophyte_share_including_mixed_Th_lifeforms"],2/3)
    def test_blank_nativeness_is_unknown_not_native(self):
        s=trait_shares({"AN","MX","PE"},fake_traits())
        self.assertAlmostEqual(s["explicitly_labelled_alien_share_of_all_recorded_taxa"],1/3)
        self.assertAlmostEqual(s["fraction_unclassified_Nativeness_status_of_all_recorded_taxa"],2/3)
        self.assertEqual(trait_shares(set(),fake_traits()),{"no_taxa":True})
    def test_exact_same_side_BS_minus_BT_cover_and_trait_delta(self):
        r=source_pair("s1","L",12,72,.5,.4)
        self.assertAlmostEqual(pair_measure(r,"cover"),-60)
        self.assertAlmostEqual(pair_measure(r,"richness"),-5)
        self.assertAlmostEqual(pair_measure(r,"annual_strict"),.2)
        self.assertAlmostEqual(pair_measure(r,"therophyte_strict"),.2)
        self.assertAlmostEqual(pair_measure(r,"explicit_alien_label"),-.1)
    def test_each_railway_location_averages_left_right_before_macro_summary(self):
        xs=[]
        for i in range(1,90):
            xs.append(source_pair(f"site{i}","L",10,70))
            if i<84:
                xs.append(source_pair(f"site{i}","R",20,80))
        self.assertEqual(len(xs),172)
        locations=locations_independent(xs)
        self.assertEqual(len(locations),89)
        self.assertEqual(locations[0]["n_matched_same_side_BS_BT_pairs"],2)
        self.assertEqual(location_summary(locations,"cover")["n_locations_BT_greater"],89)
        self.assertEqual(location_summary(locations,"annual_strict")[
            "n_locations_BS_greater"],89)
    def test_three_locations_per_segment_required_and_no_pseudoindependent_sides(self):
        xs=[source_pair(f"site{i}","L",10,70) for i in range(1,90)]
        locations=locations_independent(xs)
        report=section_summary(locations,"cover",3)
        self.assertEqual(report["n_eligible_original_rail_segments_NOT_countries"],1)
        self.assertEqual(report["n_segments_BT_greater"],1)
        self.assertLess(report["median_rail_segment_equal_mean_BS_minus_BT"],0)
    def test_zone_cover_above100_kept_source_units_not_claimed_percentage(self):
        self.assertEqual(float_field(200,"cover",allow200=True),200.)
        with self.assertRaisesRegex(ValueError,"numeric source"):
            float_field(-1,"cover")
        r=source_pair("s1","L",35,200)
        self.assertTrue(r["source_cover_value_above_100_on_either_side"])
        self.assertEqual(pair_measure(r,"cover"),-165)
    def test_only_true_original_source_joins_and_taxonomic_richness(self):
        locs=[{"location_id":f"S{i}","track_section_number":i%15,
               "date":44780,"track_type":"SINGLE"} for i in range(92)]
        zones=[]
        observed=[]
        for i in range(92):
            for side in ("L","R"):
                if len(zones)>=482:break
                zones.append({"location_id":f"S{i}","zone_id":side+"BT",
                              "area":"BT","area_side":side,"richness":1,"cover":50})
                observed.append({"location_id":f"S{i}","zone_id":side+"BT","EPPO_code":"X"})
        # Need EXACT source 482 zone rows. Duplicate taxon in one zone
        # can be valid and must not increase reported richness.
        while len(zones)<482:
            i=(len(zones)-184)%92
            idx=len(zones)
            key=f"ZONE{idx}"
            zones.append({"location_id":f"S{i}","zone_id":key,
                          "area":"SZ","richness":1,"cover":1})
            observed.append({"location_id":f"S{i}","zone_id":key,"EPPO_code":"X"})
        for k in range(6331-len(observed)):
            observed.append(dict(observed[k%len(observed)]))
        traits=[{"EPPO_code":"X","Life.form":"Th","Lifespan.categories":"annual",
                  "Nativeness.status":"invasive alien species"}]
        traits += [{"EPPO_code":f"T{i}","Life.form":None,
                    "Lifespan.categories":None,"Nativeness.status":None}
                   for i in range(311)]
        code=[{"EPPO_code":z["EPPO_code"]} for z in traits]
        _,_,_,_,_,qa=validate_original_keys({
            "locations":locs,"zones":zones,"species":observed,
            "traits":traits,"species_codebook":code})
        self.assertEqual(qa["n_original_species_zone_occurrence_rows"],6331)
        self.assertGreater(qa["n_original_duplicate_species_EPPO_rows_within_field_zone_deduplicated"],0)
        self.assertTrue(qa["all_482_original_zone_richness_equals_unique_species_codes"])

if __name__=="__main__":unittest.main()
