"""Synthetic US state 2022 FAF crop OD / 2018–25 recorded country GRIIS
flora source-only regression tests. NEVER treat ecological associations
as freight seed transport, selection, wind or heritable adaptation.
"""
import sys,io,unittest
from pathlib import Path
from types import SimpleNamespace
from html import escape
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))

from analyze_FAF6_GRIIS_US_state_OD_alignment_v83 import (
    MIN_STATES,K,COHORTS,METRICS,parse_census_centroids,geodesic_km,
    validated_botanical_states,original_crop_OD,state_pair_rows,
    partial_rank_rail)

def flora_fixture():
    source={"complete":True,"date_period":"2018_2025",
            "n_requested_state_fips_codes":51,"states":{}}
    for i in range(1,52):
        fips=f"{i:02d}"
        counts={str(1000+j):15 for j in range(5)}
        source["states"][fips]={
            "official_Census_state_name":f"State {i}",
            "GBIF_2018_2025_GRIIS_country_taxa_source":{
                "source_eligible":i<=49,
                "reason_not_qualified":[] if i<=49 else ["source-status"],
                "source_GRIIS_speciesKey_indexed_counts":counts,
                "source_all_Plantae_indexed_records":3000+100*i,
            }}
    return source

def tigertable(source):
    rows=[]
    for code,v in source["states"].items():
        i=int(code)
        rows.append({"STATE":code,"NAME":v["official_Census_state_name"],
         "STUSAB":f"S{i}","CENTLAT":30+i*.25,
         "CENTLON":-120+i*.5,"AREALAND":100000000+i*123456,
         "REGION":(i%4)+1,"DIVISION":(i%9)+1})
    columns=list(rows[0])
    header="<tr>"+"".join("<th>"+escape(key)+"</th>" for key in columns)+"</tr>"
    body="".join("<tr>"+"".join("<td>"+escape(str(row[key]))+"</td>" for key in columns)+"</tr>" for row in rows)
    return "<html><body><table>"+header+body+"</table></body></html>"

def freight_fixture():
    return {
        "complete":True,"source_hash_verified_vs_prior_v79":True,
        "source_trace":{"previous_verified_v79_total_and_modes_reconciled":True},
        "source_year":2022,
        "raw_official_state_OD_2022_crop_estimated_tonnage":{
            "RAIL":[{"o":"01","d":"02","estimated_thousand_short_tons":500},
                    {"o":"02","d":"01","estimated_thousand_short_tons":50}],
            "TRUCK":[{"o":"01","d":"02","estimated_thousand_short_tons":200},
                     {"o":"03","d":"02","estimated_thousand_short_tons":300}]
        }}

class FAF6GriisStateCrossSource(unittest.TestCase):
    def test_predeclared_responses_and_noncausal_51_state_source(self):
        self.assertEqual(K,20)
        self.assertEqual(MIN_STATES,40)
        self.assertEqual(len(COHORTS),3)
        self.assertEqual(METRICS,("rarefied20_expected_inclusion_overlap",
                                  "raw_recorded_GRIIS_Jaccard"))
    def test_real_Census_TIGER_fips_2025_state_geography_join_contract(self):
        src=flora_fixture()
        geo=parse_census_centroids(tigertable(src),src["states"])
        self.assertEqual(len(geo),51)
        self.assertTrue(geo["01"]["land_sq_m"]>0)
        self.assertGreater(geodesic_km(geo["01"],geo["02"]),0)
    def test_status_ineligible_state_not_imputed_empty_plant_list(self):
        src=flora_fixture()
        yes,no=validated_botanical_states(src)
        self.assertEqual(len(yes),49)
        self.assertEqual(set(no),{"50","51"})
        self.assertEqual(sum(yes["01"]["counts"].values()),75)
    def test_directional_economic_freight_is_symmetrized_only_for_symmetric_flora(self):
        flora=flora_fixture()
        src,excluded=validated_botanical_states(flora)
        geo=parse_census_centroids(tigertable(flora),flora["states"])
        rows=state_pair_rows(src,geo,freight_fixture())
        p=next(q for q in rows if q["a"]=="01" and q["b"]=="02")
        self.assertEqual(p["rail_estimated_thousand_short_tons"],550)
        self.assertEqual(p["truck_estimated_thousand_short_tons"],200)
        self.assertTrue(0<=p["responses"]["rarefied20_expected_inclusion_overlap"]<=1)
        self.assertEqual(len(rows),49*48//2)
    def test_source_flora_with_no_variation_does_not_invent_railway_effect(self):
        flora=flora_fixture()
        src,excluded=validated_botanical_states(flora)
        geo=parse_census_centroids(tigertable(flora),flora["states"])
        rows=state_pair_rows(src,geo,freight_fixture())
        self.assertIsNone(partial_rank_rail(rows,METRICS[0]))
    def test_FA6_tonnes_must_not_be_doubled_or_imputed_as_actual_rail_tracks(self):
        freight=freight_fixture()
        modes=original_crop_OD(freight)
        self.assertEqual(modes["RAIL"][("01","02")],500)
        freight["source_year"]=2017
        with self.assertRaisesRegex(ValueError,"source base year"):
            original_crop_OD(freight)
if __name__=="__main__":unittest.main()
