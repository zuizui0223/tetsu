"""v9.4 synthetic tests: exact per-species iNat leaveout and fixed
GBIF provider geography. None simulate measured alien seed transport,
freight carriage, roadside cutting, wind or evolving populations.
"""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_US_FAF6_firstrecord_inat_leaveout_v94 import (
    GROUPS,WINDOWS,MIN_SAME_PAIR_QUALITY_STATES,
    MIN_NON_INAT_SOURCE_ALLPLANT,MIN_NON_INAT_SOURCE_TAXA,
    split_state_source,source_quality_states,variant_temporal,
    same_pair_country_results)

def make_fixture():
    orig={"complete":True,
        "source_periods":{"before_2022":"2018,2021","after_2022":"2023,2025"},
        "source_2022_shipment_year_is_EXCLUDED_from_both_plant_windows":True,
        "states":{}}
    inat={"complete":True,"source_original_GBIF_qualified_states":48,
          "states":{}}
    for i in range(1,52):
        code=f"{i:02d}"
        pre={"1001":50,"1002":40,"1003":35,"1004":35,
             "1005":30,"1006":30,"1007":30}
        late=dict(pre,**{"1008":30,"1009":20})
        orig["states"][code]={
            "two_year_source":{"two_period_quality":i<=48},
            "source_windows":{
                "before_2022":{"n_all_Plantae_indexed_records":1000,
                               "n_US_country_GRIIS_indexed_records":sum(pre.values()),
                               "US_country_GRIIS_positive_speciesKey_record_counts":pre},
                "after_2022":{"n_all_Plantae_indexed_records":1200,
                              "n_US_country_GRIIS_indexed_records":sum(late.values()),
                              "US_country_GRIIS_positive_speciesKey_record_counts":late}
            }}
        if i<=48:
            ni=(pre if i>40 else {k:v-15 for k,v in pre.items()})
            nj=(late if i>40 else {k:v-15 for k,v in late.items()})
            inat["states"][code]={
               "complete":True,
               "source_result":{"two_GBIF_year_source_periods":{
                   "before_2022":{
                      "iNaturalist_source_country_GRIIS_taxon_keys":ni,
                      "iNaturalist_original_all_Plantae_record_source_count":990 if i>40 else 800,
                      "iNaturalist_source_country_GRIIS_record_count":sum(ni.values())},
                   "after_2022":{
                      "iNaturalist_source_country_GRIIS_taxon_keys":nj,
                      "iNaturalist_original_all_Plantae_record_source_count":1190 if i>40 else 950,
                      "iNaturalist_source_country_GRIIS_record_count":sum(nj.values())},
               }}}
    return orig,inat

class ExactGBIFProviderDecomposition(unittest.TestCase):
    def test_frozen_all_inat_other_variants_and_quality_preselection(self):
        self.assertEqual(GROUPS,("all_GBIF_original","iNaturalist_only",
                        "all_OTHER_GBIF_except_iNaturalist"))
        self.assertEqual(WINDOWS,("before_2022","after_2022"))
        self.assertEqual(MIN_SAME_PAIR_QUALITY_STATES,40)
        self.assertEqual(MIN_NON_INAT_SOURCE_ALLPLANT,100)
        self.assertEqual(MIN_NON_INAT_SOURCE_TAXA,5)
    def test_exact_per_taxon_source_record_counts_by_subtraction(self):
        original,inat=make_fixture()
        src=split_state_source(original,inat,"01","after_2022")
        all_c,n_all=src["all_GBIF_original"]
        i_c,n_inat=src["iNaturalist_only"]
        o_c,n_oth=src["all_OTHER_GBIF_except_iNaturalist"]
        self.assertEqual(n_all,n_inat+n_oth)
        for key,count in all_c.items():
            self.assertEqual(i_c.get(key,0)+o_c.get(key,0),count)
        self.assertEqual(len(o_c),9)
        self.assertEqual(n_oth,250)
    def test_original_provider_identity_source_limits_geography_to40_before_freight(self):
        orig,inat=make_fixture()
        source,eligible,missing=source_quality_states(orig,inat)
        self.assertEqual(len(source),48)
        self.assertEqual(len(eligible),40)
        self.assertEqual(len(missing),8)
        self.assertNotIn("41",eligible)
        for win in WINDOWS:
            self.assertEqual(len(source["01"][win]),3)
    def test_variant_rebuild_only_keeps_eligible_original_country_states(self):
        orig,inat=make_fixture()
        src,eligible,missing=source_quality_states(orig,inat)
        x=variant_temporal(orig,src,eligible,"all_OTHER_GBIF_except_iNaturalist")
        self.assertEqual(x["n_qualified_source_states_gate_40"],True)
        self.assertTrue(x["states"]["01"]["two_year_source"]["two_period_quality"])
        self.assertFalse(x["states"]["45"]["two_year_source"]["two_period_quality"])
        self.assertEqual(x["states"]["01"]["source_windows"]["before_2022"][
                         "n_US_country_GRIIS_indexed_records"],15*7)
        self.assertEqual(x["states"]["01"]["source_windows"]["after_2022"][
                         "n_US_country_GRIIS_observed_speciesKey_categories"],9)
    def test_provider_filtered_counts_cannot_exceed_all_source(self):
        orig,inat=make_fixture()
        inat["states"]["01"]["source_result"]["two_GBIF_year_source_periods"][
            "before_2022"]["iNaturalist_source_country_GRIIS_taxon_keys"]["1001"]=300
        with self.assertRaisesRegex(ValueError,"don't conserve"):
            split_state_source(orig,inat,"01","before_2022")
    def test_source_record_sources_never_called_40_country_replicates(self):
        out=same_pair_country_results([])
        self.assertIsNone(out["all_temporal_source_state_directions"][
            "source_responses"]["expected_20record_source_linked_first_record_fraction"]["rail_rho"])
        self.assertTrue(out["all_temporal_source_state_directions"]["NO_pairwise_independent_replicates"])
if __name__=="__main__":unittest.main()
