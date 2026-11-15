"""v16.2 synthetic arithmetic tests ONLY; no Czech/Slovenian source
plants or physical trains are simulated as observed field evidence.
"""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"tests"))
from test_Slovenia_rail_road_heldout_identity_v151 import artificial_opportunity
from test_Czech_2024_centroid_geographic_identity_v154 import fake_source_sites
from test_Czech_2023_real_road_rail_heldout_identity_v156 import synthetic_pairs
from analyze_Slovenia_rail_road_heldout_species_identity_v151 import model_records
from analyze_Czech_2023_true_rail_road_heldout_taxon_identity_v156 import pair_records
from analyze_two_country_original_four_guild_railroad_v162 import (
    guild_sets,replace_response,group_cv,source_references,
    RICHNESS_CONTROLS,MAIN_GROUPS,FROZEN_REFERENCE
)

def synthetic_slovenia():
    pairs,sites=artificial_opportunity()
    guild={}
    for s in sites:
        for t in s["known"]:
            guild[t]=("strict_annual_wind","strict_annual_other",
                      "other_known_wind","other_known_other")[t%4]
    rows,_=model_records(pairs,sites,True)
    return sites,guild,rows

def synthetic_czech():
    sites,guild=fake_source_sites()
    # Add two synthetic guild species to every site so each source-original
    # 4-group ontology is covered by field site species rather than labels only.
    for i,s in enumerate(sites):
        y="sy_"+str(i%7)
        z="sz_"+str(i%7)
        guild[y]="th_other"
        guild[z]="nonth_wind"
        s["known"]=set(s["known"])|{y,z}
    rows,_=pair_records(synthetic_pairs(sites),sites,guild,"both_Pladias_known")
    return sites,guild,rows

class SourceAuthorFunctionalSubsetTests(unittest.TestCase):
    def test_si_four_guild_codes_partition_actual_source_known_taxa_once(self):
        sites,guild,rows=synthetic_slovenia()
        g=guild_sets(sites,guild,"SI")
        self.assertEqual(set(g),set(MAIN_GROUPS))
        self.assertEqual(sum(v["source_site_guild_incidences"] for v in g.values()),
                         sum(len(s["known"]) for s in sites))
        self.assertEqual(len(rows),1770)
        self.assertTrue(all(v["observed_taxa"]>0 for v in g.values()))

    def test_cz_four_guild_codes_use_therophyte_not_slovenian_strict_annual(self):
        sites,guild,rows=synthetic_czech()
        g=guild_sets(sites,guild,"CZ")
        self.assertEqual(g["annual_wind"]["source_label"],"th_wind")
        self.assertEqual(g["nonannual_other_label"] if "nonannual_other_label" in g else
                         g["nonannual_nonwind"]["source_label"],"nonth_other")
        self.assertEqual(sum(v["source_site_guild_incidences"] for v in g.values()),
                         sum(len(s["known"]) for s in sites))
        self.assertTrue(all(v["n_original_site_with_guild_taxa"]==35 for v in g.values()))

    def test_not_identifiable_if_guild_is_absent_in_all_but_one_station(self):
        sites,guild,base=synthetic_czech()
        group=[set() for _ in sites]
        group[0]={"fake_wind_diaspore"}
        rows,support=replace_response("CZ",base,sites,group)
        self.assertEqual(support["source_guild_pair_comparisons_with_at_least_one_original_species"],
                         len(sites)-1)
        self.assertFalse(group_cv("CZ",rows)["identified"])

    def test_do_not_make_false_jaccard_zero_for_two_empty_guild_sites(self):
        sites,guild,base=synthetic_czech()
        group=[{"fake"} if i<32 else set() for i in range(35)]
        rows,qc=replace_response("CZ",base,sites,group)
        self.assertEqual(qc["source_pairs_excluded_because_both_sites_have_zero_original_guild_species"],3)
        self.assertTrue(all(x["jaccard"] in (0.,1.) for x in rows))
        self.assertTrue(group_cv("CZ",rows)["identified"])

    def test_original_source_SI_and_CZ_pair_support_keep_geo_and_route_exact(self):
        for study,fn in (("SI",synthetic_slovenia),("CZ",synthetic_czech)):
            sites,guild,base=fn()
            group=guild_sets(sites,guild,study)["annual_wind"]["site_sets"]
            sub={}
            for mode in RICHNESS_CONTROLS:
                rows,qc=replace_response(study,base,sites,group,mode)
                sub[mode]=rows
                self.assertTrue(qc["unobserved_original_species_NOT_inserted_as_absences"])
                self.assertTrue(all("jaccard" in x for x in rows))
                if study=="CZ":
                    self.assertTrue(group_cv(study,rows)["identified"])
                else:
                    self.assertTrue(group_cv(study,rows)["identified"])
            self.assertEqual([(x["i"],x["j"]) for x in sub[RICHNESS_CONTROLS[0]]],
                             [(x["i"],x["j"]) for x in sub[RICHNESS_CONTROLS[1]]])
            for a,b in zip(*sub.values()):
                self.assertEqual(a["jaccard"],b["jaccard"])
                self.assertEqual(a["features"]["geo_ln"],b["features"]["geo_ln"])
                self.assertTrue("features" in b)

    def test_strict_real_author_source_cohort_gate_does_not_accept_synthetic(self):
        sites,guild,base=synthetic_slovenia()
        with self.assertRaisesRegex(ValueError,"cohort drifted"):
            source_references("SI",{},base,sites)
        self.assertEqual(FROZEN_REFERENCE["CZ"]["n_supported_site_pairs"],594)
        self.assertEqual(FROZEN_REFERENCE["SI"]["n_supported_site_pairs"],3896)

    def test_missing_or_changed_original_trait_ontology_fails_closed(self):
        sites,guild,base=synthetic_czech()
        bad=dict(guild)
        bad[next(iter(bad))]="unverified_hypothetical_wind"
        with self.assertRaisesRegex(ValueError,"ontology changed"):
            guild_sets(sites,bad,"CZ")
        with self.assertRaisesRegex(ValueError,"unknown"):
            replace_response("XX",base,sites,[set() for _ in sites])

if __name__=="__main__":unittest.main()
