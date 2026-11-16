"""v14.3 independent original Czech 2024 station macro source math tests.

All synthetic station fixtures in this file are correctness checks only.
Actual biological results only from hash-pinned publisher 39-station S1/S2.
"""
from __future__ import annotations
import sys,unittest,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Czech_2024_39station_functional_identity_v143 import (
    GROUPS,czech_guild,jaccard,function_sim,site_pairs,swap_null,
    read_original_czech_publisher,current_sites
)


class GenuineCzechStationGuildControls(unittest.TestCase):
    def test_exact_published_plant_taxonomy_only_not_genus_guessing(self):
        self.assertEqual(czech_guild("Ailanthus altissima",
                 {"Ailanthus altissima":"Epilobium"},
                 {"Ailanthus altissima":"phanerophyte"}),"nonth_wind")
        self.assertEqual(czech_guild("Brassica napus",
                 {"Brassica napus":"Allium"},
                 {"Brassica napus":"therophyte"}),"th_other")
        self.assertIsNone(czech_guild("unknown",
                 {"Brassica napus":"Allium"},
                 {"Brassica napus":"therophyte"}))

    def test_functional_similarity_can_be_one_when_taxonomic_identity_zero(self):
        guild={"a":"th_wind","b":"th_wind","c":"nonth_other","d":"nonth_other"}
        a={"a","c"};b={"b","d"}
        self.assertEqual(jaccard(a,b),0.)
        self.assertAlmostEqual(function_sim(a,b,guild),1.)

    def test_original_stretch_as_5_not_39_independent_country_replicates(self):
        species={"a":"th_wind","b":"th_wind",
                 "c":"nonth_other","d":"nonth_other"}
        sites=[]
        for i in range(39):
            taxa=(("a","c") if i%3==0 else
                  ("b","d") if i%3==1 else ("a","d"))
            sites.append({
              "site":f"S{i}","section":f"block{i%5}",
              "known":set(taxa),"all_original_taxa":set(taxa),
              "area_ha":float(i+1)})
        obs=site_pairs(sites,species)
        self.assertEqual(sum(obs["taxonomic_pladias_known"][k]
                             for k in ("n_pairs_same","n_pairs_other")),741)
        self.assertTrue(obs["functional_four_coarse_guild_overlap"]["gap"] is not None)
        d=swap_null(sites,species,draws=24,burnin=100,thin=90,seed=4)
        self.assertTrue(d["identified"])
        self.assertEqual(d["draws_autocorrelated"],24)
        self.assertEqual(d["swap_proposals"],100+24*90)
        self.assertTrue(d["source_original_station_size_and_species_frequency_margins_preserved"])

    def test_source_no_trait_unknown_and_name_exclusion(self):
        flora={"aa":{"Species one","Species two"},
               "bb":{"Species one","Species three"}}
        meta={"aa":{"stretch":"one","original_source_station_total_area_ha":5},
              "bb":{"stretch":"two","original_source_station_total_area_ha":6}}
        d={"Species one":"Epilobium","Species three":"Phragmites"}
        l={"Species one":"therophyte","Species three":"hemicryptophyte"}
        sites,known,qc=current_sites(flora,meta,d,l)
        self.assertEqual(known,{"Species one":"th_wind","Species three":"nonth_wind"})
        self.assertEqual(qc["n_original_incidence_with_both_exact_Pladias_traits"],3)
        self.assertEqual(qc["n_original_contemporary_site_taxon_incidence_all_species"],4)
        self.assertEqual(qc["actual_original_species_and_site_trait_coverage"],.75)

    def test_fake_original_stations_refused_by_publisher_sha(self):
        with tempfile.TemporaryDirectory() as t:
            a=Path(t)/"s1.xlsx";b=Path(t)/"s2.xlsx"
            a.write_bytes(b"not the original author paper")
            b.write_bytes(b"nor the corrected S2")
            with self.assertRaisesRegex(ValueError,"SHA/bytes"):
                read_original_czech_publisher(a,b)


if __name__=="__main__":unittest.main()
