"""v28.3 synthetic source-incidence null QA; no original author species/locations."""
import random
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_Czech_35station_distance_matched_species_prevalence_null_v283 as m
import audit_Czech_2024_true_35_stations_stretch_geodistance_matched_v282 as matching


def original_shaped_synthetic_stations():
    sites=[]
    for section in range(5):
        for i in range(7):
            sites.append({
                "source_site":f"{section}-{i}",
                "source_stretch":str(section),
                "lat":49.+section*.06,"lon":14.+i*.07,
                "area_ha":float(i+3),
                "taxa":{"common",f"stretch{section}",f"plant_{i}",
                        f"cross{(i+section)%3}"}
            })
    return sites


class OriginalCzechTaxonSourceFrequencyNullQA(unittest.TestCase):
    def test_exact_source_all_and_original_stretch_plant_frequency_margins(self):
        sites=original_shaped_synthetic_stations()
        original=[set(z["taxa"]) for z in sites]
        for scope in m.SCOPE:
            groups=m.source_groups(sites,scope)
            before=m.exact_counts(original,groups)
            rng=random.Random(35)
            from audit_RAIL_ROAD_source_frequency_fixed_reference_null_v270 import swaps
            changed=[set(x) for x in original]
            for group in groups:
                temp=[changed[i] for i in group]
                swaps(temp,rng,200)
            self.assertEqual(m.exact_counts(changed,groups),before)
        self.assertEqual(len(m.source_groups(sites,m.SCOPE[1])),5)
        with self.assertRaises(ValueError):
            m.source_groups(sites,"try_posthoc_habitat_with_small_positive_data")

    def test_fixed_source_geographic_matches_and_true_source_stretch_gap(self):
        sites=original_shaped_synthetic_stations()
        match,qa=matching.candidates(sites,160,10,None)
        self.assertEqual(len(match),35)
        original=[set(x["taxa"]) for x in sites]
        obs,bygroup=m.fixed_source_matched_gap(original,sites,match)
        self.assertEqual(len(bygroup),5)
        self.assertAlmostEqual(obs,matching.aggregate(match,qa)[
            "real_stretch_equal_same_minus_other_jaccard"])

    def test_mcmc_no_fake_independent_station_pairs_or_evolutionary_effect(self):
        sites=original_shaped_synthetic_stations()
        source,_=matching.candidates(sites,160,10,None)
        for scope in m.SCOPE:
            x=m.one_source_null(sites,source,scope,seed=988,
                                draws=50,burnin=75,thin=12)
            self.assertEqual(x["n_source_correlated_MCMC_draws"],100)
            self.assertTrue(x["original_source_richness_and_species_frequency_exactly_fixed"])
            self.assertEqual(len(x["two_independent_chain_acceptance_and_mean"]),2)
            self.assertTrue(-1<=x["original_observed_minus_frequency_null"]<=1)
            self.assertEqual(x["source_real_original_station_cohort"],35)

    def test_minimum_draw_gate(self):
        sites=original_shaped_synthetic_stations()
        match,_=matching.candidates(sites,160,10,None)
        with self.assertRaises(ValueError):
            m.one_source_null(sites,match,m.SCOPE[0],seed=5,draws=10)


if __name__=="__main__":
    unittest.main()
