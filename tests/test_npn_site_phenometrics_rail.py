"""Synthetic site-phenometric screening checks; do not imply train biology."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_npn_site_phenometrics_rail import summarize_state_year


def record(site,species,phase="Open flowers",year=2024):
    return {
        "site_id":site,"species_id":species,
        "phenophase_description":phase,"year":year,
        "first_yes_date":f"{year}-05-10",
    }


class SitePhenoTests(unittest.TestCase):
    def test_shared_taxon_requires_real_near_and_far_records(self):
        bins={1:"near",2:"far",3:"middle_excluded",4:"near"}
        raw=[record(1,82),record(2,82),record(1,82),
             record(3,97),record(4,97)]
        out=summarize_state_year(raw,bins,"MD",2024)
        self.assertEqual(out["n_species_with_positive_phenometric_both_near_far"],1)
        self.assertEqual(out["species_eligible_for_later_individual_status_lookups"][0]["npn_species_id"],82)
        self.assertEqual(out["n_registered_sites_with_positive_Open_flowers_near"],2)
        self.assertEqual(out["n_registered_sites_with_positive_Open_flowers_far"],1)
        self.assertEqual(out["quality"]["duplicate_site_species_rows"],1)
        self.assertEqual(out["quality"]["intermediate_rail_distance"],1)

    def test_registered_but_unobserved_does_not_count(self):
        bins={1:"near",2:"far"}
        out=summarize_state_year([record(2,82)],bins,"MD",2024)
        self.assertEqual(out["n_registered_sites_with_positive_Open_flowers_near"],0)
        self.assertEqual(out["n_registered_sites_with_positive_Open_flowers_far"],1)
        self.assertEqual(out["n_species_with_positive_phenometric_both_near_far"],0)

    def test_missing_or_unknown_ids_excluded(self):
        out=summarize_state_year([
            record(1,None),record(90,82),record(1,82)
        ],{1:"near"},"MD",2024)
        self.assertEqual(out["quality"]["invalid_site_or_species_id"],1)
        self.assertEqual(out["quality"]["site_not_in_registered_state_frame"],1)


if __name__=="__main__":
    unittest.main()
