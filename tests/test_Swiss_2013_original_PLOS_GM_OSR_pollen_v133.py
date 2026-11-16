"""v13.3 original PLOS 2014 Basel 2013 GM plants vs seed pollen direct experimental genetics."""
from __future__ import annotations
import sys,unittest,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Swiss_2013_original_PLOS_GM_OSR_pollen_v133 import (
    read_original_leaf_source_table,verify_actual_PLOS_fulltext_source,
    STATION,PORT,run,
)

CSV=Path(__file__).resolve().parents[1]/"data"/"Swiss_Schulze2014_PLOS_2013_rail_port_GM_OSR_leaf_samples_v133.csv"


class Swiss2013GeneticSourceCohort(unittest.TestCase):
    def test_source_original_GT73_outcomes_and_other_separate_pollen_outcrossing(self):
        r=read_original_leaf_source_table(CSV)
        self.assertEqual(len(r),2)
        by={x["original_site"]:x for x in r}
        self.assertEqual(by[STATION]["source_original_OSR_leaf_samples"],279)
        self.assertEqual(by[STATION]["source_confirmed_GT73_OSR_leaf_samples"],64)
        self.assertIsNone(by[STATION][
            "source_glufosinate_bar_positive_OSR_leaf_samples_if_reported"])
        self.assertEqual(by[PORT]["source_original_OSR_leaf_samples"],2291)
        self.assertEqual(by[PORT]["source_confirmed_GT73_OSR_leaf_samples"],106)
        self.assertEqual(by[PORT][
            "source_glufosinate_bar_positive_OSR_leaf_samples_if_reported"],22)
        self.assertAlmostEqual(by[STATION][
            "source_fraction_GT73_positive_among_original_OSR_leaf_samples"],64/279)
        self.assertAlmostEqual(by[PORT][
            "source_fraction_GT73_positive_among_original_OSR_leaf_samples"],106/2291)

    def test_fake_original_2013_PLOS_data_not_accepted(self):
        txt=CSV.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"false-PLOS.csv"
            f.write_text(txt.replace("2013,279,64,","2013,279,65,"),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"Published original author"):
                read_original_leaf_source_table(f)

    def test_original_PLOS_article_source_requires_true_empirical_sentences(self):
        source="""<html><body>
        10.1371/journal.pone.0114477
        Unexpected Diversity of Feral Genetically Modified Oilseed Rape
        We sampled 279 OSR plants.
        64 out of 279 OSR leaf samples
        106 out of 2291 OSR leaf samples
        Furthermore, 22 OSR leaf samples tested positive for the bar gene.
        Two out of 48 seed pod samples from non-GM OSR mother plants.
        </body></html>"""
        source+=" "*50000
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"plos.html"
            f.write_text(source,encoding="utf-8")
            qc=verify_actual_PLOS_fulltext_source(f)
            self.assertTrue(qc[
                "original_publisher_HTML_downloaded_and_actual_source_number_phrases_verified"])
            f.write_text(source.replace("106 out of 2291","105 out of 2291"),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"exact source biological samples"):
                verify_actual_PLOS_fulltext_source(f)

    def test_2020_gov_source_receipt_needed_to_avoid_overlapping_2013_spliced_cohort(self):
        with self.assertRaisesRegex(ValueError,"Require verified original independent"):
            run("none",CSV,{"complete":False})


if __name__=="__main__":unittest.main()
