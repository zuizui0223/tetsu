"""v319 tests: original RCT published means never synthetic plot replication."""
import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Rivera2014_rail_topsoil_richness_cover_interactions_v319 import value,summarize,parse_original_table2

class OriginalSpainRailTopsoilRCT(unittest.TestCase):
  def fixture(self):
    return {"author_paper_DOI":"10.1371/journal.pone.0101413",
            "original_source_doi_verified":True,
            "true_original_article_sha256":"fixture_only",
            "original_source_tables":[{
              "label":"Table 2",
              "source_author_caption":"Mean values for total plant cover (%) and species richness in spring 2010 and 2011",
              "rows":[
                ["Period","Treatment","Plant cover","Richness","Equitability"],
                ["2010","Control","40.8±3.5 b","7.2±1","0.77±0.04a"],
                ["","T10","21.7±2.0a","6.9±0.8","0.84±0.05ab"],
                ["","T30","17.3±2.1a","6.2±1.3","0.85±0.09b"],
                ["2011","Control","62.0±3.3","14.6±1.63a","0.83±0.03"],
                ["","T10","52.2±3.8","17.8±1.95b","0.81±0.03"],
                ["","T30","61.7±3.2","19.5±1.89b","0.82±0.03"]
              ]}]}
  def test_source_means_and_SE_with_significance_letters(self):
    x=value("17.8±1.95b")
    self.assertAlmostEqual(x["mean"],17.8)
    self.assertAlmostEqual(x["SE_from_original_author"],1.95)
    self.assertEqual(x["original_published_significance_letter"],"b")
    with self.assertRaises(ValueError):value("17.8")
  def test_randomized_2010_11_source_C_vs_T10_T30_means(self):
    d=parse_original_table2(self.fixture())
    self.assertEqual(len(d),6)
    self.assertEqual(d[2011,"T30"]["Richness"]["mean"],19.5)
  def test_descriptive_DID_while_published_p_value_not_new(self):
    j=summarize(self.fixture())
    x=j["outcomes_descriptive_changes_not_independent_plot_effects"]
    self.assertAlmostEqual(x["Plant cover"]["T10"]["descriptive_did_mean_change_minus_control_change"],9.3)
    self.assertAlmostEqual(x["Plant cover"]["T30"]["descriptive_did_mean_change_minus_control_change"],23.2)
    self.assertAlmostEqual(x["Richness"]["T10"]["descriptive_did_mean_change_minus_control_change"],3.5)
    self.assertAlmostEqual(x["Richness"]["T30"]["descriptive_did_mean_change_minus_control_change"],5.9)
    self.assertEqual(j["n_original_treatment_plots"],15)
    self.assertTrue(j["no_source_randomization_p_values_or_new_95_CI_reestimated"])
  def test_missing_rct_treatment_can_never_be_interpolated(self):
    d=self.fixture();d["original_source_tables"][0]["rows"].pop()
    with self.assertRaises(ValueError):summarize(d)
  def test_wrong_source_DOI_cannot_be_analyzed(self):
    d=self.fixture();d["author_paper_DOI"]="10.0000/wrong"
    with self.assertRaises(ValueError):summarize(d)
if __name__=="__main__":unittest.main()
