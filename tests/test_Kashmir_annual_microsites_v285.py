import importlib.util
import pathlib
import tempfile
import unittest

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "analyze_Kashmir_annual_microsites_v285.py"
spec = importlib.util.spec_from_file_location("kashmir_annual_v285", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class AnnualMicrohabitatTests(unittest.TestCase):
    def setUp(self):
        self.traits = {0:"annual",1:"perennial",2:"biennial"}
        self.plots = {}
        for site in ("Station A", "Station B", "Between 1", "Between 2"):
            for year in (2014,2017):
                self.plots[site, year, "Track"]={0,3}
                self.plots[site, year, "Embankment"]={0,1}
                self.plots[site, year, "Interior Plot 1"]={1,2}
                self.plots[site, year, "Interior Plot 2"]={0,1,2}

    def test_strict_annual_known_denominator_and_unknown_bounds(self):
        a = module.annual_stat({0,1,3},self.traits)
        self.assertEqual(a["annual_fraction"],.5)
        self.assertEqual(a["minimum_fraction"],1/3)
        self.assertEqual(a["maximum_fraction"],2/3)

    def test_strict_coverage_gate_excludes_unknown_track(self):
        r = module.summarize(self.plots,self.traits,min_coverage=.7,bootstrap=20)
        self.assertEqual(r["complete_four_zones_both_years_sites"],4)
        self.assertEqual(r["comparisons"]["Track__minus__Interior Plot 1"]["all"]["eligible_original_sites"],0)

    def test_site_level_aggregation_two_years_count_once(self):
        r = module.summarize(self.plots,self.traits,min_coverage=.5,bootstrap=50)
        x = r["comparisons"]["Track__minus__Interior Plot 1"]["all"]
        self.assertEqual(x["eligible_original_sites"],4)
        self.assertAlmostEqual(x["two_year_mean"],1.0)
        self.assertEqual(x["sites_positive"],4)
        self.assertEqual(x["per_year"],{"2014":1.,"2017":1.})
        self.assertEqual(x["site_cluster_bootstrap_technical_percentiles_not_field_CI"],[1.,1.])

    def test_full_zone_common_site_not_imputed(self):
        self.plots["Station A", 2017, "Interior Plot 2"] = None
        r=module.summarize(self.plots,self.traits,min_coverage=.5,bootstrap=0)
        self.assertEqual(r["complete_four_zones_both_years_sites"],3)
        self.assertEqual(r["site_group"],{"station":1,"between":2})

    def test_unmatched_trait_name_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p=pathlib.Path(d)/"traits.csv"
            p.write_text("original_taxon,life_span,trait_source\nNot in source,annual,Test\n",encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"unmatched"):
                module.life_span_crosswalk(p,["Species A","Species B"])

    def test_missing_provenance_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p=pathlib.Path(d)/"traits.csv"
            p.write_text("original_taxon,life_span,trait_source\nSpecies A,annual,\n",encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"trait_source"):
                module.life_span_crosswalk(p,["Species A"])

if __name__ == "__main__":
    unittest.main()
