import importlib.util
import itertools
import pathlib
import unittest

SCRIPT=pathlib.Path(__file__).resolve().parents[1]/"scripts"/"analyze_Kashmir_global_unknown_trait_partial_id_v288.py"
spec=importlib.util.spec_from_file_location("kashmir_partial_v288",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class ExactPartialIdentificationTests(unittest.TestCase):
    def setUp(self):
        self.traits={0:"annual",1:"perennial"}
        # ids 2,3 are unknown globally (not a different life history per plot).
        self.plots={}
        for site in ["Station A","Station B","Between 1","Between 2"]:
            for year in (2014,2017):
                self.plots[site,year,"Track"] = (
                    {0,2} if site=="Station A" else
                    {1,3} if site=="Station B" else
                    {0,2,3} if site=="Between 1" else {1,2})
                self.plots[site,year,"Embankment"]={0,1,2}
                self.plots[site,year,"Interior Plot 1"]= (
                    {1,3} if site!="Between 2" else {0,1,3})
                self.plots[site,year,"Interior Plot 2"]={0,1}

    def brute_force(self,site_values):
        vals=[]
        for x2,x3 in itertools.product((0,1),repeat=2):
            vals.append(sum(base+coefs.get(2,0)*x2+coefs.get(3,0)*x3
                            for base,coefs in site_values)/len(site_values))
        return min(vals),max(vals)

    def test_exact_global_labels_matches_bruteforce(self):
        site_values=[mod.one_site_linear_contrast(
            self.plots,self.traits,4,s,"Track","Interior Plot 1",
            {2014:.5,2017:.5}) for s in mod.common_source_sites(self.plots)]
        ans=mod.partial_identification(site_values)
        actual=self.brute_force(site_values)
        self.assertAlmostEqual(ans["sharp_lower"],actual[0])
        self.assertAlmostEqual(ans["sharp_upper"],actual[1])
        self.assertEqual(ans["source_sites"],4)

    def test_zero_unknowns_reduces_to_identified_point(self):
        classify={0:"annual",1:"perennial",2:"perennial",3:"annual"}
        observed=mod.analyze(self.plots,classify,4)
        ans=observed["comparisons"]["Track__minus__Interior Plot 1"]["all"]["two_year_mean"]
        self.assertAlmostEqual(ans["sharp_lower"],ans["sharp_upper"])
        self.assertEqual(ans["unlabelled_taxa_with_nonzero_global_coefficient"],0)

    def test_missing_original_zone_drops_whole_site(self):
        self.plots["Station A",2017,"Interior Plot 2"]=None
        out=mod.analyze(self.plots,self.traits,4)
        self.assertEqual(out["same_all_zone_two_year_original_sites"],3)
        self.assertEqual(out["station_source_sites"],1)
        self.assertEqual(out["between_source_sites"],2)

    def test_within_site_year_change_uses_same_year_site(self):
        self.plots["Station A",2017,"Track"]={1}
        ans=mod.analyze(self.plots,self.traits,4)
        year14=ans["comparisons"]["Track__minus__Interior Plot 1"]["all"]["2014"]
        year17=ans["comparisons"]["Track__minus__Interior Plot 1"]["all"]["2017"]
        change=ans["comparisons"]["Track__minus__Interior Plot 1"]["all"]["2017_minus_2014"]
        # same scenario with all unknown taxa nonannual
        self.assertAlmostEqual(
            change["all_unlabelled_nonannual_scenario"],
            year17["all_unlabelled_nonannual_scenario"]-year14["all_unlabelled_nonannual_scenario"])

    def test_global_consistency_can_be_narrower_than_independent_unknown_in_every_plot(self):
        plotvalues=[(0,{2:0.5}), (0,{2:-0.5})]
        bounds=mod.partial_identification(plotvalues)
        self.assertEqual(bounds["sharp_lower"],0)
        self.assertEqual(bounds["sharp_upper"],0)

    def test_empty_plot_errors(self):
        with self.assertRaisesRegex(ValueError,"Zero-richness"):
            mod.zone_fractions(set(),self.traits,4)

    def test_all_taxon_indices_validated(self):
        with self.assertRaisesRegex(ValueError,"outside genuine headers"):
            mod.zone_fractions({100},self.traits,4)

    def test_all_strata_same_cohort(self):
        ans=mod.analyze(self.plots,self.traits,4)
        base=ans["comparisons"]["Track__minus__Interior Plot 1"]
        self.assertEqual(base["all"]["2014"]["source_sites"],4)
        self.assertEqual(base["station"]["2014"]["source_sites"],2)
        self.assertEqual(base["between"]["2014"]["source_sites"],2)

if __name__=="__main__":
    unittest.main()
