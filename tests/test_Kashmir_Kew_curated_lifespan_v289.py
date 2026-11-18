import csv
import importlib.util
import pathlib
import tempfile
import unittest

SCRIPT=pathlib.Path(__file__).resolve().parents[1]/"scripts"/"analyze_Kashmir_Kew_curated_strict_annual_v289.py"
spec=importlib.util.spec_from_file_location("KashmirKewCuratedV289",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class VerifiedKewAnnualTests(unittest.TestCase):
    def setUp(self):
        self.names=["Poa annua","Cynodon dactylon (L.) Pers.","Breea arvensis (L.) Less."]
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.global_csv=pathlib.Path(self.temp.name)/"global.csv"
        self.curated_csv=pathlib.Path(self.temp.name)/"kew.csv"
        with self.global_csv.open("w",newline="") as fp:
            w=csv.DictWriter(fp,fieldnames=["original_taxon","life_span","trait_source"])
            w.writeheader()
            w.writerow({"original_taxon":"Poa annua","life_span":"annual","trait_source":"Global2023"})
        with self.curated_csv.open("w",newline="") as fp:
            fields=["original_taxon","life_span","accepted_name","synonym_evidence_url",
                    "life_span_evidence_url","trait_source"]
            w=csv.DictWriter(fp,fieldnames=fields)
            w.writeheader()
            for x,accepted in [
                 ("Cynodon dactylon (L.) Pers.","Cynodon dactylon (L.) Pers."),
                 ("Breea arvensis (L.) Less.","Cirsium arvense (L.) Scop.")]:
                w.writerow({"original_taxon":x,"life_span":"perennial",
                    "accepted_name":accepted,
                    "synonym_evidence_url":"https://powo.science.kew.org/taxon/names",
                    "life_span_evidence_url":"https://powo.science.kew.org/taxon/names",
                    "trait_source":"TestKew"})

    def test_verified_two_additions_only(self):
        d,m=mod.curated_merge(self.names,self.global_csv,self.curated_csv)
        self.assertEqual(len(d),3)
        self.assertEqual([x["life_span"] for x in d],["annual","perennial","perennial"])
        self.assertEqual(m["new_independent_accepted_name_taxa"],2)

    def test_no_overwrite_of_global_trait_labels(self):
        with self.global_csv.open("a") as fp:
            fp.write("Cynodon dactylon (L.) Pers.,annual,ConflictingSource\n")
        with self.assertRaisesRegex(ValueError,"Never overwrite"):
            mod.curated_merge(self.names,self.global_csv,self.curated_csv)

    def test_missing_scientific_provenance_rejected(self):
        data=self.curated_csv.read_text()
        data=data.replace("https://powo.science.kew.org/taxon/names","https://nonbotanical.example/names")
        self.curated_csv.write_text(data)
        with self.assertRaisesRegex(ValueError,"Kew scientific"):
            mod.curated_merge(self.names,self.global_csv,self.curated_csv)

    def test_analytic_partial_id_upper_bound_decreases_when_positive_unknown_is_confirmed_perennial(self):
        baseline=mod.partial_identification([(0.0,{1:.14,2:-.04})])
        verified=mod.partial_identification([(0.0,{2:-.04})])
        self.assertAlmostEqual(baseline["sharp_upper"],.14)
        self.assertAlmostEqual(verified["sharp_upper"],0.0)

    def test_source_cluster_bootstrap_is_reproducible(self):
        p={}
        for site in ["Station A","Station B","Between 1"]:
            for year in (2014,2017):
                p[(site,year,"Track")]={1} if site=="Station A" else {0,1}
                p[(site,year,"Embankment")]={0,1}
                p[(site,year,"Interior Plot 1")]={0,2}
                p[(site,year,"Interior Plot 2")]={0,1,2}
        a=mod.technical_source_bootstrap(p,{0:"annual",1:"perennial"},3,
                                         list(mod.common_source_sites(p)),
                                         "Track","Interior Plot 1",100,seed=5)
        b=mod.technical_source_bootstrap(p,{0:"annual",1:"perennial"},3,
                                         list(mod.common_source_sites(p)),
                                         "Track","Interior Plot 1",100,seed=5)
        self.assertEqual(a,b)
        self.assertGreaterEqual(a["fraction_bootstrap_sharp_upper_negative"],0)
        self.assertLessEqual(a["fraction_bootstrap_sharp_upper_negative"],1)

if __name__=="__main__":
    unittest.main()
