import importlib.util
import pathlib
import unittest

SCRIPT = pathlib.Path(__file__).resolve().parents[1]/"scripts"/"analyze_Kashmir_global_lifecycle_source_v287.py"
spec=importlib.util.spec_from_file_location("life_cycle_v287",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class GlobalLifeCycleTests(unittest.TestCase):
    def test_header_inference(self):
        s,l,diag=mod.infer_source_columns(["Species_Name","Life_Cycle","Herbaceous"])
        self.assertEqual(s,"Species_Name")
        self.assertEqual(l,"Life_Cycle")

    def test_ambiguous_source_headers_fail_closed(self):
        s,l,diag=mod.infer_source_columns(["Species","Scientific_Name","Life_Cycle"])
        self.assertIsNone(s)

    def test_lowercase_accepted_names_are_still_binomial_species(self):
        donor,_=mod.parse_life_cycles([{"Accepted Name":"poa annua","LifeCycle":"Annual"}],"Accepted Name","LifeCycle")
        self.assertEqual(donor[0]["binomial"],"poa annua")

    def test_categorical_values_are_not_guessed(self):
        rows=[{"Species":"Poa annua","Life_Cycle":"annual"},
              {"Species":"Poa pratensis","Life_Cycle":"perennial"},
              {"Species":"Poa alpina","Life_Cycle":"1"}]
        with self.assertRaisesRegex(ValueError,"manual protocol"):
            mod.parse_life_cycles(rows,"Species","Life_Cycle")

    def test_exact_annual_perennial_are_valid(self):
        rows=[{"Species":"Poa annua","Life_Cycle":"Annual"},
              {"Species":"Poa pratensis","Life_Cycle":"Perennial"}]
        d,info=mod.parse_life_cycles(rows,"Species","Life_Cycle")
        self.assertEqual([x["span"] for x in d],["annual","perennial"])

    def test_csv_source_minimum_size_gate(self):
        with self.assertRaisesRegex(ValueError,"100,000"):
            mod.read_rows(b"Species,Life_Cycle\nPoa annua,Annual\n","test.csv")

if __name__=="__main__":unittest.main()
