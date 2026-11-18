import importlib.util
import pathlib
import unittest

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "analyze_Kashmir_external_alpine_lifespan_v286.py"
spec = importlib.util.spec_from_file_location("kashmir_external_v286", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class ExternalTraitCrosswalkTests(unittest.TestCase):
    def test_unambiguous_exact_and_binomial(self):
        orig=["Poa annua","Capsella bursa-pastoris L.","Bromus sp."]
        donor=[{"name":"Poa annua","binomial":"poa annua","synonym":"","span":"annual"},
               {"name":"Capsella bursa-pastoris","binomial":"capsella bursa-pastoris","synonym":"","span":"annual"}]
        out,meta=mod.crosswalk(orig,donor,True)
        self.assertEqual(len(out),2)
        self.assertEqual(meta["match_tiers"],{"exact_source_name":1,"unambiguous_species_binomial":1})

    def test_contradictory_donor_life_spans_rejected(self):
        out,_=mod.crosswalk(["Poa annua"],
            [{"name":"Poa annua","binomial":"poa annua","span":"annual"},
             {"name":"Poa annua","binomial":"poa annua","span":"perennial"}],True)
        self.assertEqual(out,[])

    def test_rank_mismatch_and_genus_alone_not_matched(self):
        donor=[{"name":"Poa annua","binomial":"poa annua","span":"annual"}]
        out,_=mod.crosswalk(["Poa sp.","Poa pratensis"],donor,True)
        self.assertEqual(out,[])

    def test_external_table_from_rows(self):
        rows=[
           ["title"],
           ["Scientific Name","Common Synonym","Family","Life Span"],
        ]+[["Species%d alpha"%i,"", "Family", "Annual" if i%2 else "Perennial"] for i in range(159)]
        taxa,meta=mod.detect_life_span_table(rows)
        self.assertEqual(len(taxa),159)
        self.assertEqual(meta["external_lifespan_counts"],{"perennial":80,"annual":79})

    def test_unknown_life_span_not_invented(self):
        rows=[["Scientific Name","Life Span"]]+[["Poa annua","Unknown"] for _ in range(159)]
        with self.assertRaisesRegex(ValueError,"Non-standard"):
            mod.detect_life_span_table(rows)

    def test_species_binomial_with_author_is_explicitly_not_fuzzy(self):
        self.assertEqual(mod.binomial("Capsella bursa-pastoris (L.) Medik."),
                         "capsella bursa-pastoris")
        self.assertIsNone(mod.binomial("Unknown Grass"))
        self.assertIsNone(mod.binomial("Poa spp."))

if __name__=="__main__":
    unittest.main()
