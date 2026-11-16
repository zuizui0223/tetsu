"""v11.8: Pladias author source file must be actual parsed species traits, not login HTML."""
import csv,io,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_Pladias_public_diaspore_lifeform_v118 import (
    csv_rows,inspect_source,run,DISP_CODES
)


class PublicPladiasSameSourceTraitTests(unittest.TestCase):
    def test_original_strategy_archetypes_identical_source_labels(self):
        self.assertEqual(len(DISP_CODES),9)
        self.assertIn("Epilobium",DISP_CODES)
        self.assertIn("Phragmites",DISP_CODES)
        self.assertIn("Wolffia",DISP_CODES)

    def test_html_login_page_never_treated_as_biology(self):
        with self.assertRaisesRegex(ValueError,"HTML"):
            csv_rows(("<html>"+("Something\n"*100)+"</html>").encode("utf8"))

    def test_source_disjunct_files_and_trait_strategies(self):
        def src(trait):
            r=io.StringIO()
            writer=csv.writer(r,delimiter="\t")
            writer.writerow(["taxon","trait"])
            for i in range(60):
                writer.writerow([f"Plant{i}",("Epilobium" if i%2 else "Allium") if trait=="strategy"
                                 else ("therophyte" if i%2 else "hemicryptophyte")])
            return r.getvalue().encode("utf8")
        with tempfile.TemporaryDirectory() as d:
            a=Path(d)/"disp.tsv";b=Path(d)/"life.tsv"
            a.write_bytes(src("strategy"));b.write_bytes(src("life"))
            x=run(a,b)
        self.assertTrue(x["complete"])
        self.assertEqual(x["strategy_original_source_metadata"]["n_original_species_trait_source_rows"],60)
        self.assertEqual(x["strategy_original_source_metadata"]["strategy_archetype_token_counts_all_fields"]["Epilobium"],30)
        self.assertTrue(x["no_Czech_rail_trait_effect_fitted"])

    def test_same_file_not_two_independent_traits(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"source.tsv"
            w=io.StringIO();writer=csv.writer(w,delimiter="\t")
            writer.writerow(["species","strategy"])
            for i in range(70):writer.writerow([str(i),"Allium"])
            path.write_bytes(w.getvalue().encode("utf8"))
            with self.assertRaisesRegex(ValueError,"Same published traits file"):
                run(path,path)

if __name__=="__main__":unittest.main()
