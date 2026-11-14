"""v12.6: only original published CzechVeg class prefix if anchored in original source."""
from __future__ import annotations
import sys,unittest,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_CzechVeg_original_vegetation_class_v126 import (
    source_canonical_syntaxon_class,convert_source_syntaxon_meta,
    matched_class_summary,canonical_roots_from_original_full_publisher_csv,run,
)

def raw(arm,study="real paper",region="18a",year=2012,area=16,syntax="XDA01"):
    return {
        "strict":arm,"broad":arm,"study":study,
        "source_region":region,"survey_year":year,"area":area,
        "SYNTAXON":syntax
    }


class OriginalSourceCanonicalClassTests(unittest.TestCase):
    def test_truly_published_class_to_alliance_association_prefix(self):
        original={"XD","TD","XB"}
        self.assertEqual(source_canonical_syntaxon_class("XD",original),"XD")
        self.assertEqual(source_canonical_syntaxon_class("XDA",original),"XD")
        self.assertEqual(source_canonical_syntaxon_class("XDA01",original),"XD")
        self.assertEqual(source_canonical_syntaxon_class("TDA",original),"TD")
        self.assertEqual(source_canonical_syntaxon_class("TDA01a",original),"TD")

    def test_unknown_or_unanchored_codes_never_guess_source_class(self):
        for x in (None,"","NA","UNKNOWN","X", "XYA01", "XD_A01",
                  "xd", "XDB01st", "12345"):
            self.assertIsNone(source_canonical_syntaxon_class(x,{"XD","TD","XB"}))

    def test_exact_syntaxon_different_but_same_official_class_is_valid_control(self):
        src={
            "a":raw("rail",syntax="XDA01"),
            "b":raw("road",syntax="XDB02"),
            "c":raw("rail",syntax="XDC03"),
            "d":raw("road",syntax="XD"),
            "e":raw("road",syntax="TDA01"),
        }
        x,q=convert_source_syntaxon_meta(src,{"XD","TD"})
        base=matched_class_summary(x,"strict",
             "same_study_region_decade_area_class",2)
        self.assertEqual(base["n_original_same_source_study_vegetation_class_comparison_blocks"],1)
        self.assertEqual(base["n_original_rail_plots_on_class_comparable_support"],2)
        self.assertEqual(base["n_original_road_plots_on_class_comparable_support"],2)
        self.assertEqual(base["matched_source_class_codes"],["XD"])
        self.assertEqual(q["candidate_with_verified_original_source_class_code"],5)

    def test_original_full_metadata_NA_is_not_a_vegetation_class_root(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/"plot_source.csv"
            with f.open("w",encoding="utf-8",newline="") as o:
                o.write("RELEVE_NR,SYNTAXON\n")
                for i in range(117739):
                    o.write(f"{i},{'XD' if i==0 else 'XDA01' if i==1 else 'NA'}\n")
            roots,quality=canonical_roots_from_original_full_publisher_csv(f)
        self.assertEqual(roots,{"XD"})
        self.assertEqual(quality["n_independent_source_class_code_root_labels"],1)
        self.assertEqual(quality["n_source_syntaxon_nonmissing_all_original"],2)

    def test_original_publication_year_area_region_not_silently_pooled(self):
        src={
            "a":raw("rail"),
            "b":raw("road",area=64),
            "c":raw("road",year=1999),
            "d":raw("road",region="68"),
            "e":raw("road",study="a different source"),
        }
        x,q=convert_source_syntaxon_meta(src,{"XD"})
        tight=matched_class_summary(x,"strict",
               "same_study_region_decade_area_class",1)
        self.assertEqual(tight["n_original_same_source_study_vegetation_class_comparison_blocks"],0)
        loose=matched_class_summary(x,"strict",
               "same_study_region_class",1)
        self.assertEqual(loose["n_original_same_source_study_vegetation_class_comparison_blocks"],1)

    def test_source_genuine_publication_receipts_required(self):
        with self.assertRaisesRegex(ValueError,"v12.4"):
            run({},"bad",{},{"complete":False})


if __name__=="__main__":unittest.main()
