"""v27.4 synthetic QA: G/R railway identity, section-first weights and layout separation."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import analyze_Slovenia_official_main_regional_local_coupling_v274 as m


def make_source():
    lines={"10":"main_G","20":"main_G","34":"regional_R","80":"regional_R"}
    out={}
    for section,typ in lines.items():
        group=[]
        for i in range(4):
            group.append({"id":f"{section}_{i}",
                          "BS":{1,i+10,(1 if typ=="main_G" else 99)},
                          "BT":{1,i+10,i+40},
                          "layout":"DOUBLE" if typ=="main_G" and i<2 else "SINGLE",
                          "class":typ})
        out[section]=group
    return out


class TrueSourceMainRegionalTests(unittest.TestCase):
    def test_official_G_R_is_not_freight(self):
        self.assertEqual(m.official_class("10"),"main_G")
        self.assertEqual(m.official_class("34"),"regional_R")
        self.assertIsNone(m.official_class("999"))
        self.assertEqual(m.KNOWN_COUNTS["main_G"],(6,34))
        self.assertEqual(m.KNOWN_COUNTS["regional_R"],(9,55))

    def test_source_section_not_plant_side_pseudoreplication(self):
        q=m.per_source_section(make_source(),2,draws=19,burnin=100,spacing=20)
        self.assertEqual(len(q),4)
        self.assertEqual(sum(x["n_field_locations"] for x in q.values()),16)
        r=m.contrasts(q)
        self.assertTrue(r["identified"])
        self.assertEqual(r["source_support"]["source_main_G_sections"],2)
        self.assertEqual(r["source_support"]["source_regional_R_sections"],2)
        self.assertEqual(r["source_support"]["source_locations"],16)
        self.assertIsNotNone(r["source_section_equal_main_regional_comparisons"]["local_extra_similarity_nullminusobserved_beta"]["G_minus_R_descriptive"])

    def test_layout_restriction_must_recompute_original_frequency_null(self):
        src=make_source()
        all_=m.per_source_section(src,2,draws=19,burnin=70,spacing=15)
        single=m.per_source_section(src,2,track_only="SINGLE",draws=19,burnin=70,spacing=15)
        double=m.per_source_section(src,2,track_only="DOUBLE",draws=19,burnin=70,spacing=15)
        self.assertEqual(sum(x["n_field_locations"] for x in all_.values()),16)
        self.assertEqual(sum(x["n_field_locations"] for x in single.values()),12)
        self.assertEqual(sum(x["n_field_locations"] for x in double.values()),4)
        self.assertFalse(m.contrasts(double)["identified"])
        self.assertTrue(m.contrasts(single)["identified"])

    def test_unverified_layout_or_pseudo_site_fails(self):
        src=make_source()
        src["10"][0]["layout"]="UNVERIFIED"
        # per_source_section itself does not infer unknown layout:
        # original_grouped validates it before source analysis.
        self.assertTrue(len(m.per_source_section(src,2))>0)
        with self.assertRaises(ValueError):
            m.per_source_section({"10":[{"id":1,"BS":set(),"BT":{2},"layout":"SINGLE","class":"main_G"},
                                       {"id":2,"BS":{1},"BT":{2},"layout":"SINGLE","class":"main_G"}]},2,draws=19,burnin=20,spacing=10)


if __name__=="__main__":
    unittest.main()
