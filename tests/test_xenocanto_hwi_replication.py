import unittest
from scan_xenocanto_hwi_replication import year_of,binomial,observer,esc,full_avonet_match
class XC(unittest.TestCase):
    def test_date_and_names(self):
        self.assertEqual(year_of("2019-02-01","2019"),2019)
        self.assertIsNone(year_of("2019-02-01","2020"))
        self.assertIsNone(year_of("2019/2020",""))
        self.assertIsNone(observer("A; B"))
        self.assertEqual(observer("  A  B "),observer("a b"))
        self.assertIsNotNone(binomial(lambda k:{"scientificName":"Parus major","taxonRank":"species"}.get(k,"")))
        self.assertIsNone(binomial(lambda k:{"scientificName":"Parus cf major","taxonRank":"species"}.get(k,"")))
        self.assertEqual(esc("\\t"),"\t")
    def test_fixed_hwi_full_source_names(self):
        def row(name,hwi,h="1",m="1"):
            return {"Species1":name,"Hand-Wing.Index":hwi,"Habitat.Density":h,"Migration":m}
        rows=[row("Parus major","20.999"),row("Turdus merula","21"),
              row("Corvus corax","35",h=""),row("Pica pica",""),
              row("Passer domesticus","101")]
        out,audit=full_avonet_match(rows)
        self.assertEqual(out["parus major"],("low_HWI","h1_m1_low"))
        self.assertEqual(out["turdus merula"],("high_HWI","h1_m1_high"))
        self.assertEqual(out["corvus corax"],("high_HWI",))
        self.assertEqual(out["pica pica"],())
        self.assertEqual(out["passer domesticus"],())
        self.assertEqual(audit["HWI_threshold"],21.0)
        self.assertFalse(audit["identical_to_original_iNat_matched_pool"])
        self.assertEqual(audit["mapping_SHA256"],full_avonet_match(rows[::-1])[1]["mapping_SHA256"])
        with self.assertRaises(ValueError):full_avonet_match(rows+[row(" PARUS MAJOR ","30")])
if __name__=="__main__":unittest.main()
