import unittest
from scan_xenocanto_hwi_replication import year_of,binomial,observer,esc
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
if __name__=="__main__":unittest.main()
