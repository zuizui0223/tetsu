import unittest
from audit_world_NHM_rail_trait_match import norm,number,life_class,match,counts
class TestTraits(unittest.TestCase):
    def test_numeric_strict(self):
        self.assertEqual(number("2.5"),2.5)
        self.assertIsNone(number("2-5"))
        self.assertIsNone(number("0"))
        self.assertIsNone(number(["2"]))
    def test_exact_name_ambiguity(self):
        rows=[{"taxon":"Abc  def","life form":"annual","habit":"herb"},{"taxon":"Abc xyz","life form":"perennial"},{"taxon":"abc XYZ","life form":"annual"}]
        j,qa,_=match(rows,{"abc def":[1],"abc xyz":[2]})
        self.assertEqual(set(j),{1})
        self.assertEqual(j[1]["life_class"],"annual")
        self.assertEqual(qa["source_ambiguous_normalized_names"],1)
        self.assertEqual(counts({1,2},j)["annual"],1)
    def test_conflict(self):
        self.assertIsNone(life_class({"life form":"annual","habit":"tree"}))
if __name__=="__main__":unittest.main()
