import unittest
from scan_splotopen_stature import year_of,summarize,PANELS
from analyze_splotopen_stature import fit
class Splot(unittest.TestCase):
    def test_date(self):
        self.assertEqual(year_of("05-07-1998"),1998)
        self.assertEqual(year_of("1998-07-05"),1998)
        self.assertIsNone(year_of("1998"))
    def test_species_union_and_gate(self):
        ss={f"s{i}" for i in range(20)}
        life={s:"non-woody" for s in ss};height={s:1. for s in ss}
        sel={str(i):((1,2),2001,1.,"d",False) for i in range(5)}
        by={k:ss.copy() for k in sel}
        out,qa=summarize(sel,by,life,height)
        for p in PANELS:
            self.assertEqual(out[p][0]["known_height"],20)
            self.assertEqual(out[p][0]["n_plots"],5)
            self.assertEqual(out[p][0]["trait_value"],1.)
        del by["0"]
        out,_=summarize(sel,by,life,height)
        self.assertEqual(out[PANELS[0]],[])
        self.assertNotEqual(fit([])["status"],"estimated")
if __name__=="__main__":unittest.main()
