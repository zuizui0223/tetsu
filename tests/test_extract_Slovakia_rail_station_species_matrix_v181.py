import importlib.util
from pathlib import Path
import unittest

PATH = Path(__file__).resolve().parents[1] / "scripts/extract_Slovakia_rail_station_species_matrix_v181.py"
spec=importlib.util.spec_from_file_location("slovak_v181", PATH)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class FakePage:
    number=4
    def __init__(self, words): self._words=words
    def get_text(self, mode, sort=False):
        if mode!="words": raise ValueError(mode)
        return self._words

def w(x,y,s):
    return (x,y,x+len(s)*4.5,y+10,s,0,0,0)

class TestSlovakiaOriginalQA(unittest.TestCase):
    def test_source_is_pdf_original_not_mock_taxa(self):
        self.assertIn("Májeková",m.__doc__)
        self.assertEqual(len(m.PAGE_INDICES),7)
        self.assertEqual(sum(m.EXPECTED.values()),511)
        self.assertFalse(m.extract.__doc__ or False)
    def test_publisher_headings_positions(self):
        centers=[100,150,220,270]
        a=[]
        for x,s in zip(centers,["historic","recent","historic","recent"]):
            a.append(w(x-13,5,s))
        lines=m.merge_words_to_rows(a)
        out=m.detect_header_centres(lines)
        self.assertEqual(len(out),4)
        self.assertLess(out[0],out[1])
    def test_absence_not_biological_nondetection(self):
        line=[w(10,80,"Echium"),w(47,80,"vulgare"),
              w(100,80,"x"),w(150,80,"b-c"),w(220,80,"x"),w(270,80,"c")]
        page=FakePage(line)
        result,diagnostics=m.parse_table_page(page,[102,153,222,272])
        self.assertEqual(len(result),1)
        self.assertEqual(result[0]["values"],["x","b-c","x","c"])
        self.assertFalse(diagnostics)
    def test_table_gate_not_relaxed(self):
        self.assertEqual(m.UNION_EXPECTED["Kosice"],207)
        self.assertEqual(m.EXPECTED["Spisska_Nova_Ves_recent_2012_2013"],191)

if __name__=="__main__": unittest.main()
