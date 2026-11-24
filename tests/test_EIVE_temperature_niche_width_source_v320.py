"""v320 read-only EIVE OOXML schema parser synthetic tests."""
import sys
import unittest
from pathlib import Path
from xml.etree import ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import audit_EIVE_temperature_niche_width_source_v320 as m

class OoxmlChecks(unittest.TestCase):
    def test_cell_numbering(self):
        self.assertEqual([m.colnum(x) for x in ("A1","Z2","AA3","AB4")],[0,25,26,27])
    def test_bad_cell_address_rejected(self):
        with self.assertRaises(ValueError):m.colnum("1A")
    def test_inline_and_shared_string(self):
        uri=m.NS["s"]
        a=ET.fromstring(f'<c xmlns="{uri}" r="C1" t="s"><v>0</v></c>')
        b=ET.fromstring(f'<c xmlns="{uri}" r="C2" t="inlineStr"><is><t>Plant</t></is></c>')
        self.assertEqual(m.cell_value(a,["A species"]),"A species")
        self.assertEqual(m.cell_value(b,[]),"Plant")
    def test_error_cell_rejected(self):
        uri=m.NS["s"]
        a=ET.fromstring(f'<c xmlns="{uri}" r="D4" t="e"><v>#VALUE!</v></c>')
        with self.assertRaises(ValueError):m.cell_value(a,[])
if __name__=="__main__":unittest.main()
