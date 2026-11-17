import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch
MOD=Path(__file__).resolve().parents[1]/'scripts/audit_MIREN_road_reference_scale_v255.py'
spec=importlib.util.spec_from_file_location('m',MOD);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def source():
    regions={}
    for region in ['A','B']:
        roads=[]
        for road in range(2):
            roads.append({'road':f'{region}{road}','sites':4,'raw':[[.3],[.4],[.5]],'equal':[[.6],[.7],[.8]]})
        regions[region]=roads
    return {'regions':regions,'summary':{'raw_jaccard_1minus2':{'n':2,'mean':-.1},'raw_jaccard_1minus3':{'n':2,'mean':-.2},'raw_jaccard_2minus3':{'n':2,'mean':-.1},'equal_jaccard_1minus2':{'n':2,'mean':-.1},'equal_jaccard_1minus3':{'n':2,'mean':-.2},'equal_jaccard_2minus3':{'n':2,'mean':-.1}}}
class Test(unittest.TestCase):
    def test_identity_and_regions(self):
        s=m.make_panel(source(),0,'jaccard')
        self.assertEqual(s['n_regions'],2)
        self.assertEqual(s['n_roads'],4)
        self.assertEqual(s['n_original_sites'],16)
        self.assertAlmostEqual(s['raw']['regional_mean_reference_only_middle_minus_interior'],-.1)
        self.assertEqual(s['equal']['strictly_monotone_regions'],2)
    def test_source_mean_gate(self):
        s=source();s['summary']['equal_jaccard_1minus3']['mean']=-.25
        with self.assertRaises(ValueError):m.make_panel(s,0,'jaccard')
    def test_source_road_sites_gate(self):
        s=source();s['regions']['B'][0]['sites']=3
        with self.assertRaises(ValueError):m.make_panel(s,0,'jaccard')
if __name__=='__main__':unittest.main()
