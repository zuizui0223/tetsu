import random,sys,unittest,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'scripts'))
from validate_Czech_39station_area_matched_functional_divergence_v237 import (
  area_bands,source_margins,means_and_LOO,shuffle_null
)
class TestCzechOriginalAreaStratifiedNull(unittest.TestCase):
    def setUp(self):
        rng=random.Random(20261116)
        names=[f't{i:03d}' for i in range(120)]
        self.values={x:(rng.gauss(0,1),rng.gauss(0,1)) for x in names}
        sizes=(5,7,7,9,11)
        self.sites=[];n=0
        for m,num in enumerate(sizes):
            for local in range(num):
                self.sites.append({
                   'station':f'st{n:02d}','stretch':f'section{m}',
                   'taxa':set(rng.sample(names,30)),
                   'source_unequal_station_area_ha':.2+local*.35+m*.03})
                n+=1
    def test_five_actual_stretches_10_area_bins_min2(self):
        bands,qc=area_bands(self.sites)
        self.assertEqual(len(bands),10)
        self.assertEqual(len(qc),5)
        self.assertEqual(sum(len(v) for v in bands.values()),39)
        self.assertTrue(all(len(v)>=2 for v in bands.values()))
    def test_original_taxon_margins_preserved_under_area_groups(self):
        bands,_=area_bands(self.sites)
        before=source_margins(self.sites,bands)
        null,chains,qc=shuffle_null(self.sites,self.values,bands,
                                   draws=2,burnin=4500,gap=1000)
        self.assertEqual(len(null),4)
        self.assertEqual(before,source_margins(self.sites,bands))
        self.assertTrue(all(math.isfinite(q['all_39station_pair_functional_centroid_distance']) for q in null))
    def test_independent_stretch_drop_changes_site_pairs(self):
        q=means_and_LOO(self.sites,self.values)
        self.assertEqual(len(q['leave_one_source_stretch_out_beta']),5)
        self.assertGreater(q['all_39station_pair_functional_centroid_distance'],0)
        self.assertTrue(all(x>0 for x in q['leave_one_source_stretch_out_beta'].values()))
    def test_station_area_identical_bin_boundary_rejected(self):
        # Source stretch 0 has five stations and a 2/3 rank split.
        self.sites[1]['source_unequal_station_area_ha']=self.sites[2]['source_unequal_station_area_ha']
        with self.assertRaises(ValueError):area_bands(self.sites)
    def test_invalid_station_area_never_assigned_control_zero(self):
        self.sites[4]['source_unequal_station_area_ha']=None
        with self.assertRaises(ValueError):area_bands(self.sites)
if __name__=='__main__':unittest.main()
