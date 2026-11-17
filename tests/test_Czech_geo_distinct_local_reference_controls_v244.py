import math,sys,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"scripts"))
import validate_Czech_rail_station_local_reference_geo_cells_unique_v244 as mod
from validate_Czech_rail_station_local_reference_geo_cells_unique_v244 import (
    approximate_coord_cell,exact_max_stations_with_cell_uniqueness
)
class TestExactOriginalCzechCoordinateDisjointness(unittest.TestCase):
    def data(self,overlap=False):
        sites=[];geo={}
        for i in range(16):
            ref=[]
            for p in range(8):
                # Source GPS 4decimal cells distinct within station but
                # the same original 3decimal rounded cell can occur
                # multiple times for dense field plots.
                lat=48.45+i*.052+p*.00013
                lon=13.01+i*.087+p*.00013
                if overlap:lat=48.45+i*.052;lon=13.01+i*.087
                rid=f"source_{i}_original_plot_{p}"
                geo[rid]={'lat':lat,'lon':lon,'year':2022,'study':'a' if p<5 else 'b'}
                ref.append({'id':rid,'study':'a' if p<5 else 'b',
                            'd':.7+p*.21,'taxa':{'Genus speciesA'}})
            sites.append({'id':f'rail_{i}','stretch':f'orig_stretch_{i%4}',
                          'control_options':ref,
                          'rail_source_trait_taxa':{'Genus speciesA'}})
        return sites,geo
    def test_source_grid_bin_4decimal_and_3decimal_real_WGS84(self):
        sites,geo=self.data()
        first=sites[0]['control_options'][0]
        a=approximate_coord_cell(first,geo,4)
        b=approximate_coord_cell(first,geo,3)
        self.assertEqual(a,(48.45,13.01))
        self.assertEqual(b,(48.45,13.01))
        self.assertNotEqual(approximate_coord_cell(sites[0]['control_options'][1],geo,4),a)
    def test_original_global_source_cell_constraint_allows_true_unique_synthetic_coordinates(self):
        sites,geo=self.data()
        with patch.object(mod,'SOURCE_KNOWN_MAX_STATIONS',16):
            alloc,qa=exact_max_stations_with_cell_uniqueness(sites,geo,4,timeout_s=12)
        self.assertIsNotNone(alloc)
        self.assertEqual(qa['status'],'EXACT_SOURCE_MILP_GEO_UNIQUE_OPTIMAL')
        self.assertEqual(qa['n_source_independent_releve_plot_IDs'],5*len(alloc))
        used=[approximate_coord_cell(p,geo,4) for rows in alloc.values() for p in rows]
        self.assertEqual(len(used),len(set(used)))
        self.assertEqual(len(alloc),16)
    def test_all_same_field_position_within_station_cannot_count_five_distinct_plots(self):
        sites,geo=self.data(overlap=True)
        with patch.object(mod,'SOURCE_KNOWN_MAX_STATIONS',16):
            alloc,qa=exact_max_stations_with_cell_uniqueness(sites,geo,4,timeout_s=12)
        self.assertIsNone(alloc)
        self.assertEqual(qa['status'],'MIP_NOT_CERTIFIED_STOP')
    def test_negative_or_bad_field_original_source_geolocation_rejected(self):
        sites,geo=self.data()
        geo[sites[0]['control_options'][0]['id']]['lat']=0
        with self.assertRaises(ValueError):approximate_coord_cell(sites[0]['control_options'][0],geo,4)
if __name__=="__main__":unittest.main()
