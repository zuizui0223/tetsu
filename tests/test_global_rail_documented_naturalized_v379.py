import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import geopandas as gpd
from shapely.geometry import box
from pyproj import Transformer
from analyze_global_rail_documented_naturalized_v379 import cell_region,within_region,RICH_BASE,RICH_BUILT,X_START,Y_START
class TestV379(unittest.TestCase):
    def test_unique_whole_cell_not_center_and_guard(self):
        inv=Transformer.from_crs("EPSG:6933","EPSG:4326",always_xy=True)
        key=(174,74);w=100000
        a,b=inv.transform(X_START+key[0]*w,Y_START+key[1]*w)
        c,d=inv.transform(X_START+(key[0]+1)*w,Y_START+(key[1]+1)*w)
        g=gpd.GeoDataFrame({"region_id":["A"]},geometry=[box(a-.1,b-.1,c+.1,d+.1)],crs=4326)
        self.assertEqual(cell_region(key,w,g,inv,0),"A")
        self.assertIsNone(cell_region(key,w,g,inv,.3))
        g2=gpd.GeoDataFrame({"region_id":["A","B"]},geometry=[box(a-.1,b-.1,c+.1,d+.1)]*2,crs=4326)
        self.assertIsNone(cell_region(key,w,g2,inv,0))
    def test_within_region_removes_between_region_difference(self):
        rows=[]
        for reg,intercept in [("a",.1),("b",.4),("c",.2)]:
            for j in range(20):
                r={k:1. for k in RICH_BASE+(RICH_BUILT,"log_species")}
                r.update(region=reg,continent="Asia",block=(j//2,0),log_rail_km=j/10,
                         known_fraction=intercept+.05*j/10)
                rows.append(r)
        x=within_region(rows)
        self.assertAlmostEqual(x["coefficient_pp_per_log1p_rail_km"],5.,places=8)
        for r in rows:r["log_observations"]=r["log_rail_km"]
        self.assertIsNone(within_region(rows)["coefficient_pp_per_log1p_rail_km"])
    def test_sparse_returns_unsupported(self):
        self.assertIsNone(within_region([])["coefficient_pp_per_log1p_rail_km"])
if __name__=="__main__":unittest.main()
