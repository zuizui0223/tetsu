import unittest
from pyproj import Transformer
from scan_global_grip4_road_crossings import feature_edges,W
from analyze_global_hwi_road_crossings import fit_road,ACTIVE,INACTIVE
from collections import Counter
from scan_fullworld_2025_iNat_plant_macro_grid_v369 import X_START,Y_START
class RoadTests(unittest.TestCase):
    def test_actual_crossing(self):
        inv=Transformer.from_crs(6933,4326,always_xy=True)
        ident=Transformer.from_crs(4326,4326,always_xy=True)
        fw=Transformer.from_crs(4326,6933,always_xy=True)
        xx=X_START+87*W; yy=Y_START+37*W+W/2
        pts=[inv.transform(xx-1000,yy),inv.transform(xx+1000,yy)]
        g={"type":"LineString","coordinates":pts}
        self.assertEqual(feature_edges(g,ident,fw,Counter()),{((86,37),(87,37))})
        self.assertFalse(feature_edges({"type":"LineString","coordinates":pts[:1]+[inv.transform(xx-2000,yy)]},ident,fw,Counter()))
    def test_paired_direct_contrast(self):
        rows=[]
        for block in range(25):
            for i in range(8):
                a=i&1;b=(i>>1)&1;c=(i>>2)&1
                rows.append({"a":(block,i*2),"b":(block,i*2+1),"continent":"test","block":block,ACTIVE:a,INACTIVE:b,"cross_highway":c,"y":.03*a-.01*b+.005*c+block*.001})
        r=fit_road(rows,"y",False,(ACTIVE,INACTIVE,"cross_highway"),"cross_highway")
        self.assertEqual(r["status"],"estimated")
        self.assertAlmostEqual(r["active_minus_road"]["estimate"],.025)
        self.assertAlmostEqual(r["active_minus_inactive"]["estimate"],.04)
        for row in rows:row["cross_highway"]=1
        self.assertNotEqual(fit_road(rows,"y",False,(ACTIVE,INACTIVE,"cross_highway"),"cross_highway")["status"],"estimated")
if __name__=="__main__":unittest.main()
