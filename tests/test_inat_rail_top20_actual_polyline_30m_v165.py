"""v16.5 synthetic GIS/API metadata-only test fixtures: no photo or person IDs
are real data and no invented plant count can be published as observed.
"""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_20_original_rail_plants_inat_30m_two_regions_v165 import (
    validate_source_candidates,original_pbf_authentic,exact_taxon_or_descendant,
    observation_quality,load_complete_photo_observations,fixed_query,
    summarize_species_geo,ACCURACY_MAX_M,RAIL_DISTANCE_M,
    RailwayPolylineSource,STUDY_REGIONS
)

def fake_observed_record(species="Poa annua",observed_date="2024-05-21",accuracy=5,lon=14.,lat=48.):
    return {"id":13,"taxon":{"name":species},"obscured":False,"geoprivacy":"open",
      "captive":False,"positional_accuracy":accuracy,"observed_on":observed_date,
      "geojson":{"coordinates":[lon,lat]},
      "user":{"id":17},"photos":[{"license_code":"cc-by"}],
      "annotations":[{"controlled_attribute_id":12,"controlled_value_id":13}]}

def fake_rail_at(lon=14.,lat=48.):
    from shapely.geometry import LineString
    from shapely.strtree import STRtree
    from pyproj import Transformer
    obj=object.__new__(RailwayPolylineSource)
    obj.transformer=Transformer.from_crs(4326,32633,always_xy=True)
    x,y=obj.transformer.transform(lon,lat)
    obj.lines=[LineString([(x-100,y),(x+100,y)])]
    obj.tree=STRtree(obj.lines)
    return obj

class INatGeoPolyline30mTests(unittest.TestCase):
    def test_nearest_polyline_distance_not_nearest_osm_way_vertex_distance(self):
        rail=fake_rail_at()
        d,cell=rail.distance_and_cell(14.,48.)
        self.assertAlmostEqual(d,0.,places=6)
        self.assertEqual(len(cell),2)
        from pyproj import Transformer
        inv=Transformer.from_crs(32633,4326,always_xy=True)
        x,y=rail.transformer.transform(14.,48.)
        nlon,nlat=inv.transform(x,y+22)
        d,_=rail.distance_and_cell(nlon,nlat)
        self.assertAlmostEqual(d,22.,places=5)
        self.assertLessEqual(d,30.)

    def test_quality_filter_positional_accuracy_and_obscuration(self):
        original=fake_observed_record()
        self.assertEqual(observation_quality(original,"Poa annua")[1],"pass")
        o=fake_observed_record(accuracy=31)
        self.assertEqual(observation_quality(o,"Poa annua")[1],"missing_or_above_30m_public_coordinate_accuracy")
        o=fake_observed_record()
        o["geoprivacy"]="obscured"
        self.assertEqual(observation_quality(o,"Poa annua")[1],"obscured_or_private")
        o=fake_observed_record();o["captive"]=True
        self.assertEqual(observation_quality(o,"Poa annua")[1],"captive_or_cultivated")
        o=fake_observed_record();o["photos"]=[]
        self.assertEqual(observation_quality(o,"Poa annua")[1],"no_photos")
        o=fake_observed_record();o["geojson"]["coordinates"]=[0.,float("nan")]
        self.assertEqual(observation_quality(o,"Poa annua")[1],"missing_or_invalid_public_coordinates")

    def test_only_exact_original_taxon_or_curated_accepted_synonym(self):
        self.assertTrue(exact_taxon_or_descendant("Poa annua","Poa annua"))
        self.assertTrue(exact_taxon_or_descendant("Erigeron canadensis","Conyza canadensis"))
        self.assertTrue(exact_taxon_or_descendant("Poa annua subsp. annua","Poa annua"))
        self.assertFalse(exact_taxon_or_descendant("Poa trivialis","Poa annua"))
        self.assertFalse(exact_taxon_or_descendant("Erigeron annuus","Conyza canadensis"))
        self.assertFalse(exact_taxon_or_descendant("Poa","Poa annua"))

    def test_certain_within30_requires_coordinate_error_circle_inside(self):
        rail=fake_rail_at()
        q=[fake_observed_record(accuracy=10),fake_observed_record(accuracy=30),
           fake_observed_record(accuracy=5,species="Poa trivialis")]
        result=summarize_species_geo(q,"Poa annua",rail)
        self.assertEqual(result["n_point_near_30m_strict_accuracy"],2)
        self.assertEqual(result["n_location_certain_within_30m"],2)
        self.assertEqual(result["rejected_original_photo_quality_reasons"]["accepted_taxon_mismatch_unresolved_synonym_or_complex"],1)
        self.assertTrue(result["whole_plant_height_photo_validation_not_performed"])

    def test_unobservable_total_over_API_medium_limit_is_unknown_not_zero(self):
        class FakeClient:
            def get(self,params):
                return {"total_results":9000,"results":[]}
        meta,rows=load_complete_photo_observations(FakeClient(),"Poa annua",STUDY_REGIONS["CZ"]["bbox"])
        self.assertFalse(meta["complete"])
        self.assertEqual(meta["total"],9000)
        self.assertIsNone(rows)
        self.assertIsNone(summarize_species_geo(None,"Poa annua",None)["n_point_near_30m_strict_accuracy"])

    def test_no_fake_20_species_source_if_rail_occupancy_not_verified(self):
        with self.assertRaisesRegex(ValueError,"unverified"):
            validate_source_candidates({"version":"fake","original_selected_20_species":[]})

    def test_exact_source_geofabrik_original_hash_required(self):
        import tempfile
        with tempfile.TemporaryDirectory() as p:
            f=Path(p)/"fake2023.osm.pbf";f.write_bytes(b"not authentic")
            with self.assertRaisesRegex(ValueError,"size changed"):
                original_pbf_authentic(f,"CZ")

    def test_full_download_pages_source_deduplicated_with_true_total(self):
        class FakeClient:
            def get(self,p):
                if p["page"]==1:return {"total_results":2,
                     "results":[fake_observed_record()]}
                raise AssertionError("Only one page requested; intentionally incomplete census must fail")
        with self.assertRaisesRegex(ValueError,"Incomplete public"):
            load_complete_photo_observations(FakeClient(),"Poa annua",STUDY_REGIONS["SI"]["bbox"])

if __name__=="__main__":unittest.main()
