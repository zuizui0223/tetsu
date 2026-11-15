"""v15.7 synthetic test fixtures are NEVER actual Preslia botanical observations."""
import sys, unittest, math, tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from analyze_Czech_2024_real_archaeophyte_neophyte_geographic_identity_v157 import (
    read_source_status,group_sets,source_class_pairs,source_subgroup_gap,
    conditional_group_null
)
from analyze_Czech_2024_real_35station_centroid_geography_identity_v154 import (
    conditional_source_gap
)

def artificial():
    sites=[]
    status={}
    for i in range(35):
        taxa=set()
        for j in range(5):
            name="arch_"+str((i+j)%9);taxa.add(name);status[name]="arch"
            name="neo_"+str((i*3+j)%10);taxa.add(name);status[name]="neo"
        taxa.add("unlabelled_"+str(i))
        status["unlabelled_"+str(i)]=None
        sites.append({"section":str(i%5),"site":str(i),"all_original_taxa":taxa})
    pairs=[]
    for i in range(35):
        for j in range(i+1,35):
            pairs.append({"i":i,"j":j,"same_stretch":sites[i]["section"]==sites[j]["section"],
                          "geo_bin":2,"function_bin":3,"area_log_diff_bin":0})
    return sites,status,pairs

class CzechAlienOriginalStudyIntegrity(unittest.TestCase):
    def test_original_explicit_archaeophyte_neophyte_are_not_native(self):
        sites,status,pairs=artificial()
        d=group_sets(sites,status)
        self.assertEqual(set(d),{"archaeophyte","neophyte","explicit_aliens_union"})
        self.assertTrue(all("unlabelled" not in x for a in d.values()
                            for s in a["source_site_sets"] for x in s))
        self.assertEqual(d["explicit_aliens_union"]["n_source_site_taxon_incidence"],
                         sum(d[x]["n_source_site_taxon_incidence"] for x in ("archaeophyte","neophyte")))
        self.assertEqual(d["explicit_aliens_union"]["n_actual_original_site_stations_with_class"],35)

    def test_drop_double_empty_source_jaccard_as_missing_not_zero(self):
        sites,status,pairs=artificial()
        x=[set() for _ in sites];x[0]={"arch_0"};x[2]={"arch_0"}
        p=source_class_pairs(pairs,x)
        self.assertEqual(len(p),len(sites)-1+(len(sites)-2))
        self.assertFalse(any(not (x[e["i"]] or x[e["j"]]) for e in p))
        raw=source_subgroup_gap(p,x)
        self.assertTrue(raw["identified"])
        self.assertAlmostEqual(raw["mean_jaccard_same"],0) # same section pairs share none

    def test_source_matched_null_preserves_site_class_richness_and_original_species_frequency(self):
        sites,status,pairs=artificial()
        d=group_sets(sites,status)
        x=d["archaeophyte"]["source_site_sets"]
        original=conditional_source_gap(pairs,sets=x)
        self.assertTrue(original["identified"])
        n=conditional_group_null(sites,pairs,x,nchains=1,draws=5,burnin=100,thin=100)
        self.assertTrue(n["identified"])
        self.assertEqual(n["n_correlated_MCMC_draws"],5)
        self.assertTrue(n["station_subgroup_richness_and_taxon_total_frequency_preserved"])
        tighter=conditional_group_null(sites,pairs,x,nchains=1,draws=5,burnin=100,thin=100,within_stretch=True)
        self.assertTrue(tighter["within_original_source_stretch_taxon_frequency_also_fixed"])
        self.assertEqual(tighter["n_correlated_MCMC_draws"],5)

    def test_actual_original_xlsx_source_hash_is_mandatory(self):
        with tempfile.TemporaryDirectory() as td:
            fake=Path(td)/"original.xlsx";fake.write_bytes(b"not genuine")
            with self.assertRaisesRegex(ValueError,"hash mismatch"):
                read_source_status(fake)

if __name__=="__main__":unittest.main()
