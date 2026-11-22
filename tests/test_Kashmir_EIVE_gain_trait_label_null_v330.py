"""v330 synthetic tests: one consistent species trait shuffle across source sites."""
import sys,unittest,random
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import test_Kashmir_EIVE_gain_trait_label_null_v330 as m

NAMES=["genus alpha","genus beta","genus gamma","other delta"]
EIVE={s:{"position":float(i+3),"width":float(i+1),"systems":3}
      for i,s in enumerate(NAMES)}
def fake():
    plots={}
    for i in range(1,17):
        site=f"Station {i}"
        for year in (2014,2017):
            plots[site,year,"Track"]={0,1} if year==2014 else {1,2}
            plots[site,year,"Interior Plot 2"]={0,1}
    for i in range(1,16):
        site=f"Between {i}"
        for year in (2014,2017):
            plots[site,year,"Track"]={0,1} if year==2014 else {1,2}
            plots[site,year,"Interior Plot 2"]={0,1}
    return plots

class TaxonNull(unittest.TestCase):
    def test_real_field_site_count_not_species_record_count(self):
        q=m.site_cohorts(fake(),NAMES,EIVE,"Interior Plot 2",1)
        self.assertEqual(len(q),31)
        self.assertEqual(sum(x["station_or_between"]=="station" for x in q),16)
    def test_stricter_gate_wont_miscode_low_richness_as_zero(self):
        self.assertEqual(m.site_cohorts(fake(),NAMES,EIVE,"Interior Plot 2",3),[])
    def test_source_taxon_pool(self):
        pool=m.validated_source_taxa(fake(),NAMES,EIVE)
        self.assertEqual(pool,sorted(NAMES[:3]))
    def test_observed_gain_identity(self):
        c=m.site_cohorts(fake(),NAMES,EIVE,"Interior Plot 2",1)
        v={s:EIVE[s]["position"] for s in NAMES}
        stats=m.turnover_stats(c,v)
        self.assertAlmostEqual(
             stats["reference_total_minus_track_total"],
             stats["reference_gain_minus_track_gain"]+
             stats["reference_loss_minus_track_loss"])
    def test_permutation_preserves_per_site_taxon_and_fixed_distribution(self):
        pool=sorted(NAMES[:3]);v={s:EIVE[s]["position"] for s in NAMES}
        gr={"genus":pool}
        x=m.permuted_values(pool,v,gr,random.Random(2),"global")
        self.assertEqual(sorted(x[s] for s in pool),sorted(v[s] for s in pool))
    def test_within_genus_permutation_preserves_genus_mean(self):
        pool=NAMES;v={s:EIVE[s]["position"] for s in NAMES}
        gr={"genus":NAMES[:3],"other":[NAMES[3]]}
        x=m.permuted_values(pool,v,gr,random.Random(2),"within_genus")
        self.assertEqual(sum(x[s] for s in NAMES[:3]),
                         sum(v[s] for s in NAMES[:3]))
        self.assertEqual(x[NAMES[3]],v[NAMES[3]])
    def test_exact_no_habitat_randomization_because_only_traits_shuffle(self):
        q=m.site_cohorts(fake(),NAMES,EIVE,"Interior Plot 2",1)
        base={s:EIVE[s]["position"] for s in NAMES}
        got=m.null(q,NAMES[:3],base,"within_genus",draws=19,seed=4)
        self.assertEqual(got["n_source_sites"],31)
        self.assertEqual(got["null_statistics"]["reference_gain_minus_track_gain"]["n_draws"],19)
        self.assertTrue(0<got["null_statistics"]["reference_gain_minus_track_gain"]["technical_two_sided_fraction_plus_one"]<=1)
if __name__=="__main__":unittest.main()
