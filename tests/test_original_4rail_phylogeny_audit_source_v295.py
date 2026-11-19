"""v29.5 synthetic unit tests. No field outcomes are claimed using fixtures."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import prepare_Europe_4original_rail_phylogeny_audit_source_v295 as m

def fake_studies():
    names=[f'GenusA species{chr(97+i)}' for i in range(20)]
    traits={key: {taxon:(float(i)-5.,float(i)*.1)
                  for i,taxon in enumerate(names)} for key in m.STUDIES}
    cohorts={}
    for j,k in enumerate(m.STUDIES):
        cohorts[k]=[{
            'station':f'source_{k}_station_{i}',
            'group':str(i%3),
            'taxa':set(names[(ii+2*i+j)%20] for ii in range(9)),
        } for i in range(m.N[k])]
    return cohorts,traits

class OriginalPhylogenyAudit(unittest.TestCase):
    def test_171_true_source_field_units_not_equal_stations(self):
        cohorts,traits=fake_studies()
        rows,qa=m.encode_original_source_rows(cohorts,traits)
        self.assertEqual(sum(v['true_original_field_sites'] for v in qa.values()),171)
        self.assertEqual(len(rows),171*9)
        self.assertEqual(set(qa),set(m.STUDIES))
        self.assertNotIn('source_plant_taxa',str(qa))
        for row in rows:
            self.assertEqual(set(row),{'study','site_num','group_num',
                       'taxon','zlog_canopy_height','zlog_seed_mass'})

    def test_no_richness_imputation_missing_real_trait(self):
        cohorts,traits=fake_studies()
        missing=next(iter(cohorts[m.STUDIES[0]][0]['taxa']))
        del traits[m.STUDIES[0]][missing]
        with self.assertRaises(ValueError):
            m.encode_original_source_rows(cohorts,traits)

    def test_no_invalid_binomial_source_species(self):
        cohorts,traits=fake_studies()
        cohort=cohorts[m.STUDIES[0]][0]
        cohort['taxa'].add('not_binomial')
        traits[m.STUDIES[0]]['not_binomial']=(1.,2.)
        with self.assertRaises(ValueError):
            m.encode_original_source_rows(cohorts,traits)

    def test_true_original_paper_group_not_split_by_country(self):
        cohorts,traits=fake_studies()
        rows,qa=m.encode_original_source_rows(cohorts,traits)
        self.assertEqual(qa['HR_BA_2017_one_study']['true_original_field_sites'],15)
        self.assertEqual(qa['SI_Slovenia_100m_track']['true_original_field_sites'],89)
        self.assertEqual(len(qa),4)

    def test_source_study_sites_are_not_synthetic_sitepair_replicates(self):
        cohorts,traits=fake_studies()
        rows,qa=m.encode_original_source_rows(cohorts,traits)
        self.assertEqual(len({(x['study'],x['site_num']) for x in rows}),171)
        self.assertEqual(len({(x['study'],x['site_num'],x['taxon']) for x in rows}),len(rows))

    def test_source_population_or_study_drift_fails_closed(self):
        cohorts,traits=fake_studies()
        del cohorts[m.STUDIES[1]][0]
        with self.assertRaises(ValueError):
            m.encode_original_source_rows(cohorts,traits)

if __name__=='__main__':unittest.main()
