"""v28.9 synthetic schema/permissions tests; NO synthetic plants used as empirical evidence."""
import io
import sys
import unittest
import zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import audit_Figshare_Sweden_nonrail_infrastructure_source_opportunity_v289 as m

class SchemaOnlyTests(unittest.TestCase):
    def test_study_id_predeclared_dois(self):
        self.assertEqual(set(m.ARTICLES),{13286141,22178651})
        with self.assertRaises(ValueError):
            m.article_audit(-5)

    def test_wide_scientific_binomials_detected(self):
        header=['plot_id','habitat']+[f'Genus species{chr(97+i)}' for i in range(20)]
        res=m.plant_columns(header)
        self.assertTrue(res['plausible_wide_site_by_taxon'])
        self.assertTrue(res['site_identifier_column_present'])

    def test_long_occurrence_format_detected(self):
        a=m.plant_columns(['quadrat','species','cover','year'])
        self.assertTrue(a['plausible_long_site_by_taxon'])
        self.assertTrue(a['plant_occurrence_column_present'])

    def test_richness_only_not_false_taxon_matrix(self):
        a=m.plant_columns(['landscape','habitat','plant_richness','n_native','Latitude'])
        self.assertFalse(a['plausible_wide_site_by_taxon'])
        self.assertFalse(a['plausible_long_site_by_taxon'])

    def test_tabular_schema_never_serializes_source_rows(self):
        source=b'plot,species,cover\nSENSITIVE-LOCATION-123,Plantus examplea,1\n'
        s=m.summarize_table(source,'.csv')
        self.assertEqual(s['n_rows_excluding_first_line'],1)
        self.assertNotIn('SENSITIVE-LOCATION-123',str(s))
        self.assertNotIn('Plantus examplea',str(s))

    def test_xlsx_original_shape_without_rows(self):
        import openpyxl
        book=openpyxl.Workbook()
        ws=book.active;ws.title="Plants"
        ws.append(['site_id','species','presence'])
        ws.append(['PRIVATE-LOCATION','Taxon speciesa',1])
        buf=io.BytesIO();book.save(buf)
        s=m.summarize_table(buf.getvalue(),'.xlsx')
        self.assertEqual(s['n_sheets'],1)
        self.assertEqual(s['sheets'][0]['source_max_row'],2)
        self.assertTrue(s['sheets'][0]['header_structure']['plausible_long_site_by_taxon'])
        self.assertNotIn('PRIVATE-LOCATION',str(s))

    def test_zip_safely_contains_schema_only(self):
        buf=io.BytesIO()
        with zipfile.ZipFile(buf,'w') as z:
            z.writestr('original/sensitive.csv','site_id,species,cover\nPRIVATE99,Plantus examplea,3')
        r=m.inspect_file(buf.getvalue(),'publisher_original.zip')
        self.assertEqual(r['member_count'],1)
        self.assertEqual(r['members'][0]['table']['n_rows_excluding_first_line'],1)
        self.assertNotIn('PRIVATE99',str(r))

if __name__=='__main__':unittest.main()
