# pylint: disable=missing-docstring
from resolwe.test import tag_process
from resolwe_bio.utils.test import with_resolwe_host, KBBioProcessTestCase


class EnrichmentProcessorTestCase(KBBioProcessTestCase):

    @with_resolwe_host
    @tag_process('goenrichment')
    def test_go_enrichment_dicty(self):
        with self.preparation_stage():
            inputs = {'src': 'ontology_dicty_cropped.obo.gz'}
            ontology = self.run_process('upload-obo', inputs)

            inputs = {'src': 'gaf_dicty_cropped.gz', 'source': 'DICTYBASE', 'species': 'Dictyostelium discoideum'}
            annotation = self.run_process('upload-gaf', inputs)

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'pval_threshold': 1,
            'source': 'DICTYBASE',
            'species': 'Dictyostelium discoideum',
            'genes': ['DDB_G0277589', 'DDB_G0286855', 'DDB_G0267640']
        }

        enrichment = self.run_process('goenrichment', inputs)
        self.assertEqual(len(enrichment.process_warning), 0)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms_dicty.json.gz')

    @with_resolwe_host
    @tag_process('goenrichment')
    def test_go_enrichment(self):
        with self.preparation_stage():
            inputs = {'src': 'ontology_mus_cropped.obo.gz'}
            ontology = self.run_process('upload-obo', inputs)

            inputs = {'src': 'gaf_mgi_cropped.gz', 'source': 'MGI', 'species': 'Mus musculus'}
            gaf = self.run_process('upload-gaf', inputs)

        inputs = {
            'ontology': ontology.pk,
            'gaf': gaf.pk,
            'pval_threshold': 1,
            'genes': ['193202', '56535'],
            'source': 'NCBI',
            'species': 'Mus musculus'
        }

        enrichment = self.run_process('goenrichment', inputs)
        self.assertEqual(len(enrichment.process_warning), 1)
        self.assertEqual(enrichment.process_warning[0], "Not all features could be mapped.")
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms_mouse.json.gz')
