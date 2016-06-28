# pylint: disable=missing-docstring
from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase


class EnrichmentProcessorTestCase(BioProcessTestCase):

    def test_go_enrichment_dicty(self):
        inputs = {'src': 'ontology_dicty_cropped.obo.gz'}
        ontology = self.run_processor('upload-obo', inputs)

        inputs = {'src': 'gaf_dicty_cropped.gz'}
        annotation = self.run_processor('upload-gaf', inputs)

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'pval_threshold': 1,
            'genes': ['DDB_G0277589', 'DDB_G0286855', 'DDB_G0267640']}

        enrichment = self.run_processor('goenrichment-bcm', inputs)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms_dicty.json.gz')

        inputs = {'src': 'purpureum_ortholog-10-28-2014.cropped.txt.gz'}
        orthologues = self.run_processor('upload-orthologues', inputs)

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'orthologues': orthologues.pk,
            'pval_threshold': 1,
            'genes': ['DPU_G0074602', 'DDB_G0286855', 'DPU_G0074318']}

        enrichment = self.run_processor('goenrichment-bcm', inputs)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms_dicty.json.gz')

    def test_go_enrichment_mouse(self):
        inputs = {'src': 'ontology_mus_cropped.obo.gz'}
        ontology = self.run_processor('upload-obo', inputs)

        inputs = {'src': 'gaf_mgi_cropped.gz'}
        annotation = self.run_processor('upload-gaf', inputs)

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'pval_threshold': 1,
            'genes': ['MGI:1929646', 'MGI:107486']}
        enrichment = self.run_processor('goenrichment-bcm', inputs)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms_mouse.json.gz')
