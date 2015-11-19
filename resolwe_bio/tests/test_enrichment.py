# pylint: disable=missing-docstring
from .utils import ProcessTestCase
import unittest


class EnrichmentProcessorTestCase(ProcessTestCase):
    @unittest.skip("test data not ready")
    def test_go_enrichment_v2(self):
        inputs = {'src': 'ontology.obo.gz'}
        ontology = self.run_processor('import:upload:ontology', inputs)

        inputs = {'src': 'gene_association.dictyBase.cropped.gz'}
        annotation = self.run_processor('import:upload:gaf', inputs)

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'genes': ['DDB_G0272813', 'DDB_G0288677', 'DDB_G0285417',
                      'DDB_G0267442', 'DDB_G0268480', 'DDB_G0283279']}

        enrichment = self.run_processor('goenrichment:bcm-2-0-0', inputs)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms.json.gz')

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'genes': ['DDB_G0267640', 'DDB_G0279331', 'DDB_G0289651', 'DDB_G0281087']
        }

        enrichment = self.run_processor('goenrichment:bcm-2-0-0', inputs)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms_2.json.gz')

        inputs = {'src': 'purpureum_ortholog-10-28-2014.cropped.txt.gz'}
        orthologues = self.run_processor('import:upload:orthologues', inputs)

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'orthologues': orthologues.pk,
            'genes': ['DDB_G0272813', 'DDB_G0288677', 'DDB_G0285417',
                      'DDB_G0267442', 'DPU_G0053558', 'DPU_G0071398']}

        enrichment = self.run_processor('goenrichment:bcm-2-0-0', inputs)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms.json.gz')

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'orthologues': orthologues.pk,
            'genes': ['DDB_G0267640', 'DDB_G0279331', 'DDB_G0289651', 'DDB_G0281087']
        }

        enrichment = self.run_processor('goenrichment:bcm-2-0-0', inputs)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms_2.json.gz')

    @unittest.skip("test data not ready")
    def test_go_enrichment_mouse(self):
        inputs = {'src': 'ontology.obo.gz'}
        ontology = self.run_processor('import:upload:ontology', inputs)

        inputs = {'src': 'gene_association.mgi.cropped.gz'}
        annotation = self.run_processor('import:upload:gaf', inputs)

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'genes': ['MGI:101765', 'MGI:102956', 'MGI:1925584',
                      'MGI:1933126', 'MGI:109606', 'MGI:109183']}
        enrichment = self.run_processor('goenrichment:bcm-2-0-0', inputs)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms_mouse.json.gz')
