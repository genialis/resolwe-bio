# pylint: disable=missing-docstring
from .base import BaseProcessorTestCase
from .utils import PreparedData


class EnrichmentProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_go_enrichment(self):
        inputs = {'src': 'ontology.obo.gz'}
        ontology = self.run_processor('import:upload:ontology', inputs)

        inputs = {'src': 'gene_association.dictyBase.gz'}
        annotation = self.run_processor('import:upload:gaf', inputs)

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'genes': ['DDB_G0272813', 'DDB_G0288677', 'DDB_G0285417',
                      'DDB_G0267442', 'DDB_G0268480', 'DDB_G0283279']}
        enrichment = self.run_processor('goenrichment:bcm-1-0-0', inputs)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms.json')

    def test_go_enrichment_v2(self):
        inputs = {'src': 'ontology.obo.gz'}
        ontology = self.run_processor('import:upload:ontology', inputs)

        inputs = {'src': 'gene_association.dictyBase.gz'}
        annotation = self.run_processor('import:upload:gaf', inputs)

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'genes': ['DDB_G0272813', 'DDB_G0288677', 'DDB_G0285417',
                      'DDB_G0267442', 'DDB_G0268480', 'DDB_G0283279']}
        enrichment = self.run_processor('goenrichment:bcm-2-0-0', inputs)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms_2.json')

    def test_go_enrichment_mouse(self):
        inputs = {'src': 'ontology.obo.gz'}
        ontology = self.run_processor('import:upload:ontology', inputs)

        inputs = {'src': 'gene_association.mgi.gz'}
        annotation = self.run_processor('import:upload:gaf', inputs)

        inputs = {
            'ontology': ontology.pk,
            'gaf': annotation.pk,
            'genes': ['MGI:101765', 'MGI:102956', 'MGI:1925584',
                      'MGI:1933126', 'MGI:109606', 'MGI:109183']}
        enrichment = self.run_processor('goenrichment:bcm-2-0-0', inputs)
        self.assertJSON(enrichment, enrichment.output['terms'], '', 'go_enriched_terms_mouse.json')
