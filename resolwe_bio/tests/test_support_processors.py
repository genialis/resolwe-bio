from .base import BaseProcessorTestCase
from .utils import PreparedData
from server.models import Data


class CompatibilityProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_reference_compatibility(self):
        mapping = self.prepare_bam()
        genome = self.prepare_genome('sp_test.fasta')
        annotation = self.prepare_annotation()

        inputs = {'reference': genome.pk, 'bam': mapping.pk, 'annot': annotation.pk}
        compatibility_test = self.run_processor('reference_compatibility', inputs, Data.STATUS_DONE)
        self.assertFiles(compatibility_test, 'report_file', 'sp_test_compatibility_report.txt')
