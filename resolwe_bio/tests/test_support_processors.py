from .base import BaseProcessorTestCase
from server.models import Data


class CompatibilityProcessorTestCase(BaseProcessorTestCase):
    def prepair_bam(self):
        inputs = {'src': 'sp_test.bam'}
        reads = self.run_processor('import:upload:mapping-bam', inputs, Data.STATUS_DONE)
        self.assertFiles(reads, 'bam', 'sp_test.bam')
        self.assertFiles(reads, 'bai', 'sp_test.bam.bai')
        return reads

    def prepair_genome(self):
        inputs = {'src': 'sp_test.fasta'}
        genome = self.run_processor('import:upload:genome-fasta', inputs, Data.STATUS_DONE)
        self.assertFiles(genome, 'fasta', 'sp_test_genome.fasta.gz', gzipped=True)
        return genome

    def prepair_annotation(self):
        inputs = {'src': 'sp_test.gtf'}
        annotation = self.run_processor('import:upload:annotation-gtf', inputs, Data.STATUS_DONE)
        self.assertFiles(annotation, 'gtf', 'sp_test.gtf')
        return annotation

    def test_reference_compatibility(self):
        mapping = self.prepair_bam()
        genome = self.prepair_genome()
        annotation = self.prepair_annotation()

        inputs = {'reference': genome.pk, 'bam': mapping.pk, 'annot': annotation.pk}
        compatibility_test = self.run_processor('reference_compatibility', inputs, Data.STATUS_DONE)
        self.assertFiles(compatibility_test, 'report_file', 'sp_test_compatibility_report.txt')
