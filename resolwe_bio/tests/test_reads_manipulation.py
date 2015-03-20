from .base import BaseProcessorTestCase
from .utils import PreparedData


class ReadsProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_merge_reads(self):
        reads = self.prepare_reads()
        reads2 = self.prepare_reads()

        inputs = {
            'reads_1': reads.pk,
            'reads_2': reads2.pk}
        merged_reads = self.run_processor('reads:merge', inputs)
        self.assertFiles(merged_reads, 'fastq', 'paired_end_forward.fastq.gz')
