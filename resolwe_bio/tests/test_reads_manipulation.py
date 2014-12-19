from .base import BaseProcessorTestCase


class ReadsProcessorTestCase(BaseProcessorTestCase):
    def prepair_reads(self):
        inputs = {'src': 'reads.fastq.gz'}
        reads = self.run_processor('import:upload:reads-fastq', inputs)
        self.assertDone(reads)
        self.assertFields(reads, 'bases', 35)
        return reads

    def test_merge_reads(self):
        reads = self.prepair_reads()
        reads2 = self.prepair_reads()

        inputs = {'reads_1': reads.pk, 'reads_2': reads2.pk}
        merged_reads = self.run_processor('reads:merge', inputs)
        self.assertDone(merged_reads)
        self.assertFiles(merged_reads, 'fastq', 'paired_end_forward.fastq.gz')
