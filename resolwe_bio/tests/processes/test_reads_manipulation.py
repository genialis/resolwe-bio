# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class ReadsProcessorTestCase(BioProcessTestCase):

    def test_merge_reads(self):
        reads = self.prepare_reads()
        reads2 = self.prepare_reads()

        inputs = {
            'reads_1': reads.pk,
            'reads_2': reads2.pk}
        merged_reads = self.run_process('reads-merge', inputs)
        self.assertFiles(merged_reads, 'fastq', ['paired_end_forward.fastq.gz'], compression='gzip')
        self.assertFiles(merged_reads, 'fastq2', ['paired_end_reverse.fastq.gz'], compression='gzip')
        self.assertFields(merged_reads, "fastqc_url", [{'file': 'fastqc/fw_reads_fastqc/fastqc_report.html',
                                                        'refs': ['fastqc/fw_reads_fastqc'],
                                                        'size': 311414}])
        self.assertFields(merged_reads, "fastqc_url2", [{'file': 'fastqc/rw_reads_fastqc/fastqc_report.html',
                                                         'refs': ['fastqc/rw_reads_fastqc'],
                                                         'size': 311414}])
