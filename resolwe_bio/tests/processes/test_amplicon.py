# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class AmpliconProcessorTestCase(BioProcessTestCase):

    def test_bwa_trim(self):
        inputs = {
            'src1': ['56GSID_10k_mate1.fastq.gz'],
            'src2': ['56GSID_10k_mate2.fastq.gz']}
        reads = self.run_process('upload-fastq-paired', inputs)

        genome = self.run_process('upload-genome', {'src': 'hs_b37_chr2_small.fasta.gz'})
        master_file = self.run_process('upload-master-file', {'src': '56G_masterfile_test.txt'})

        inputs = {
            'master_file': master_file.id,
            'genome': genome.id,
            'reads': reads.id
        }
        bwa_trim = self.run_process('align-bwa-trim', inputs)

        self.assertFile(bwa_trim, 'stats', 'bwa_trim_stats.txt')
