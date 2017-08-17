# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe_bio.utils.test import BioProcessTestCase


class ImportProcessorTestCase(BioProcessTestCase):
    def external_test_sra(self):
        # single-end reads from Polyak RNA-seq demo dataset
        inputs = {'sra_accession': 'SRR1661332', 'prefetch': False, 'max_spot_id': 1}
        self.run_process('import-sra', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        import_sra = Data.objects.last()
        self.assertFiles(import_sra, 'fastq', ['SRR1661332.fastq.gz'], compression='gzip')
        self.assertFields(import_sra, 'fastqc_url', [{'file': 'fastqc/SRR1661332_fastqc/fastqc_report.html',
                                                      'refs': ['fastqc/SRR1661332_fastqc']}])

        # paired-end reads from Zoghbi RNA-seq demo dataset
        inputs = {'sra_accession': 'SRR2124780', 'prefetch': False, 'max_spot_id': 1}
        self.run_process('import-sra', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        import_sra = Data.objects.last()
        self.assertFiles(import_sra, 'fastq', ['SRR2124780.fastq.gz'], compression='gzip')
        self.assertFields(import_sra, 'fastqc_url', [{'file': 'fastqc/SRR2124780_1_fastqc/fastqc_report.html',
                                                      'refs': ['fastqc/SRR2124780_1_fastqc']}])
        self.assertFields(import_sra, 'fastqc_url2', [{'file': 'fastqc/SRR2124780_2_fastqc/fastqc_report.html',
                                                       'refs': ['fastqc/SRR2124780_2_fastqc']}])
