# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase


class ImportProcessorTestCase(BioProcessTestCase):
    @tag_process('import-sra')
    def external_test_sra(self):
        # single-end read from Polyak RNA-seq demo dataset
        inputs = {
            'sra_accession': 'SRR1661332',
            'advanced': {
                'max_spot_id': 1,
            },
        }
        self.run_process('import-sra', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        import_sra = Data.objects.last()
        self.assertFiles(import_sra, 'fastq', ['SRR1661332.fastq.gz'], compression='gzip')
        del import_sra.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        self.assertFields(import_sra, 'fastqc_url', [{'file': 'fastqc/SRR1661332_fastqc/fastqc_report.html',
                                                      'refs': ['fastqc/SRR1661332_fastqc']}])

        # paired-end read from Zoghbi RNA-seq demo dataset
        inputs = {
            'sra_accession': 'SRR2124780',
            'advanced': {
                'max_spot_id': 1,
            },
        }
        self.run_process('import-sra', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        import_sra = Data.objects.last()
        self.assertFiles(import_sra, 'fastq', ['SRR2124780.fastq.gz'], compression='gzip')
        del import_sra.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        self.assertFields(import_sra, 'fastqc_url', [{'file': 'fastqc/SRR2124780_1_fastqc/fastqc_report.html',
                                                      'refs': ['fastqc/SRR2124780_1_fastqc']}])
        del import_sra.output['fastqc_url2'][0]['total_size']  # Non-deterministic output.
        self.assertFields(import_sra, 'fastqc_url2', [{'file': 'fastqc/SRR2124780_2_fastqc/fastqc_report.html',
                                                       'refs': ['fastqc/SRR2124780_2_fastqc']}])
