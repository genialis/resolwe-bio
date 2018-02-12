# pylint: disable=missing-docstring
from resolwe.flow.models import Secret
from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase


class FilesToReadsTestCase(BioProcessTestCase):
    @tag_process('basespace-file-import', 'files-to-fastq-single')
    def external_test_files_fq_single(self):
        with self.preparation_stage():
            # Token with limited scope pre-obtained from dedicated BaseSpace testing app.
            handle = Secret.objects.create_secret('9bdf059c759a429f8af52ca084130060', self.admin)

            import_inputs_1 = {'file_id': '9461130722', 'access_token_secret': {'handle': handle}}
            basespace_import_1 = self.run_process('basespace-file-import', import_inputs_1)

            import_inputs_2 = {'file_id': '9461121664', 'access_token_secret': {'handle': handle}}
            basespace_import_2 = self.run_process('basespace-file-import', import_inputs_2)

        reads = self.run_process('files-to-fastq-single', {'src': [basespace_import_1.pk, basespace_import_2.pk]})

        self.assertFiles(reads, 'fastq', ['Test_S1_L001_R1_001.fastq.gz',
                                          'Test_S1_L002_R1_001.fastq.gz'], compression='gzip')
        del reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        del reads.output['fastqc_url'][1]['total_size']  # Non-deterministic output.
        self.assertFields(reads, 'fastqc_url', [{'file': 'fastqc/Test_S1_L001_R1_001_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/Test_S1_L001_R1_001_fastqc']},
                                                {'file': 'fastqc/Test_S1_L002_R1_001_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/Test_S1_L002_R1_001_fastqc']}])
