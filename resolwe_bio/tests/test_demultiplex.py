# pylint: disable=missing-docstring
from .utils import ProcessTestCase


class UploadProcessorTestCase(ProcessTestCase):
    def test_bam_upload(self):
        inputs = {
            'reads': 'pool24.read1.small.qseq.bz2',
            'reads2': 'pool24.read3.small.qseq.bz2',
            'barcodes': 'pool24.read2.small.qseq.bz2',
            'annotation': 'pool24.tsv'
        }
        obj = self.run_processor('import:upload:multiplexed-paired-end', inputs)
        self.assertFields(obj, 'matched', '5 reads (7.04 %)')
        self.assertFields(obj, 'notmatched', '51 reads (71.83 %)')
        self.assertFields(obj, 'badquality', '13 reads (18.31 %)')
        self.assertFields(obj, 'skipped', '0 reads (0.00 %)')

        inputs = {
            'reads': 'pool24.read3.small.qseq.bz2',
            'reads2': 'pool24.read1.small.qseq.bz2',
            'barcodes': 'pool24.read2.small.qseq.bz2',
            'annotation': 'pool24.tsv'
        }
        obj = self.run_processor('import:upload:multiplexed-paired-end', inputs)
        self.assertFields(obj, 'matched', '5 reads (7.04 %)')
        self.assertFields(obj, 'notmatched', '51 reads (71.83 %)')
        self.assertFields(obj, 'badquality', '13 reads (18.31 %)')
        self.assertFields(obj, 'skipped', '0 reads (0.00 %)')

        inputs = {
            'reads': 'pool24.read1.small.qseq.bz2',
            'reads2': 'pool24.read3.small.qseq.bz2',
            'barcodes': 'pool24.read2.small.qseq.bz2',
            'annotation': 'pool24.2.tsv'
        }
        obj = self.run_processor('import:upload:multiplexed-paired-end', inputs)
        self.assertFields(obj, 'matched', '5 reads (7.04 %)')
        self.assertFields(obj, 'notmatched', '51 reads (71.83 %)')
        self.assertFields(obj, 'badquality', '13 reads (18.31 %)')
        self.assertFields(obj, 'skipped', '0 reads (0.00 %)')

        inputs = {
            'reads': 's_1_1.qseq.small.txt.bz2',
            'reads2': 's_1_3.qseq.small.txt.bz2',
            'barcodes': 's_1_2.qseq.small.txt.bz2',
            'annotation': 's_1.tsv'
        }
        obj = self.run_processor('import:upload:multiplexed-paired-end', inputs)
        self.assertFields(obj, 'matched', '12 reads (11.88 %)')
        self.assertFields(obj, 'notmatched', '80 reads (79.21 %)')
        self.assertFields(obj, 'badquality', '8 reads (7.92 %)')
        self.assertFields(obj, 'skipped', '0 reads (0.00 %)')
