# pylint: disable=missing-docstring
from .base import BaseProcessorTestCase
from .utils import PreparedData


class UploadProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_bam_upload(self):
        inputs = {
            'reads': 'pool24.read1.small.qseq.bz2',
            'reads2': 'pool24.read3.small.qseq.bz2',
            'barcodes': 'pool24.read2.small.qseq.bz2',
            'annotation': 'pool24.tsv'
        }
        obj = self.run_processor('import:upload:multiplexed-paired-end', inputs)
        self.assertFields(obj, 'matched', '437 reads (62.34 %)')
        self.assertFields(obj, 'notmatched', '132 reads (18.83 %)')
        self.assertFields(obj, 'badquality', '131 reads (18.69 %)')
        self.assertFields(obj, 'skipped', '0 reads (0.00 %)')
