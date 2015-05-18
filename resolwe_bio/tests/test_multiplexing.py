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
        self.run_processor('import:upload:multiplexed-paired-end', inputs)
