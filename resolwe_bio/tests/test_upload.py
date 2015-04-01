# pylint: disable=missing-docstring
import mongoengine

from server.models import Data

from .base import BaseProcessorTestCase
from .utils import PreparedData


class UploadProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_bam_upload(self):
        inputs = {"src": "alignment_name_sorted.bam"}
        upload_bam = self.run_processor("import:upload:mapping-bam", inputs)
        self.assertFiles(upload_bam, 'bam', 'alignment_position_sorted.bam')
        self.assertFiles(upload_bam, 'bai', 'alignment_bam_upload_index.bai')

        inputs = {"src": "alignment_invalid_unsorted.bam"}
        upload_bam = self.run_processor("import:upload:mapping-bam", inputs, Data.STATUS_ERROR)

        inputs = {"src": "alignment_position_sorted.bam", "src2": "alignment_bam_upload_index.bai"}
        self.assertRaises(mongoengine.ValidationError, self.run_processor, "import:upload:mapping-bam-indexed", inputs)

        inputs = {"src": "alignment_position_sorted.bam", "src2": "alignment_bam_upload_index.bam.bai"}
        upload_bam = self.run_processor("import:upload:mapping-bam-indexed", inputs, Data.STATUS_ERROR)
        self.assertFields(upload_bam, 'proc.error', 'BAI should have the same name as BAM with .bai extension')

        inputs = {"src": "alignment_position_sorted.bam", "src2": "alignment_position_sorted.bam.bai"}
        upload_bam = self.run_processor("import:upload:mapping-bam-indexed", inputs)
        self.assertFiles(upload_bam, 'bam', 'alignment_position_sorted.bam')
        self.assertFiles(upload_bam, 'bai', 'alignment_position_sorted.bam.bai')
