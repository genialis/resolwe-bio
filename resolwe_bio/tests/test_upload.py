# pylint: disable=missing-docstring
from server.models import Data

from .base import BaseProcessorTestCase
from .utils import PreparedData


class UploadProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_bam_upload(self):
        inputs = {"src": "name_sorted.bam"}
        upload_bam = self.run_processor("import:upload:mapping-bam", inputs)
        self.assertFiles(upload_bam, 'bai', 'bam_upload_index.bai')

        inputs = {"src": "invalid_unsorted.bam"}
        upload_bam = self.run_processor("import:upload:mapping-bam", inputs, Data.STATUS_ERROR)
