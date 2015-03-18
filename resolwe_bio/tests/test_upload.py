from .base import BaseProcessorTestCase
from server.models import Data


class UploadProcessorTestCase(BaseProcessorTestCase):
    def test_bam_upload(self):
        inputs = {"src": "name_sorted.bam"}
        upload_bam = self.run_processor("import:upload:mapping-bam", inputs, Data.STATUS_DONE)
        self.assertFiles(upload_bam, 'bai', 'bam_upload_index.bai')

        inputs = {"src": "invalid_unsorted.bam"}
        upload_bam = self.run_processor("import:upload:mapping-bam", inputs, Data.STATUS_ERROR)
