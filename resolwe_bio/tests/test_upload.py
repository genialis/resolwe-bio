from .base import BaseProcessorTestCase


class UploadProcessorTestCase(BaseProcessorTestCase):
    def test_bam_upload(self):
        inputs = {"src": "name_sorted.bam"}
        upload_bam = self.run_processor("import:upload:mapping-bam", inputs)

        self.assertDone(upload_bam)
        self.assertFiles(upload_bam, 'bai', 'bam_upload_index.bai')

        inputs = {"src": "invalid_unsorted.bam"}
        upload_bam = self.run_processor("import:upload:mapping-bam", inputs)
        self.assertError(upload_bam)
