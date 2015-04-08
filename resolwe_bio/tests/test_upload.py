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

    def test_upload_expression(self):
        inputs = {"exp_type": "TPM"}
        expressions_1 = self.run_processor("import:upload:expression", inputs, Data.STATUS_ERROR)

        inputs = {"exp": "00Hr_tpm.tab.gz", "rc": "00Hr_rc.tab.gz"}
        expressions_2 = self.run_processor("import:upload:expression", inputs, Data.STATUS_ERROR)

        inputs = {"rc": "00Hr_rc.tab.gz"}
        expressions_3 = self.run_processor("import:upload:expression", inputs)
        self.assertFiles(expressions_3, "rc", "00Hr_rc.tab.gz")
        self.assertFiles(expressions_3, 'exp', '00Hr_tpm.tab.gz')

        inputs = {"exp": "00Hr_tpm.tab.gz", "exp_type": "TPM"}
        expressions_4 = self.run_processor("import:upload:expression", inputs)
        self.assertFiles(expressions_4, "rc", "00Hr_rc_2.tab.gz")
        self.assertFiles(expressions_4, 'exp', '00Hr_tpm_2.tab.gz')

        inputs = {"rc": "00Hr_rc.tab.gz", "exp": "00Hr_tpm.tab.gz", "exp_type": "TPM"}
        expressions_5 = self.run_processor("import:upload:expression", inputs)

        self.assertFields(expressions_5, 'exp_type', 'TPM')
        self.assertFiles(expressions_5, 'exp', 'expression_tpm.tab.gz')
        self.assertFiles(expressions_5, 'rc', 'expression_rc.tab.gz')
        self.assertJSON(expressions_5, expressions_5.output['exp_json'], '', 'expression.json.gz')
