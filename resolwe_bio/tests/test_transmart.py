# pylint: disable=missing-docstring
import unittest
from .utils import ProcessTestCase


class TranSMARTProcessorTestCase(ProcessTestCase):
    @unittest.skip("processor connects to remote tranSMART server")
    def test_import(self):
        self.assertDataCount(0)
        inputs = {'exps': 'transmart_log_exp.xlsx'}
        annotation = self.run_processor('import:web:transmart:expressions', inputs)
        self.assertDataCount(5)
        self.assertFields(annotation, 'expset_type', 'Log2')
        self.assertFileExists(annotation, 'expset')

    @unittest.skip("processor connects to remote tranSMART server")
    def test_import_with_annotation(self):
        self.keep_all()
        inputs = {
            'exps': '/studies/GSE22148/concepts/Expression/Affymetrix%20Human%20Genome%20U133%20Plus%202.0%20Array/Small',
            # Larga data
            # 'exps': '/studies/GSE22148/concepts/Expression/Affymetrix%20Human%20Genome%20U133%20Plus%202.0%20Array/Lung',
            'ann': '/studies/GSE22148/concepts/Subjects/Annotations/Batch/B;/studies/GSE22148/concepts/Subjects/Annotations/Batch/A;/studies/GSE22148/concepts/Subjects/Annotations/Batch/C;/studies/GSE22148/concepts/Subjects/Annotations/Age',
            # Alternative annotation
            # 'ann': '/studies/GSE22148/concepts/Subjects/Annotations/Batch/B',
            'token': '6ef09b1a-b23c-471e-97c1-75d4d48eb5d6'
        }
        annotation = self.run_processor('import:web:transmart:expressions', inputs)
        self.assertFields(annotation, 'expset_type', 'Log2')
        self.assertFileExists(annotation, 'expset')
