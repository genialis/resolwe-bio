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
            'ann': '/studies/GSE22148/concepts/Subjects/Annotations/PlateID',
            # Alternative annotation
            # 'ann': '/studies/GSE22148/concepts/Subjects/Annotations/Batch/B',
            'token': '15f3a910-9c8f-4e12-a2d7-b5d9bef67e75'
        }
        annotation = self.run_processor('import:web:transmart:expressions', inputs)
        self.assertFields(annotation, 'expset_type', 'Log2')
        self.assertFileExists(annotation, 'expset')
