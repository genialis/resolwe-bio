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
        self.assertDataCount(0)
        inputs = {'exps': 'transmart_log_exp.xlsx', 'ann': 'transmart_annotation.xlsx'}
        annotation = self.run_processor('import:web:transmart:expressions', inputs)
        self.assertDataCount(5)
        self.assertFields(annotation, 'expset_type', 'Log2')
        self.assertFileExists(annotation, 'expset')
