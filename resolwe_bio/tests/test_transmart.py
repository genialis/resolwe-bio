# pylint: disable=missing-docstring
import unittest

from .utils import ProcessTestCase


class TranSMARTProcessorTestCase(ProcessTestCase):

    @unittest.skip("test data not ready")
    def test_import(self):
        self.assertDataCount(0)
        inputs = {'exps': 'transmart_log.xlsx'}
        annotation = self.run_processor('import:web:transmart:expressions', inputs)
        self.assertDataCount(7)
        self.assertFields(annotation, 'expset_type', 'Log2')
        self.assertFileExists(annotation, 'expset')

    @unittest.skip("test data not ready")
    def test_import_with_annotation(self):
        self.assertDataCount(0)
        inputs = {'exps': 'transmart_log.xlsx', 'ann': 'transmart_log_annotation.xlsx'}
        annotation = self.run_processor('import:web:transmart:expressions', inputs)
        self.assertDataCount(7)
        self.assertFields(annotation, 'expset_type', 'Log2')
        self.assertFileExists(annotation, 'expset')
