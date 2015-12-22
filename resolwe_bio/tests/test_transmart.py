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

    #@unittest.skip("processor connects to remote tranSMART server")
    def test_import_with_annotation(self):
        self.keep_all()
        inputs = {
            'exps': '/studies/GSE22148/concepts/Expression/Affymetrix%20Human%20Genome%20U133%20Plus%202.0%20Array/Lung', 
            'ann': '/studies/GSE22148/concepts/Subjects/Annotations/PlateID',
            'token': '0779d50e-822f-403e-bdf8-1f701199f3f6'
        }
        annotation = self.run_processor('import:web:transmart:expressions', inputs)
        self.assertFields(annotation, 'expset_type', 'Log2')
        self.assertFileExists(annotation, 'expset')
