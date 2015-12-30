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
            # 'ann': '/studies/GSE22148/concepts/Subjects/Annotations/Batch/B;/studies/GSE22148/concepts/Subjects/Annotations/Batch/A;/studies/GSE22148/concepts/Subjects/Annotations/Batch/C;/studies/GSE22148/concepts/Subjects/Annotations/Age',
            # Alternative annotations
            # 'ann': '/studies/GSE22148/concepts/Subjects/Annotations/Batch/B',
            'ann': '/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Non%20Allergic%20Rhinitis%20Active/0;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Non%20Allergic%20Rhinitis%20Active/NA;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/109,43;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/125,9;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/126,6;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/134,34;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/136,34;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/139,14;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/139,47;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/140,46;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/154,56;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/208,54;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/249,23;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/261,57;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/275,74;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/291,78;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/43,72;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/48,2;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/50,57;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/612,04;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/86,93;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/94,79;/studies/UBIOPRED/concepts/Subject%20History/Medical%20or%20Surgical%20History/Gerd%20Diagnosed/NA',
            'token': 'a8b7c9e2-9b7b-400c-ba4c-2c51b7989f6b'
        }
        annotation = self.run_processor('import:web:transmart:expressions', inputs)
        self.assertFields(annotation, 'expset_type', 'Log2')
        self.assertFileExists(annotation, 'expset')
