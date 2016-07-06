# pylint: disable=missing-docstring
import unittest
from resolwe_bio.utils.test import BioProcessTestCase


class MicroarrayProcessorTestCase(BioProcessTestCase):
    @unittest.skip("test data not ready")
    def test_microarray(self):
        cel_1 = self.run_processor("upload-microarray-affy", {"cel": "affy_v3_example.CEL"})
        # self.assertFile(cel_1, 'cel', 'sample1.CEL')

        cel_2 = self.run_processor("upload-microarray-affy", {"cel": "affy_v3_example_2.CEL"})
        # self.assertFile(cel_2, 'cel', 'sample2.CEL')

        # self.run_processor('microarray-affy-qc', inputs)

        inputs = {'cel': [cel_1.pk, cel_2.pk],
                  'library': 'affy',
                  'logtransform': True}

        self.run_processor('microarray-affy-qc', inputs)
