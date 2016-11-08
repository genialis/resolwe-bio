# pylint: disable=missing-docstring
from os.path import join

from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles


class MicroarrayProcessorTestCase(BioProcessTestCase):

    @skipUnlessLargeFiles('affy_test_1.CEL', 'affy_test_2.CEL')
    def test_microarray(self):
        inputs = {'cel': join('large', 'affy_test_1.CEL')}
        cel_1 = self.run_process('upload-microarray-affy', inputs)
        inputs = {'cel': join('large', 'affy_test_2.CEL')}
        cel_2 = self.run_process('upload-microarray-affy', inputs)

        inputs = {'cel': [cel_1.pk, cel_2.pk],
                  'library': 'affy',
                  'logtransform': True}

        self.run_process('microarray-affy-qc', inputs)
