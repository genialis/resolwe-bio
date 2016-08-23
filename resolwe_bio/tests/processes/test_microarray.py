# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles


class MicroarrayProcessorTestCase(BioProcessTestCase):

    @skipUnlessLargeFiles('large/affy_test_1.CEL', 'large/affy_test_2.CEL')
    def test_microarray(self):
        cel_1 = self.run_process('upload-microarray-affy', {'cel': 'large/affy_test_1.CEL'})
        cel_2 = self.run_process('upload-microarray-affy', {'cel': 'large/affy_test_2.CEL'})

        inputs = {'cel': [cel_1.pk, cel_2.pk],
                  'library': 'affy',
                  'logtransform': True}

        self.run_process('microarray-affy-qc', inputs)
