# pylint: disable=missing-docstring
from .base import BaseProcessorTestCase
from .utils import PreparedData


class MicroarrayProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_microarray(self):
        cel_1 = self.run_processor("import:upload:microarray:affy", {"cel": "sample1.CEL"})
        self.assertFiles(cel_1, 'cel', 'sample1.CEL')

        cel_2 = self.run_processor("import:upload:microarray:affy", {"cel": "sample2.CEL"})
        self.assertFiles(cel_2, 'cel', 'sample2.CEL')

        cel_3 = self.run_processor("import:upload:microarray:affy", {"cel": "sample3.CEL"})
        self.assertFiles(cel_3, 'cel', 'sample3.CEL')

        inputs = {'cel': [cel_1.pk, cel_2.pk, cel_3.pk],
                  'library': 'affy',
                  'logtransform': True}

        self.run_processor('microarray:affy:qc', inputs)
