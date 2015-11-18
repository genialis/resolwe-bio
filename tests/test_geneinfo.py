# pylint: disable=missing-docstring

from .base import BaseProcessorTestCase
from .utils import PreparedData


class GIProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_gi(self):
        inputs = {'src': 'mouse_gene_info.xls'}
        self.run_processor("import:upload:geneinfo", inputs)

        inputs = {'src': 'mouse_gene_info.txt'}
        self.run_processor("import:upload:geneinfo", inputs)
