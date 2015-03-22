# pylint: disable=missing-docstring
from server.models import Data

from .base import BaseProcessorTestCase
from .utils import PreparedData


class GIProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_gi(self):
        inputs = {'src': 'mouse_gene_info.xls'}
        gi = self.run_processor("import:upload:geneinfo", inputs)

        inputs = {'src': 'mouse_gene_info.txt'}
        gi = self.run_processor("import:upload:geneinfo", inputs)
