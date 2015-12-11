# pylint: disable=missing-docstring
from .utils import BioProcessTestCase


class GIProcessorTestCase(BioProcessTestCase):
    def test_gi(self):
        inputs = {'src': 'mouse_gene_info.txt'}
        self.run_processor("import:upload:geneinfo", inputs)
