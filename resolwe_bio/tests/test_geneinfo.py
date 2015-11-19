# pylint: disable=missing-docstring
from .utils import ProcessTestCase


class GIProcessorTestCase(ProcessTestCase):
    def test_gi(self):
        inputs = {'src': 'mouse_gene_info.txt'}
        self.run_processor("import:upload:geneinfo", inputs)
