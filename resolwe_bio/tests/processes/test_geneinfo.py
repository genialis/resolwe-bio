# pylint: disable=missing-docstring
from .utils import skipDockerFailure, BioProcessTestCase


class GIProcessorTestCase(BioProcessTestCase):

    @skipDockerFailure("Fails with: ImportError: No module named xlrd")
    def test_gi(self):
        inputs = {'src': 'mouse_gene_info.txt'}
        self.run_processor("upload-geneinfo", inputs)
