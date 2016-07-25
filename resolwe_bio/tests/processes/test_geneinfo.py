# pylint: disable=missing-docstring
from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase


class GIProcessorTestCase(BioProcessTestCase):

    def test_gi(self):
        inputs = {'src': 'mouse_gene_info.txt'}
        self.run_processor("upload-geneinfo", inputs)
