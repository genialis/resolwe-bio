# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class JSONProcessorTestCase(BioProcessTestCase):

    def test_save_json_2000(self):
        self.run_processor("test-json", {"n_genes": 2000})  # 0.3 MB

    def test_save_json_20000(self):
        self.run_processor("test-json", {"n_genes": 20000})  # 3.8 MB

    def test_save_json_200000(self):
        self.run_processor("test-json", {"n_genes": 200000})  # 39 MB

    def test_save_json_2000000(self):
        self.run_processor("test-json", {"n_genes": 2000000})  # 399 MB

    def test_save_json_5000000(self):
        self.run_processor("test-json", {"n_genes": 4500000})  # <1 GB

    def test_save_json_6000000(self):
        self.run_processor("test-json", {"n_genes": 5500000})  # >1 GB
