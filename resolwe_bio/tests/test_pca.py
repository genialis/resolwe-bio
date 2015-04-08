# pylint: disable=missing-docstring
from .base import BaseProcessorTestCase
from .utils import PreparedData


class PcaProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_pca(self):
        genome = self.prepare_genome()
        expression_1 = self.prepare_expression(f_rc='00Hr_rc.tab.gz', f_exp='00Hr_tpm.tab.gz', f_type="TPM")
        expression_2 = self.prepare_expression(f_rc='20Hr_rc.tab.gz', f_exp='20Hr_tpm.tab.gz', f_type="TPM")

        inputs = {'expressions': [expression_1.pk, expression_2.pk]}
        pca = self.run_processor('pca:1-0-0', inputs)
        self.assertJSON(pca, pca.output['pca'], 'flot.data', 'pca.json.gz')
