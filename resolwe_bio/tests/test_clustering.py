# pylint: disable=missing-docstring
from .utils import ProcessTestCase


class ClusteringProcessorTestCase(ProcessTestCase):
    def test_hc_clustering(self):
        """Cannot use assertJSON - JSON output contains ETC object IDs."""
        expression_1 = self.prepare_expression(f_rc='00Hr_rc.tab.gz', f_exp='00Hr_tpm.tab.gz', f_type="TPM")
        expression_2 = self.prepare_expression(f_rc='20Hr_rc.tab.gz', f_exp='20Hr_tpm.tab.gz', f_type="TPM")
        expression_3 = self.prepare_expression(f_rc='20Hr_rc.tab.gz', f_exp='20Hr_tpm.tab.gz', f_type="TPM")

        inputs = {'expressions': [expression_1.pk, expression_2.pk, expression_3.pk]}
        etc = self.run_processor('etc:bcm-1-0-0', inputs)

        inputs = {
            'etcs': [etc.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067098', 'DPU_G0067100']}
        self.run_processor('clustering:hierarchical:bcm-1-0-0', inputs)
