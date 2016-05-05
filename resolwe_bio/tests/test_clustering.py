# pylint: disable=missing-docstring
from .utils import skipDockerFailure, BioProcessTestCase


class ClusteringProcessorTestCase(BioProcessTestCase):

    @skipDockerFailure("Errors with: ERROR: basic:json value in exp_json not "
        "ObjectId but {u'genes': {u'DPU_G0067108': 0.0, ...}} at "
        "etc = self.run_processor('etc:bcm-1-0-0', inputs)")
    def test_hc_clustering(self):
        """Cannot use assertJSON - JSON output contains ETC object IDs."""
        expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM")
        expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")
        expression_3 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")

        inputs = {'expressions': [expression_1.pk, expression_2.pk, expression_3.pk]}
        etc = self.run_processor('etc:bcm-1-0-0', inputs)

        inputs = {
            'etcs': [etc.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067098', 'DPU_G0067100']}
        self.run_processor('clustering:hierarchical:bcm-1-0-0', inputs)
