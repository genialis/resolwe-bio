# pylint: disable=missing-docstring
from .utils import ProcessTestCase


class PcaProcessorTestCase(ProcessTestCase):
    def test_pca(self):
        expression_1 = self.prepare_expression(f_rc='00Hr_rc.tab.gz', f_exp='00Hr_tpm.tab.gz', f_type="TPM")
        expression_2 = self.prepare_expression(f_rc='20Hr_rc.tab.gz', f_exp='20Hr_tpm.tab.gz', f_type="TPM")

        inputs = {'exps': [expression_1.pk, expression_2.pk]}
        pca = self.run_processor('pca:1-0-0', inputs)
        self.assertJSON(pca, pca.output['pca'], 'flot.data', 'pca_plot.json.gz')
        self.assertJSON(pca, pca.output['pca'], 'explained_variance_ratios', 'pca_ratios.json.gz')
        self.assertJSON(pca, pca.output['pca'], 'components', 'pca_components.json.gz')

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067102', 'DPU_G0067098']
        }

        pca = self.run_processor('pca:1-0-0', inputs)
        self.assertJSON(pca, pca.output['pca'], 'flot.data', 'pca_plot_w_genes.json.gz')

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': ['DPU_G0067098', 'DPU_G0067100', 'DPU_G0067104']  # all zero
        }
        pca = self.run_processor('pca:1-0-0', inputs)

        self.assertJSON(pca, pca.output['pca'], 'flot.data', 'pca_filtered_zeros.json.gz')
        self.assertTrue(pca.output['proc']['warning'])
