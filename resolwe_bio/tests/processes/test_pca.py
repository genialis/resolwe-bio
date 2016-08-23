# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class PcaProcessorTestCase(BioProcessTestCase):

    def test_pca(self):
        expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM")
        expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")

        inputs = {'exps': [expression_1.pk, expression_2.pk]}
        pca = self.run_process('pca', inputs)
        self.assertJSON(pca, pca.output['pca'], 'flot.data', 'pca_plot.json.gz')
        self.assertJSON(pca, pca.output['pca'], 'explained_variance_ratios', 'pca_ratios.json.gz')
        self.assertJSON(pca, pca.output['pca'], 'components', 'pca_components.json.gz')

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067102', 'DPU_G0067098']
        }

        pca = self.run_process('pca', inputs)
        self.assertJSON(pca, pca.output['pca'], 'flot.data', 'pca_plot_w_genes.json.gz')

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': ['DPU_G0067098', 'DPU_G0067100', 'DPU_G0067104']  # all zero
        }
        pca = self.run_process('pca', inputs)

        self.assertJSON(pca, pca.output['pca'], 'flot.data', 'pca_filtered_zeros.json.gz')
        self.assertEqual(pca.process_warning[0], "Filtering removed all PCA attributes.")
