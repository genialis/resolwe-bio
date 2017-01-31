# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class PcaProcessorTestCase(BioProcessTestCase):

    def test_pca(self):
        expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM")
        expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")

        inputs = {'exps': [expression_1.pk, expression_2.pk]}
        pca = self.run_process('pca', inputs)
        saved_json, test_json = self.get_json('pca_plot.json.gz', pca.output['pca'])
        self.assertEqual(test_json['flot']['data'], saved_json['flot']['data'])
        self.assertEqual(test_json['explained_variance_ratios'], saved_json['explained_variance_ratios'])
        self.assertEqual(test_json['components'], saved_json['components'])

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': ['DPU_G0067098', 'DPU_G0067100', 'DPU_G0067104']  # all zero
        }

        pca = self.run_process('pca', inputs)
        self.assertEqual(pca.process_warning[0], "Filtering removed all PCA attributes.")
