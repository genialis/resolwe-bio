# pylint: disable=missing-docstring
import six

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
        self.assertEqual(
            pca.process_warning[0],
            "Expressions of all selected genes are 0. Please select different samples or genes."
        )

    def test_pca_ncbi(self):
        inputs = {'exp': 'clustering_NCBI.gz',
                  'exp_type': 'rc',
                  'exp_name': 'Expression',
                  'source': 'NCBI'}
        expression_1 = self.run_process('upload-expression', inputs)

        inputs = {'exp': 'clustering_NCBI_1.gz',
                  'exp_type': 'rc',
                  'exp_name': 'Expression',
                  'source': 'NCBI'}
        expression_2 = self.run_process('upload-expression', inputs)

        inputs = {'exp': 'clustering_NCBI_2.gz',
                  'exp_type': 'rc',
                  'exp_name': 'Expression',
                  'source': 'NCBI'}
        expression_3 = self.run_process('upload-expression', inputs)

        inputs = {'exps': [expression_1.pk, expression_2.pk, expression_3.pk, ],
                  'genes': ["1", "503538", "56934", "29974", "2", "144571", "3", "abc", "lll"]}
        pca = self.run_process('pca', inputs)
        saved_json, test_json = self.get_json('pca_plot_ncbi.json.gz', pca.output['pca'])
        six.assertCountEqual(self, test_json['zero_gene_symbols'], saved_json['zero_gene_symbols'])
