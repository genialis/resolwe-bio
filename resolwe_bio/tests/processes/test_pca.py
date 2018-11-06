# pylint: disable=missing-docstring
import six

from resolwe.test import tag_process
from resolwe_bio.utils.test import with_resolwe_host, KBBioProcessTestCase


class PcaProcessorTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process('pca')
    def test_pca(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz',
                                                   f_exp='exp_1_tpm.tab.gz',
                                                   f_type='TPM',
                                                   source='DICTYBASE',
                                                   species='Dictyostelium discoideum')
            expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz',
                                                   f_exp='exp_2_tpm.tab.gz',
                                                   f_type='TPM',
                                                   source='DICTYBASE',
                                                   species='Dictyostelium discoideum')
            expression_noname_1 = self.prepare_expression(
                f_rc='pca_exp_noname1.tab.gz',
                f_exp='pca_exp_noname1.tab.gz',
                f_type='TPM',
                source='DICTYBASE',
                species='Dictyostelium discoideum',
            )
            expression_noname_2 = self.prepare_expression(
                f_rc='pca_exp_noname2.tab.gz',
                f_exp='pca_exp_noname2.tab.gz',
                f_type='TPM',
                source='DICTYBASE',
                species='Dictyostelium discoideum',
            )

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'source': 'DICTYBASE',
            'species': 'Dictyostelium discoideum'
        }
        pca = self.run_process('pca', inputs)
        saved_json, test_json = self.get_json('pca_plot.json.gz', pca.output['pca'])
        self.assertAlmostEqualGeneric(test_json['flot']['data'], saved_json['flot']['data'])
        self.assertAlmostEqualGeneric(test_json['explained_variance_ratios'], saved_json['explained_variance_ratios'])
        self.assertAlmostEqualGeneric(test_json['components'], saved_json['components'])
        self.assertEqual(len(pca.process_warning), 0)

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': ['DPU_G0067098', 'DPU_G0067100', 'DPU_G0067104'],  # all zero
            'source': 'DICTYBASE',
            'species': 'Dictyostelium discoideum'
        }
        pca = self.run_process('pca', inputs)
        saved_json, test_json = self.get_json('pca_plot2.json.gz', pca.output['pca'])
        self.assertAlmostEqualGeneric(test_json['flot']['data'], saved_json['flot']['data'])
        self.assertAlmostEqualGeneric(test_json['explained_variance_ratios'], saved_json['explained_variance_ratios'])
        self.assertAlmostEqualGeneric(test_json['components'], saved_json['components'])
        self.assertEqual(len(pca.process_warning), 0)

        inputs = {
            'exps': [expression_1.pk],
            'source': 'DICTYBASE',
            'species': 'Dictyostelium discoideum'
        }
        pca = self.run_process('pca', inputs)
        saved_json, test_json = self.get_json('pca_plot_single_sample.json.gz', pca.output['pca'])
        self.assertAlmostEqualGeneric(test_json['flot']['data'], saved_json['flot']['data'])
        self.assertAlmostEqualGeneric(test_json['explained_variance_ratios'], saved_json['explained_variance_ratios'])
        self.assertAlmostEqualGeneric(test_json['components'], saved_json['components'])
        self.assertEqual(len(pca.process_warning), 0)

        inputs = {
            'exps': [
                expression_noname_1.pk,
                expression_noname_2.pk,
            ],
            'source': 'DICTYBASE',
            'species': 'Dictyostelium discoideum',
        }
        pca = self.run_process('pca', inputs)

    @with_resolwe_host
    @tag_process('pca')
    def test_pca_ncbi(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(f_exp='clustering_NCBI.gz',
                                                   f_type='rc',
                                                   name='Expression',
                                                   source='NCBI',
                                                   species='Homo sapiens')
            expression_2 = self.prepare_expression(f_exp='clustering_NCBI_1.gz',
                                                   f_type='rc',
                                                   name='Expression',
                                                   source='NCBI',
                                                   species='Homo sapiens')
            expression_3 = self.prepare_expression(f_exp='clustering_NCBI_2.gz',
                                                   f_type='rc',
                                                   name='Expression',
                                                   source='NCBI',
                                                   species='Homo sapiens')

        inputs = {
            'exps': [expression_1.pk, expression_2.pk, expression_3.pk, ],
            'genes': ['1', '503538', '56934', '29974', '2', '144571', '3', 'abc', 'lll'],
            'source': 'NCBI',
            'species': 'Homo sapiens'
        }
        pca = self.run_process('pca', inputs)
        saved_json, test_json = self.get_json('pca_plot_ncbi.json.gz', pca.output['pca'])
        six.assertCountEqual(self, test_json['zero_gene_symbols'], saved_json['zero_gene_symbols'])
