# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class ClusteringProcessorTestCase(BioProcessTestCase):

    def test_hc_clustering_samples_ucsc(self):
        expression_1 = self.prepare_expression(f_exp='clustering_expressions_1.tab.gz',
                                               f_type='Log2',
                                               name='Expression',
                                               source='UCSC')
        expression_2 = self.prepare_expression(f_exp='clustering_expressions_2.tab.gz',
                                               f_type='Log2',
                                               name='Expression',
                                               source='UCSC')
        expression_3 = self.prepare_expression(f_exp='clustering_expressions_3.tab.gz',
                                               f_type='Log2',
                                               name='Expression',
                                               source='UCSC')

        inputs = {'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
                  'genes': ['A1BG', 'E2f4', 'A2ML1', 'A2MP1', 'A3GALT2', 'A4GALT',
                            'AADACL4', 'AADACL3', 'AADACL2', 'AADAC'],
                  'genes_source': 'UCSC'}

        clustering = self.run_process('clustering-hierarchical-samples', inputs)

        saved_json, test_json = self.get_json('sample_cluster_data.json.gz', clustering.output['cluster'])
        self.assertEqual(test_json['linkage'], saved_json['linkage'])
        self.assertTrue('order' in test_json)

    def test_hc_clustering_genes_ucsc(self):
        expression_1 = self.prepare_expression(f_exp='clustering_expressions_1.tab.gz',
                                               f_type='Log2',
                                               name='Expression',
                                               source='UCSC')
        expression_2 = self.prepare_expression(f_exp='clustering_expressions_2.tab.gz',
                                               f_type='Log2',
                                               name='Expression',
                                               source='UCSC')
        expression_3 = self.prepare_expression(f_exp='clustering_expressions_3.tab.gz',
                                               f_type='Log2',
                                               name='Expression',
                                               source='UCSC')

        inputs = {'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
                  'genes': ['A1BG', 'E2f4', 'A2ML1', 'A2MP1', 'A3GALT2', 'A4GALT',
                            'AADACL4', 'AADACL3', 'AADACL2', 'AADAC'],
                  'genes_source': 'UCSC'}

        clustering = self.run_process('clustering-hierarchical-genes', inputs)

        saved_json, test_json = self.get_json('gene_cluster_data.json.gz', clustering.output['cluster'])
        self.assertEqual(test_json['linkage'], saved_json['linkage'])
        self.assertTrue('order' in test_json)

    def test_hc_clustering_genes_ncbi(self):
        expression_1 = self.prepare_expression(f_exp='clustering_NCBI.gz',
                                               f_type='rc',
                                               name='Expression',
                                               source='NCBI')
        expression_2 = self.prepare_expression(f_exp='clustering_NCBI_1.gz',
                                               f_type='rc',
                                               name='Expression',
                                               source='NCBI')
        expression_3 = self.prepare_expression(f_exp='clustering_NCBI_2.gz',
                                               f_type='rc',
                                               name='Expression',
                                               source='NCBI')

        inputs = {'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
                  'genes': ['1', '503538', '56934', '29974', '2', '144571', '3'],
                  'genes_source': 'NCBI'}

        clustering = self.run_process('clustering-hierarchical-genes', inputs)
        saved_json, test_json = self.get_json('gene_cluster_data_NCBI.json.gz', clustering.output['cluster'])
        self.assertEqual(test_json['linkage'], saved_json['linkage'])
        self.assertTrue('order' in test_json)

    def test_hc_clustering_samples_ncbi(self):
        expression_1 = self.prepare_expression(f_exp='clustering_NCBI.gz',
                                               f_type='rc',
                                               name='Expression',
                                               source='NCBI')
        expression_2 = self.prepare_expression(f_exp='clustering_NCBI_1.gz',
                                               f_type='rc',
                                               name='Expression',
                                               source='NCBI')
        expression_3 = self.prepare_expression(f_exp='clustering_NCBI_2.gz',
                                               f_type='rc',
                                               name='Expression',
                                               source='NCBI')

        inputs = {'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
                  'genes': ['1', '503538', '56934', '29974', '2', '144571', '3'],
                  'genes_source': 'NCBI'}

        clustering = self.run_process('clustering-hierarchical-samples', inputs)
        saved_json, test_json = self.get_json('sample_cluster_data_NCBI.json.gz', clustering.output['cluster'])
        self.assertAlmostEqualGeneric(test_json['linkage'], saved_json['linkage'])
        self.assertTrue('order' in test_json)
