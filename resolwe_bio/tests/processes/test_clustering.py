# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class ClusteringProcessorTestCase(BioProcessTestCase):

    def test_hc_clustering(self):
        """Cannot use assertJSON - JSON output contains ETC object IDs."""
        expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM")
        expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")
        expression_3 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")

        inputs = {'expressions': [expression_1.pk, expression_2.pk, expression_3.pk]}
        etc = self.run_process('etc-bcm', inputs)

        inputs = {
            'etcs': [etc.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067098', 'DPU_G0067100']}
        self.run_process('clustering-hierarchical-genes-etc', inputs)

    def test_hc_clustering_samples_ucsc(self):
        inputs = {'exp': 'clustering_expressions_1.tab.gz',
                  'exp_type': 'Log2',
                  'exp_name': 'Expression',
                  'source': 'UCSC'}
        expression_1 = self.run_process('upload-expression', inputs)

        inputs = {'exp': 'clustering_expressions_2.tab.gz',
                  'exp_type': 'Log2',
                  'exp_name': 'Expression',
                  'source': 'UCSC'}
        expression_2 = self.run_process('upload-expression', inputs)

        inputs = {'exp': 'clustering_expressions_3.tab.gz',
                  'exp_type': 'Log2',
                  'exp_name': 'Expression',
                  'source': 'UCSC'}
        expression_3 = self.run_process('upload-expression', inputs)

        inputs = {'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
                  'genes': [' A1BG E2f4 A2ML1 A2MP1 A3GALT2 A4GALT AADACL4 AADACL3 AADACL2 AADAC'],
                  'genes_source': 'UCSC'}

        clustering = self.run_process('clustering-hierarchical-samples', inputs)

        saved_json, test_json = self.get_json('sample_cluster_data.json.gz', clustering.output['cluster'])
        self.assertEqual(test_json['linkage'], saved_json['linkage'])
        self.assertTrue('order' in test_json)

    def test_hc_clustering_genes_ucsc(self):
        inputs = {'exp': 'clustering_expressions_1.tab.gz',
                  'exp_type': 'Log2',
                  'exp_name': 'Expression',
                  'source': 'UCSC'}
        expression_1 = self.run_process('upload-expression', inputs)

        inputs = {'exp': 'clustering_expressions_2.tab.gz',
                  'exp_type': 'Log2',
                  'exp_name': 'Expression',
                  'source': 'UCSC'}
        expression_2 = self.run_process('upload-expression', inputs)

        inputs = {'exp': 'clustering_expressions_3.tab.gz',
                  'exp_type': 'Log2',
                  'exp_name': 'Expression',
                  'source': 'UCSC'}
        expression_3 = self.run_process('upload-expression', inputs)

        inputs = {'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
                  'genes': [' A1BG E2f4 A2ML1 A2MP1 A3GALT2 A4GALT AADACL4 AADACL3 AADACL2 AADAC'],
                  'genes_source': 'UCSC'}

        clustering = self.run_process('clustering-hierarchical-genes', inputs)

        saved_json, test_json = self.get_json('gene_cluster_data.json.gz', clustering.output['cluster'])
        self.assertEqual(test_json['linkage'], saved_json['linkage'])
        self.assertTrue('order' in test_json)

    def test_hc_clustering_genes_ncbi(self):
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

        inputs = {'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
                  'genes': ["1", "503538", "56934", "29974", "2", "144571", "3"],
                  'genes_source': 'NCBI'}

        clustering = self.run_process('clustering-hierarchical-genes', inputs)
        saved_json, test_json = self.get_json('gene_cluster_data_NCBI.json.gz', clustering.output['cluster'])
        self.assertEqual(test_json['linkage'], saved_json['linkage'])
        self.assertTrue('order' in test_json)

    def test_hc_clustering_samples_ncbi(self):
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

        inputs = {'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
                  'genes': ["1", "503538", "56934", "29974", "2", "144571", "3"],
                  'genes_source': 'NCBI'}

        clustering = self.run_process('clustering-hierarchical-samples', inputs)
        saved_json, test_json = self.get_json('sample_cluster_data_NCBI.json.gz', clustering.output['cluster'])
        self.assertEqual(test_json['linkage'], saved_json['linkage'])
        self.assertTrue('order' in test_json)
