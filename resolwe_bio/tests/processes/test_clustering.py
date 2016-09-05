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
        self.run_process('clustering-hierarchical-bcm', inputs)

    def test_hc_clustering_samples(self):
        inputs = {'exp': 'hs_expressions_1.tab.gz', 'exp_type': 'Log2', 'exp_name': 'Expression', 'source': 'hg19'}
        expression_1 = self.run_process('upload-expression', inputs)

        inputs = {'exp': 'hs_expressions_2.tab.gz', 'exp_type': 'Log2', 'exp_name': 'Expression', 'source': 'hg19'}
        expression_2 = self.run_process('upload-expression', inputs)

        inputs = {'exp': 'hs_expressions_3.tab.gz', 'exp_type': 'Log2', 'exp_name': 'Expression', 'source': 'hg19'}
        expression_3 = self.run_process('upload-expression', inputs)

        inputs = {'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
                  'sampleids': ["1", "2", "3"]}
        clustering = self.run_process('clustering-hierarchical-samples', inputs)
        self.assertJSON(clustering, clustering.output['cluster'], 'linkage', 'sample_custer_data.json.gz')
        self.assertJSON(clustering, clustering.output['cluster'], 'samples_names', 'sample_custer_names.json.gz')
