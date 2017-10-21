# pylint: disable=missing-docstring
from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase


class ExpressionAggregatorTestCase(BioProcessTestCase):
    @tag_process('expression-aggregator')
    def test_expression_aggregator(self):
        with self.preparation_stage():
            descriptor_hs = {'sample': {'organism': 'Homo sapiens'}}
            expression1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz',
                                                  f_type='TPM', descriptor=descriptor_hs,
                                                  species='Homo sapiens', build='ens_90')
            expression2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz',
                                                  f_type='TPM', descriptor=descriptor_hs,
                                                  species='Homo sapiens', build='ens_90')

            descriptor_mm = {'sample': {'organism': 'Mus musculus'}}
            expression3 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz',
                                                  f_type='TPM', descriptor=descriptor_mm,
                                                  species='Mus musculus', build='ens_90')
            expression4 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz',
                                                  f_type='TPM', descriptor=descriptor_mm,
                                                  species='Mus musculus', build='ens_90')

        inputs = {
            'exps': [expression1.id, expression2.id, expression3.id, expression4.id],
            'group_by': 'sample.organism'
        }
        expression_aggregator = self.run_process('expression-aggregator', inputs)

        saved_matrix, test_matrix = self.get_json('exp_matrix.json.gz', expression_aggregator.output['exp_matrix'])
        self.assertEqual(saved_matrix, test_matrix)
        self.assertJSON(
            expression_aggregator, expression_aggregator.output['box_plot'], '', 'box_plot.json.gz'
        )
        self.assertJSON(
            expression_aggregator, expression_aggregator.output['log_box_plot'], '', 'log_box_plot.json.gz'
        )
        self.assertFields(expression_aggregator, 'source', 'DICTYBASE')
        self.assertFields(expression_aggregator, 'exp_type', 'TPM')
