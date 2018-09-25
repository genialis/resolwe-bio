# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process
from resolwe_bio.utils.test import with_resolwe_host, KBBioProcessTestCase


class ExpressionAggregatorTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process('expression-aggregator')
    def test_expression_aggregator(self):
        with self.preparation_stage():
            descriptor_artery = {'general': {'organ': 'artery'}}
            expression1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz',
                                                  f_type='TPM', descriptor=descriptor_artery,
                                                  species='Homo sapiens', build='ens_90')

            descriptor_blood = {'general': {'organ': 'blood'}}
            expression2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz',
                                                  f_type='TPM', descriptor=descriptor_blood,
                                                  species='Homo sapiens', build='ens_90')

            expression3 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz',
                                                  f_type='TPM', descriptor=descriptor_artery,
                                                  species='Mus musculus', build='ens_90')

            expression4 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz',
                                                  f_type='TPM', descriptor=descriptor_blood,
                                                  species='Mus musculus', build='ens_90')

        # Expect the process to fail if expression data from multiple species is used
        inputs = {
            'exps': [expression1.id, expression2.id, expression3.id],
            'group_by': 'general.organ',
        }
        expression_aggregator = self.run_process('expression-aggregator', inputs, Data.STATUS_ERROR)

        inputs['exps'] = [expression3.id, expression4.id]
        expression_aggregator_mm = self.run_process('expression-aggregator', inputs)

        inputs['exps'] = [expression1.id, expression2.id]
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
        self.assertFields(expression_aggregator, 'species', 'Homo sapiens')

        # Do not allow appending expression data of one species to the existing Aggregator objects of different species
        inputs['expr_aggregator'] = expression_aggregator_mm.id
        expression_aggregator = self.run_process('expression-aggregator', inputs, Data.STATUS_ERROR)
