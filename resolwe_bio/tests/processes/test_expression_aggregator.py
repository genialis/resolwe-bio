from resolwe.flow.models import Data
from resolwe.flow.models.annotations import (
    AnnotationField,
    AnnotationGroup,
    AnnotationValue,
)
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.models import Sample
from resolwe_bio.utils.test import KBBioProcessTestCase


class ExpressionAggregatorTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("expression-aggregator")
    def test_expression_aggregator(self):
        with self.preparation_stage():
            ann_group = AnnotationGroup.objects.get(name="biospecimen_information")
            ann_field = AnnotationField.objects.get(group=ann_group, name="organ")
            expression1 = self.prepare_expression(
                f_rc="exp_1_rc.tab.gz",
                f_exp="exp_1_tpm.tab.gz",
                f_type="TPM",
                species="Homo sapiens",
                build="ens_90",
            )
            entity = Sample.objects.get(data=expression1)
            AnnotationValue.objects.create(
                entity=entity, field=ann_field, value="artery"
            )

            expression2 = self.prepare_expression(
                f_rc="exp_2_rc.tab.gz",
                f_exp="exp_2_tpm.tab.gz",
                f_type="TPM",
                species="Homo sapiens",
                build="ens_90",
            )
            entity = Sample.objects.get(data=expression2)
            AnnotationValue.objects.create(
                entity=entity, field=ann_field, value="blood"
            )

            expression3 = self.prepare_expression(
                f_rc="exp_1_rc.tab.gz",
                f_exp="exp_1_tpm.tab.gz",
                f_type="TPM",
                species="Mus musculus",
                build="ens_90",
            )
            entity = Sample.objects.get(data=expression3)
            AnnotationValue.objects.create(
                entity=entity, field=ann_field, value="artery"
            )

            expression4 = self.prepare_expression(
                f_rc="exp_2_rc.tab.gz",
                f_exp="exp_2_tpm.tab.gz",
                f_type="TPM",
                species="Mus musculus",
                build="ens_90",
            )
            entity = Sample.objects.get(data=expression4)
            AnnotationValue.objects.create(
                entity=entity, field=ann_field, value="blood"
            )

        # Expect the process to fail if expression data from multiple species is used
        inputs = {
            "exps": [expression1.id, expression2.id, expression3.id],
            "group_by": "biospecimen_information.organ",
        }
        expression_aggregator = self.run_process(
            "expression-aggregator", inputs, Data.STATUS_ERROR
        )

        inputs["exps"] = [expression3.id, expression4.id]
        expression_aggregator_mm = self.run_process("expression-aggregator", inputs)

        # List of samples with unordered descriptors.
        inputs["exps"] = [expression2.id, expression1.id]
        expression_aggregator = self.run_process("expression-aggregator", inputs)

        saved_matrix, test_matrix = self.get_json(
            "exp_matrix.json.gz", expression_aggregator.output["exp_matrix"]
        )
        self.assertEqual(saved_matrix, test_matrix)
        self.assertJSON(
            expression_aggregator,
            expression_aggregator.output["box_plot"],
            "",
            "box_plot.json.gz",
        )
        self.assertJSON(
            expression_aggregator,
            expression_aggregator.output["log_box_plot"],
            "",
            "log_box_plot.json.gz",
        )
        self.assertFields(expression_aggregator, "source", "DICTYBASE")
        self.assertFields(expression_aggregator, "exp_type", "TPM")
        self.assertFields(expression_aggregator, "species", "Homo sapiens")

        # Do not allow appending expression data of one species to the existing Aggregator objects of different species
        inputs["expr_aggregator"] = expression_aggregator_mm.id
        expression_aggregator = self.run_process(
            "expression-aggregator", inputs, Data.STATUS_ERROR
        )
