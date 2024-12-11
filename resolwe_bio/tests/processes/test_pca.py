from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


def round_elements(element, places=2):
    """
    This function takes a list of arbitrary depth and rounds all numerical elements to the specified places.
    Used due to numerical differences in the output of the PCA process between different CPU architectures.

    :param element: a nested list of arbitrary depth
    :param places: number of decimal places to round to
    :return: a new nested list with all numerical elements rounded
    """
    if isinstance(element, list):
        return [round_elements(sublist) for sublist in element]
    elif isinstance(element, (int, float)):
        return round(element, places)
    else:
        return element


class PcaProcessorTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("pca")
    def test_pca(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(
                f_rc="exp_1_rc.tab.gz",
                f_exp="exp_1_tpm.tab.gz",
                f_type="TPM",
                source="DICTYBASE",
                species="Dictyostelium discoideum",
            )
            expression_2 = self.prepare_expression(
                f_rc="exp_2_rc.tab.gz",
                f_exp="exp_2_tpm.tab.gz",
                f_type="TPM",
                source="DICTYBASE",
                species="Dictyostelium discoideum",
            )
            expression_3 = self.prepare_expression(
                f_exp="exp_1_tpm.tab.gz",
                f_rc=None,
                f_type="TPM",
                source="DICTYBASE",
                species="Dictyostelium discoideum",
            )
            expression_4 = self.prepare_expression(
                f_exp="exp_2_tpm.tab.gz",
                f_rc=None,
                f_type="TPM",
                source="DICTYBASE",
                species="Dictyostelium discoideum",
            )
            expression_5 = self.prepare_expression(
                f_rc="exp_6_rc.tab.gz",
                f_exp="exp_6_tpm.tab.gz",
                f_type="TPM",
                source="DICTYBASE",
                species="Dictyostelium discoideum",
            )
            expression_6 = self.prepare_expression(
                f_rc="exp_7_rc.tab.gz",
                f_exp="exp_7_tpm.tab.gz",
                f_type="TPM",
                source="DICTYBASE",
                species="Dictyostelium discoideum",
            )
            expression_7 = self.prepare_expression(
                f_rc="exp_8_rc.tab.gz",
                f_exp="exp_8_tpm.tab.gz",
                f_type="TPM",
                source="DICTYBASE",
                species="Dictyostelium discoideum",
            )

        inputs = {
            "exps": [
                expression_1.pk,
                expression_2.pk,
                expression_5.pk,
                expression_6.pk,
                expression_7.pk,
            ],
            "source": "DICTYBASE",
            "species": "Dictyostelium discoideum",
        }
        pca = self.run_process("pca", inputs)
        saved_json, test_json = self.get_json("pca_plot.json.gz", pca.output["pca"])
        # returns 4 PCA components. Last component differs when testing on different systems
        # and is not tested here.
        self.assertAlmostEqualGeneric(
            round_elements(test_json["flot"]["data"][0][:3]),
            round_elements(saved_json["flot"]["data"][0][:3]),
        )
        self.assertAlmostEqualGeneric(
            round_elements(test_json["flot"]["data"][1][:3]),
            round_elements(saved_json["flot"]["data"][1][:3]),
        )
        self.assertAlmostEqualGeneric(
            round_elements(test_json["flot"]["data"][2][:3]),
            round_elements(saved_json["flot"]["data"][2][:3]),
        )
        self.assertAlmostEqualGeneric(
            round_elements(test_json["flot"]["data"][3][:3]),
            round_elements(saved_json["flot"]["data"][3][:3]),
        )
        self.assertAlmostEqualGeneric(
            round_elements(test_json["flot"]["data"][4][:3]),
            round_elements(saved_json["flot"]["data"][4][:3]),
        )
        self.assertAlmostEqualGeneric(
            round_elements(test_json["explained_variance_ratios"]),
            round_elements(saved_json["explained_variance_ratios"]),
        )
        self.assertAlmostEqualGeneric(
            round_elements(test_json["components"]),
            round_elements(saved_json["components"]),
        )
        self.assertEqual(len(pca.process_warning), 0)

        inputs["genes"] = [
            "DPU_G0067096",
            "DPU_G0067106",
            "DPU_G0067102",
            "DPU_G0067112",
        ]  # all non-zero
        pca = self.run_process("pca", inputs)
        saved_json, test_json = self.get_json("pca_plot_2.json.gz", pca.output["pca"])
        self.assertAlmostEqualGeneric(
            round_elements(test_json["flot"]["data"][0][:3]),
            round_elements(saved_json["flot"]["data"][0][:3]),
        )

        self.assertEqual(len(pca.process_warning), 0)

        inputs["low_expression_filter"] = False
        inputs["standard_scaler"] = False
        pca = self.run_process("pca", inputs)
        saved_json, test_json = self.get_json("pca_plot_3.json.gz", pca.output["pca"])
        self.assertAlmostEqualGeneric(
            round_elements(test_json["flot"]["data"][0][:3]),
            round_elements(saved_json["flot"]["data"][0][:3]),
        )

        self.assertEqual(len(pca.process_warning), 0)

        inputs = {
            "exps": [expression_3.pk, expression_4.pk],
            "source": "DICTYBASE",
            "species": "Dictyostelium discoideum",
        }
        pca = self.run_process("pca", inputs)
        saved_json, test_json = self.get_json("pca_plot_4.json.gz", pca.output["pca"])
        self.assertAlmostEqualGeneric(
            round_elements(test_json["flot"]["data"]),
            round_elements(saved_json["flot"]["data"]),
        )

        self.assertEqual(len(pca.process_warning), 1)
