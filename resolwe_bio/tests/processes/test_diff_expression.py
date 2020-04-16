from resolwe.flow.models import Data, Process
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class DiffExpProcessorTestCase(KBBioProcessTestCase):
    @tag_process("cuffdiff")
    def test_cuffdiff(self):
        with self.preparation_stage():
            inputs = {
                "src": "cuffquant 1.cxb",
                "source": "UCSC",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            cuffquant = self.run_process("upload-cxb", inputs)

            inputs = {
                "src": "cuffquant_2.cxb",
                "source": "UCSC",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            cuffquant2 = self.run_process("upload-cxb", inputs)

            annotation = self.prepare_annotation(
                fn="hg19_chr20_small.gtf.gz", source="UCSC"
            )

        inputs = {
            "case": [cuffquant.id],
            "control": [cuffquant2.id],
            "annotation": annotation.id,
        }
        cuffdiff = self.run_process("cuffdiff", inputs)
        self.assertFile(cuffdiff, "raw", "raw_cuffdiff.tab.gz", compression="gzip")
        self.assertFile(
            cuffdiff, "de_file", "de_file_cuffdiff.tab.gz", compression="gzip"
        )
        self.assertJSON(cuffdiff, cuffdiff.output["de_json"], "", "cuffdiff.json.gz")
        self.assertFields(cuffdiff, "source", "UCSC")
        self.assertFields(cuffdiff, "species", "Homo sapiens")
        self.assertFields(cuffdiff, "build", "hg19")
        self.assertFields(cuffdiff, "feature_type", "gene")

    @with_resolwe_host
    @tag_process("differentialexpression-deseq2")
    def test_deseq2_genes(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(
                f_rc="exp_1_rc.tab.gz", source="DICTYBASE"
            )
            expression_2 = self.prepare_expression(
                f_rc="exp_2_rc.tab.gz", source="DICTYBASE"
            )
            expression_3 = self.prepare_expression(
                f_rc="exp_3_rc.tab.gz", source="DICTYBASE"
            )
            expression_4 = self.prepare_expression(
                f_rc="exp 4_rc.tab.gz", source="DICTYBASE"
            )

        inputs = {
            "case": [expression_1.pk, expression_3.pk],
            "control": [expression_2.pk, expression_4.pk],
            "filter": {"min_count_sum": 0,},
        }

        diff_exp = self.run_process("differentialexpression-deseq2", inputs)

        self.assertFileExists(diff_exp, "raw")
        self.assertFile(
            diff_exp, "count_matrix", "deseq2_count_matrix.tab.gz", compression="gzip"
        )
        self.assertJSON(diff_exp, diff_exp.output["de_json"], "", "deseq2.json.gz")
        self.assertFields(diff_exp, "source", "DICTYBASE")
        self.assertFields(diff_exp, "species", "Dictyostelium discoideum")
        self.assertFields(diff_exp, "build", "dd-05-2009")
        self.assertFields(diff_exp, "feature_type", "gene")

    @with_resolwe_host
    @tag_process("differentialexpression-deseq2")
    def test_deseq2_source(self):
        with self.preparation_stage():
            expression_dictybase = self.prepare_expression(source="DICTYBASE")
            expression_ucsc = self.prepare_expression(source="UCSC")

        inputs = {"case": [expression_dictybase.pk], "control": [expression_ucsc.pk]}

        self.run_process("differentialexpression-deseq2", inputs, Data.STATUS_ERROR)

    @with_resolwe_host
    @tag_process("differentialexpression-deseq2")
    def test_deseq2_expression_error(self):
        with self.preparation_stage():
            case = self.prepare_expression(
                f_rc="deseq2_exp1.tab.gz",
                source="UCSC",
                species="Mus musculus",
                build="mm10",
            )
            control = self.prepare_expression(
                f_rc="deseq2_exp2.tab.gz",
                source="UCSC",
                species="Mus musculus",
                build="mm10",
            )

        inputs = {
            "case": [case.pk,],
            "control": [control.pk,],
        }
        deseq2 = self.run_process(
            "differentialexpression-deseq2", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            (
                "Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc, : "
                "every gene contains at least one zero, cannot compute log geometric means"
            )
        ]
        self.assertEqual(deseq2.process_error, error_msg)

    @with_resolwe_host
    @tag_process("differentialexpression-edger")
    def test_edger(self):
        with self.preparation_stage():
            inputs = {
                "rc": "exp_1_rc.tab.gz",
                "exp_name": "Expression",
                "source": "DICTYBASE",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            expression_1 = self.run_process("upload-expression", inputs)

            inputs = {
                "rc": "exp_2_rc.tab.gz",
                "exp_name": "Expression",
                "source": "DICTYBASE",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            expression_2 = self.run_process("upload-expression", inputs)

            inputs = {
                "rc": "exp_3_rc.tab.gz",
                "exp_name": "Expression",
                "source": "DICTYBASE",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            expression_3 = self.run_process("upload-expression", inputs)

            inputs = {
                "rc": "exp_4_rc.tab.gz",
                "exp_name": "Expression",
                "source": "DICTYBASE",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            expression_4 = self.run_process("upload-expression", inputs)

        inputs = {
            "case": [expression_1.id, expression_3.id],
            "control": [expression_2.id, expression_4.id],
        }

        diff_exp = self.run_process("differentialexpression-edger", inputs)
        self.assertFile(diff_exp, "raw", "diffexp_edgeR.tab.gz", compression="gzip")
        self.assertJSON(diff_exp, diff_exp.output["de_json"], "", "edgeR.json.gz")
        self.assertFields(diff_exp, "source", "DICTYBASE")
        self.assertFields(diff_exp, "species", "Dictyostelium discoideum")
        self.assertFields(diff_exp, "build", "dd-05-2009")
        self.assertFields(diff_exp, "feature_type", "gene")

    @with_resolwe_host
    @tag_process("differentialexpression-shrna")
    def test_shrna_diffexp(self):
        with self.preparation_stage():
            pf_in = "shrna_diffexp/input/"
            pf_out = "shrna_diffexp/output/"

            # Prepare parameter file.
            parameter_file = self.run_process(
                "upload-file", {"src": pf_in + "template_doDE_inputs.xlsx"}
            )

            # Prepare mock upload process.
            process = Process.objects.create(
                name="Upload shRNA expression files produced by shrna-quant.",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {"network": True,},
                    "executor": {"docker": {"image": "resolwebio/base:ubuntu-18.04",},},
                },
                data_name="shRNA expression ({{ reads|sample_name|default('?') }})",
                # data_name="shRNA expression ({{ input_data | name | default('?') }})",
                contributor=self.contributor,
                type="data:expression:shrna2quant:",
                input_schema=[{"name": "exp", "type": "basic:file:",}],
                output_schema=[
                    {"name": "exp", "type": "basic:file:",},
                    {"name": "rc", "type": "basic:file:"},
                    {"name": "exp_type", "type": "basic:string:"},
                    {"name": "source", "type": "basic:string:"},
                    {"name": "species", "type": "basic:string:"},
                    {"name": "build", "type": "basic:string:"},
                    {"name": "feature_type", "type": "basic:string:"},
                    {"name": "mapped_species", "type": "basic:file:"},
                ],
                run={
                    "language": "bash",
                    "program": r"""
            re-import {{ exp.file_temp|default(exp.file) }} {{ exp.file }} "txt" "txt" 0.1 compress

            cp "${NAME}".txt.gz "${NAME}"_trimmed_trimmed_count_matrix.txt.gz

            re-save-file exp "${NAME}"_trimmed_trimmed_count_matrix.txt.gz
            re-save-file rc "${NAME}"_trimmed_trimmed_count_matrix.txt.gz
            re-save exp_type 'RC'
            re-save source 'shRNA-gene-sequences'
            re-save species 'Homo sapiens'
            re-save build 'custom-from-file'
            re-save feature_type 'shRNA'
            # Uploading fake mapped species file, do not use in tests.
            re-save-file mapped_species "${NAME}".txt.gz
            """,
                },
            )

            sample_filenames = [
                pf_in + f"sample-{i}.txt.gz" for i in list(range(1, 13))
            ]

            list_expressions = []
            for i in sample_filenames:
                list_expressions.append(self.run_process(process.slug, {"exp": i}))

        inputs = {
            "parameter_file": parameter_file.id,
            "expression_data": [i.pk for i in list_expressions],
        }

        de_res = self.run_process("differentialexpression-shrna", inputs)

        self.assertFile(
            de_res, "deseq_results", pf_out + "deseq_results.txt.gz", compression="gzip"
        )
        self.assertFile(
            de_res, "class_results", pf_out + "class_results.txt.gz", compression="gzip"
        )
        self.assertFile(
            de_res,
            "beneficial_counts",
            pf_out + "beneficial_counts.txt.gz",
            compression="gzip",
        )
        self.assertFile(
            de_res, "lethal_counts", pf_out + "lethal_counts.txt.gz", compression="gzip"
        )
