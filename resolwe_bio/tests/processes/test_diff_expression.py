from pathlib import Path

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
            "create_sets": True,
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

        gene_set = Data.objects.filter(process__slug="upload-geneset").first()
        self.assertFile(
            gene_set,
            "geneset",
            "geneset_cuffdiff_all.tab.gz",
            compression="gzip",
        )
        self.assertFields(gene_set, "species", "Homo sapiens")
        self.assertFields(gene_set, "source", "UCSC")

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
            "create_sets": True,
            "filter_options": {
                "min_count_sum": 0,
            },
        }

        diff_exp = self.run_process("differentialexpression-deseq2", inputs)

        self.assertFileExists(diff_exp, "raw")
        self.assertFile(
            diff_exp, "count_matrix", "deseq2_count_matrix.tab.gz", compression="gzip"
        )
        self.assertFile(
            diff_exp,
            "count_matrix_normalized",
            "deseq2_count_matrix_normalized.tab.gz",
            compression="gzip",
        )
        self.assertJSON(diff_exp, diff_exp.output["de_json"], "", "deseq2.json.gz")
        self.assertFields(diff_exp, "source", "DICTYBASE")
        self.assertFields(diff_exp, "species", "Dictyostelium discoideum")
        self.assertFields(diff_exp, "build", "dd-05-2009")
        self.assertFields(diff_exp, "feature_type", "gene")

        gene_set = Data.objects.filter(process__slug="upload-geneset").last()
        self.assertFile(
            gene_set,
            "geneset",
            "geneset_deseq2_up.tab.gz",
            compression="gzip",
        )
        self.assertFields(gene_set, "species", "Dictyostelium discoideum")
        self.assertFields(gene_set, "source", "DICTYBASE")

        inputs["filter_options"]["independent"] = False
        diff_exp = self.run_process("differentialexpression-deseq2", inputs)
        self.assertFileExists(diff_exp, "raw")
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
    def test_deseq2_build(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(build="dd-05-2009")
            expression_2 = self.prepare_expression(build="dd-42-2009")

        inputs = {"case": [expression_1.pk], "control": [expression_2.pk]}

        deseq2 = self.run_process(
            "differentialexpression-deseq2", inputs, Data.STATUS_ERROR
        )

        error_msg = ["Input samples are of different Build: dd-42-2009 and dd-05-2009."]
        self.assertEqual(deseq2.process_error, error_msg)

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
            "case": [
                case.pk,
            ],
            "control": [
                control.pk,
            ],
        }
        deseq2 = self.run_process(
            "differentialexpression-deseq2", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc, : "
            "every gene contains at least one zero, cannot compute log geometric means"
        ]
        self.assertEqual(deseq2.process_error, error_msg)

    @with_resolwe_host
    @tag_process("differentialexpression-deseq2")
    def test_deseq2_nanostring(self):
        base = Path("test_nanostring_deseq2")
        inputs = base / "inputs"
        outputs = base / "outputs"
        with self.preparation_stage():
            expression = Process.objects.create(
                name="Upload nanostring expression data mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {
                        "network": True,
                    },
                    "executor": {
                        "docker": {
                            "image": "public.ecr.aws/genialis/resolwebio/base:ubuntu-22.04-14112023",
                        },
                    },
                },
                contributor=self.contributor,
                type="data:expression:nanostring:",
                entity_type="sample",
                data_name="{{ src.file }}",
                input_schema=[
                    {
                        "name": "src",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "exp",
                        "type": "basic:file:",
                    },
                    {
                        "name": "species",
                        "type": "basic:string:",
                    },
                    {
                        "name": "build",
                        "type": "basic:string:",
                    },
                    {
                        "name": "source",
                        "type": "basic:string:",
                    },
                    {
                        "name": "feature_type",
                        "type": "basic:string:",
                    },
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "txt" "txt" 0.1 compress
re-save-file exp "${NAME}".txt.gz
re-save species "Dictyostelium discoideum"
re-save build "dd-05-2009"
re-save source "DICTYBASE"
re-save feature_type "gene"
""",
                },
            )

            exp_1 = self.run_process(
                expression.slug,
                {"src": str(inputs / "exp_1_norm.txt.gz")},
            )
            exp_2 = self.run_process(
                expression.slug,
                {"src": str(inputs / "exp_2_norm.txt.gz")},
            )
            exp_3 = self.run_process(
                expression.slug,
                {"src": str(inputs / "exp_3_norm.txt.gz")},
            )
            exp_4 = self.run_process(
                expression.slug,
                {"src": str(inputs / "exp_4_norm.txt.gz")},
            )

        inputs = {
            "case": [exp_1.id, exp_3.id],
            "control": [exp_2.id, exp_4.id],
            "filter_options": {
                "min_count_sum": 0,
            },
        }

        diff_exp = self.run_process("differentialexpression-deseq2", inputs)

        self.assertFileExists(diff_exp, "raw")
        self.assertFile(
            diff_exp,
            "count_matrix",
            str(outputs / "count_matrix.tab.gz"),
            compression="gzip",
        )
        self.assertJSON(
            diff_exp,
            diff_exp.output["de_json"],
            "",
            str(outputs / "de_data_deseq.json.gz"),
        )
        self.assertFields(diff_exp, "source", "DICTYBASE")
        self.assertFields(diff_exp, "species", "Dictyostelium discoideum")
        self.assertFields(diff_exp, "build", "dd-05-2009")
        self.assertFields(diff_exp, "feature_type", "gene")

    @with_resolwe_host
    @tag_process("differentialexpression-deseq2", "differentialexpression-edger")
    def test_de_microarray(self):
        with self.preparation_stage():
            expression = Process.objects.create(
                name="Upload microarray expression data mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {
                        "network": True,
                    },
                    "executor": {
                        "docker": {
                            "image": "public.ecr.aws/genialis/resolwebio/base:ubuntu-22.04-14112023",
                        },
                    },
                },
                contributor=self.contributor,
                type="data:expression:microarray:",
                entity_type="sample",
                data_name="{{ name }}",
                input_schema=[
                    {
                        "name": "name",
                        "type": "basic:string:",
                    },
                ],
                output_schema=[
                    {
                        "name": "name",
                        "type": "basic:string:",
                    },
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-save name name
""",
                },
            )

            exp_1 = self.run_process(
                expression.slug,
                {"name": "mock_ma_expression"},
            )

        inputs = {
            "case": [exp_1.id],
            "control": [exp_1.id],
        }

        deseq2 = self.run_process(
            "differentialexpression-deseq2", inputs, Data.STATUS_ERROR
        )
        edger = self.run_process(
            "differentialexpression-edger", inputs, Data.STATUS_ERROR
        )

        error_msg = ["Microarray expressions are not supported."]
        self.assertEqual(deseq2.process_error, error_msg)
        self.assertEqual(edger.process_error, error_msg)

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
            "create_sets": True,
        }

        diff_exp = self.run_process("differentialexpression-edger", inputs)
        self.assertFile(
            obj=diff_exp,
            field_path="raw",
            fn="diffexp_edgeR.tab.gz",
            compression="gzip",
        )

        self.assertJSON(diff_exp, diff_exp.output["de_json"], "", "edgeR.json.gz")
        self.assertFields(diff_exp, "source", "DICTYBASE")
        self.assertFields(diff_exp, "species", "Dictyostelium discoideum")
        self.assertFields(diff_exp, "build", "dd-05-2009")
        self.assertFields(diff_exp, "feature_type", "gene")

        gene_set = Data.objects.filter(process__slug="upload-geneset").last()
        self.assertFile(
            gene_set,
            "geneset",
            "geneset_edgeR_up.tab.gz",
            compression="gzip",
        )
        self.assertFields(gene_set, "species", "Dictyostelium discoideum")
        self.assertFields(gene_set, "source", "DICTYBASE")
