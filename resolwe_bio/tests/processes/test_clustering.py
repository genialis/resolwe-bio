from pathlib import Path

from guardian.shortcuts import assign_perm

from resolwe.flow.models import Data, Entity, Relation, RelationPartition, RelationType
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class ClusteringProcessTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("clustering-hierarchical-samples")
    def test_sample_clustering(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(
                f_exp="clustering_ensembl_1.tab.gz",
                f_type="TPM",
                name="expression",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            wrong_source = self.prepare_expression(
                f_exp="clustering_ensembl_1.tab.gz",
                f_type="TPM",
                name="wrong source",
                source="UCSC",
                species="Homo sapiens",
            )
            wrong_species = self.prepare_expression(
                f_exp="clustering_ensembl_1.tab.gz",
                f_type="TPM",
                name="wrong species",
                source="ENSEMBL",
                species="Mus musculus",
            )
            wrong_expression_type = self.prepare_expression(
                f_exp="clustering_ensembl_1.tab.gz",
                f_type="FPKM",
                name="wrong expression type",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            wrong_feature_type = self.prepare_expression(
                f_exp="clustering_ensembl_1.tab.gz",
                f_type="TPM",
                feature_type="transcript",
                name="wrong feature type",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_2 = self.prepare_expression(
                f_exp="clustering_ensembl_2.tab.gz",
                f_type="TPM",
                name="expression 2",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_3 = self.prepare_expression(
                f_exp="clustering_ensembl_3.tab.gz",
                f_type="TPM",
                name="expression 3",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_5 = self.prepare_expression(
                f_exp="clustering_ensembl_5.tab.gz",
                f_type="TPM",
                name="expression 5",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_6 = self.prepare_expression(
                f_exp="clustering_ensembl_6.tab.gz",
                f_type="TPM",
                name="expression 6",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_7 = self.prepare_expression(
                f_exp="clustering_ensembl_7.tab.gz",
                f_type="TPM",
                name="expression 7",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_8 = self.prepare_expression(
                f_exp="clustering_int1.tab.gz",
                f_type="TPM",
                name="expression int 1",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_9 = self.prepare_expression(
                f_exp="clustering_int2.tab.gz",
                f_type="TPM",
                name="expression int 2",
                source="ENSEMBL",
                species="Homo sapiens",
            )

        inputs = {
            "exps": [
                expression_8.pk,
                expression_9.pk,
            ],
            "preprocessing": {
                "log2": False,
            },
        }
        clustering = self.run_process("clustering-hierarchical-samples", inputs)
        saved_json, test_json = self.get_json(
            "clustering_out_sample.json.gz", clustering.output["cluster"]
        )
        # correct sample IDs
        saved_json["sample_ids"] = {
            str(i): {
                "id": Entity.objects.get(
                    data=Data.objects.get(id=inputs["exps"][i])
                ).id,
            }
            for i in range(len(inputs["exps"]))
        }
        self.assertEqual(saved_json.keys(), test_json.keys())
        for key in saved_json:
            assert_method = (
                self.assertAlmostEqualGeneric if key == "linkage" else self.assertEqual
            )
            assert_method(saved_json[key], test_json[key])
        self.assertEqual(len(clustering.process_warning), 0)

        inputs = {
            "exps": [
                expression_1.pk,
                wrong_source.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        warning_msg = [
            "All expression data must be annotated by the same genome database."
        ]
        error_msg = [
            (
                "Sample 'expression' has 'ENSEMBL' gene IDs, "
                "while sample 'wrong source' has 'UCSC' gene IDs."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                wrong_species.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        warning_msg = ["All expressions must be of the same Species."]
        error_msg = [
            (
                "Sample 'expression' is 'Homo sapiens', "
                "while sample 'wrong species' is 'Mus musculus'."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                wrong_expression_type.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        warning_msg = ["All expressions must be of the same Expression type."]
        error_msg = [
            (
                "Expression 'expression' has 'TPM' expression type, "
                "while sample 'wrong expression type' has 'FPKM' expression type."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                wrong_feature_type.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        warning_msg = ["All expressions must be of the same Feature type."]
        error_msg = [
            (
                "Expression 'expression' has 'gene' feature type, "
                "while sample 'wrong feature type' has 'transcript' feature type."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
            ],
            "preprocessing": {
                "genes": ["gene"],
                "source": "UCSC",
                "species": "Homo sapiens",
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        warning_msg = [
            "Selected genes must be annotated by the same genome database as all expression files."
        ]
        error_msg = [
            (
                "Gene IDs are from 'UCSC' database, "
                "while sample 'expression' has gene IDs from 'ENSEMBL' database."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
            ],
            "preprocessing": {
                "genes": ["gene"],
                "source": "ENSEMBL",
                "species": "Mus musculus",
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        warning_msg = [
            "Selected genes must be from the same species as all expression files."
        ]
        error_msg = [
            "Selected genes are 'Mus musculus', while expression 'expression' is 'Homo sapiens'."
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Select at least two samples to compute hierarchical clustering of samples."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_2.pk,
            ],
            "preprocessing": {
                "genes": [
                    "gene_1",
                ],
                "source": "ENSEMBL",
                "species": "Homo sapiens",
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            (
                "Select at least two genes to compute hierarchical clustering of samples "
                "with correlation distance metric or use Euclidean distance metric."
            )
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_2.pk,
            ],
            "preprocessing": {
                "genes": [
                    "gene_1",
                    "gene_2",
                ],
                "source": "ENSEMBL",
                "species": "Homo sapiens",
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        error_msg = ["None of the selected genes are present in all samples."]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_3.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        error_msg = ["The selected samples do not have any common genes."]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_2.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            (
                "The selected samples contain only one common gene (DEFB125). At "
                "least two common genes are required to compute hierarchical clustering of "
                "samples with correlation distance metric. Select a different set of samples "
                "or use Euclidean distance metric."
            )
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_2.pk,
            ],
            "preprocessing": {
                "genes": [
                    "ENSG00000178591",
                    "gene_2",
                ],
                "source": "ENSEMBL",
                "species": "Homo sapiens",
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            (
                "Only one of the selected genes (DEFB125) is present in all samples "
                "but at least two such genes are required to compute hierarchical clustering "
                "of samples with correlation distance metric. Select more genes or use "
                "Euclidean distance metric."
            )
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_6.pk,
                expression_6.pk,
            ],
            "preprocessing": {
                "z_score": False,
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            (
                "All of the selected samples have constant expression across genes. "
                "Hierarchical clustering of samples cannot be computed."
            )
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_7.pk,
            ],
        }
        clustering = self.run_process("clustering-hierarchical-samples", inputs)

        inputs = {
            "exps": [
                expression_3.pk,
                expression_6.pk,
            ],
            "preprocessing": {
                "z_score": False,
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-samples", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Only one of the selected samples (expression 3) has a non-constant expression "
            "across genes. However, hierarchical clustering of samples cannot be computed "
            "with just one sample."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_3.pk,
                expression_5.pk,
                expression_6.pk,
                expression_6.pk,
            ],
            "preprocessing": {
                "z_score": False,
            },
        }
        clustering = self.run_process("clustering-hierarchical-samples", inputs)
        warning_msg = [
            (
                "2 of the selected samples (expression 6, expression 6) have constant "
                "expression across genes. Those samples are excluded from the computation of "
                "hierarchical clustering of samples with correlation distance metric."
            ),
            (
                "3 genes (C20orf96, ZCCHC3, NRSN2-AS1) are present in some but not all of the "
                "selected samples. Those genes are excluded from the computation of hierarchical "
                "clustering of samples."
            ),
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        saved_json, test_json = self.get_json(
            "sample_cluster_data.json.gz", clustering.output["cluster"]
        )
        self.assertEqual(test_json["linkage"], saved_json["linkage"])

    @with_resolwe_host
    @tag_process("clustering-hierarchical-genes")
    def test_gene_clustering(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(
                f_exp="clustering_ensembl_1.tab.gz",
                f_type="TPM",
                name="expression",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            wrong_source = self.prepare_expression(
                f_exp="clustering_ensembl_1.tab.gz",
                f_type="TPM",
                name="wrong source",
                source="UCSC",
                species="Homo sapiens",
            )
            wrong_species = self.prepare_expression(
                f_exp="clustering_ensembl_1.tab.gz",
                f_type="TPM",
                name="wrong species",
                source="ENSEMBL",
                species="Mus musculus",
            )
            wrong_expression_type = self.prepare_expression(
                f_exp="clustering_ensembl_1.tab.gz",
                f_type="FPKM",
                name="wrong expression type",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            wrong_feature_type = self.prepare_expression(
                f_exp="clustering_ensembl_1.tab.gz",
                f_type="TPM",
                feature_type="transcript",
                name="wrong feature type",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_2 = self.prepare_expression(
                f_exp="clustering_ensembl_2.tab.gz",
                f_type="TPM",
                name="expression 2",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_3 = self.prepare_expression(
                f_exp="clustering_ensembl_3.tab.gz",
                f_type="TPM",
                name="expression 3",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_4 = self.prepare_expression(
                f_exp="clustering_ensembl_4.tab.gz",
                f_type="TPM",
                name="expression 4",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_5 = self.prepare_expression(
                f_exp="clustering_ensembl_5.tab.gz",
                f_type="TPM",
                name="expression 5",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_7 = self.prepare_expression(
                f_exp="clustering_ensembl_7.tab.gz",
                f_type="TPM",
                name="expression 7",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_8 = self.prepare_expression(
                f_exp="clustering_int1.tab.gz",
                f_type="TPM",
                name="expression int 1",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_9 = self.prepare_expression(
                f_exp="clustering_int2.tab.gz",
                f_type="TPM",
                name="expression int 2",
                source="ENSEMBL",
                species="Homo sapiens",
            )

        inputs = {
            "exps": [
                expression_8.pk,
                expression_9.pk,
            ],
            "preprocessing": {
                "log2": False,
            },
        }
        clustering = self.run_process("clustering-hierarchical-genes", inputs)
        saved_json, test_json = self.get_json(
            "clustering_out_genes.json.gz", clustering.output["cluster"]
        )
        self.assertEqual(saved_json.keys(), test_json.keys())
        for key in saved_json:
            assert_method = (
                self.assertAlmostEqualGeneric if key == "linkage" else self.assertEqual
            )
            assert_method(saved_json[key], test_json[key])
        self.assertEqual(len(clustering.process_warning), 0)

        inputs = {
            "exps": [
                expression_1.pk,
                wrong_source.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        warning_msg = [
            "All expression data must be annotated by the same genome database."
        ]
        error_msg = [
            (
                "Sample 'expression' has 'ENSEMBL' gene IDs, "
                "while sample 'wrong source' has 'UCSC' gene IDs."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                wrong_species.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        warning_msg = ["All expressions must be of the same Species."]
        error_msg = [
            (
                "Sample 'expression' is 'Homo sapiens', "
                "while sample 'wrong species' is 'Mus musculus'."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                wrong_expression_type.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        warning_msg = ["All expressions must be of the same Expression type."]
        error_msg = [
            (
                "Expression 'expression' has 'TPM' expression type, "
                "while sample 'wrong expression type' has 'FPKM' expression type."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                wrong_feature_type.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        warning_msg = ["All expressions must be of the same Feature type."]
        error_msg = [
            (
                "Expression 'expression' has 'gene' feature type, "
                "while sample 'wrong feature type' has 'transcript' feature type."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
            ],
            "preprocessing": {
                "genes": ["gene"],
                "source": "UCSC",
                "species": "Homo sapiens",
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        warning_msg = [
            "Selected genes must be annotated by the same genome database as all expression files."
        ]
        error_msg = [
            (
                "Gene IDs are from 'UCSC' database, "
                "while sample 'expression' has gene IDs from 'ENSEMBL' database."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
            ],
            "preprocessing": {
                "genes": ["gene"],
                "source": "ENSEMBL",
                "species": "Mus musculus",
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        warning_msg = [
            "Selected genes must be from the same species as all expression files."
        ]
        error_msg = [
            "Selected genes are 'Mus musculus', while expression 'expression' is 'Homo sapiens'."
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_2.pk,
            ],
            "preprocessing": {
                "genes": [
                    "gene label",
                ],
                "species": "Homo sapiens",
                "source": "ENSEMBL",
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Select at least two genes to compute hierarchical clustering of genes."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            (
                "Select at least two samples to compute hierarchical clustering of genes with "
                "correlation distance metric or use Euclidean distance metric."
            ),
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_3.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        error_msg = ["The selected samples do not have any common genes."]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_3.pk,
            ],
            "preprocessing": {
                "genes": [
                    "gene_1",
                    "gene_2",
                ],
                "source": "ENSEMBL",
                "species": "Homo sapiens",
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        error_msg = ["None of the selected genes are present in all samples."]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_2.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            (
                "The selected samples contain only one common gene (DEFB125). At "
                "least two common genes are required to compute hierarchical clustering of "
                "genes with correlation distance metric. Select a different set of samples "
                "or use Euclidean distance metric."
            )
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_2.pk,
            ],
            "preprocessing": {
                "genes": [
                    "ENSG00000178591",
                    "gene_2",
                ],
                "source": "ENSEMBL",
                "species": "Homo sapiens",
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            (
                "Only one of the selected genes (DEFB125) is present in all samples "
                "but at least two such genes are required to compute hierarchical clustering "
                "of genes with correlation distance metric. Select more genes or use "
                "Euclidean distance metric."
            )
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_7.pk,
            ],
            "preprocessing": {
                "genes": [
                    "ENSG00000185982",
                    "ENSG00000125903",
                    "ENSG00000186458",  # this gene is missing in both expression files
                ],
                "source": "ENSEMBL",
                "species": "Homo sapiens",
            },
        }
        clustering = self.run_process("clustering-hierarchical-genes", inputs)
        warning_msg = [
            (
                "1 of the selected genes (DEFB132) is missing in at least one of the "
                "selected samples. This gene is excluded from the computation of "
                "hierarchical clustering of genes."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        saved_json, test_json = self.get_json(
            "gene_cluster_filtered.json.gz", clustering.output["cluster"]
        )
        self.assertEqual(test_json["linkage"], saved_json["linkage"])

        inputs = {
            "exps": [
                expression_1.pk,
                expression_1.pk,
            ],
            "preprocessing": {
                "z_score": False,
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            (
                "All of the selected genes have constant expression across samples. "
                "Hierarchical clustering of genes cannot be computed."
            )
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_3.pk,
                expression_4.pk,
            ],
            "preprocessing": {
                "z_score": False,
            },
        }
        clustering = self.run_process(
            "clustering-hierarchical-genes", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Only one of the selected genes (DEFB132) has a non-constant expression "
            "across samples. However, hierarchical clustering of genes cannot be computed "
            "with just one gene."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "exps": [
                expression_1.pk,
                expression_7.pk,
            ],
        }
        clustering = self.run_process("clustering-hierarchical-genes", inputs)

        inputs = {
            "exps": [
                expression_3.pk,
                expression_5.pk,
            ],
            "preprocessing": {
                "z_score": False,
            },
        }
        clustering = self.run_process("clustering-hierarchical-genes", inputs)
        warning_msg = [
            (
                "3 genes (C20orf96, NRSN2-AS1, ZCCHC3) are present in "
                "some but not all of the selected samples. Those genes are excluded from the "
                "computation of hierarchical clustering of genes."
            )
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        saved_json, test_json = self.get_json(
            "gene_cluster_data.json.gz", clustering.output["cluster"]
        )
        self.assertEqual(test_json["linkage"], saved_json["linkage"])

    fixtures = ["relationtypes.yaml"]

    @with_resolwe_host
    @tag_process("clustering-hierarchical-etc")
    def test_etc_clustering(self):
        base = Path("test_etc_clustering")
        inputs = base / "inputs"
        outputs = base / "outputs"
        with self.preparation_stage():
            exp_1 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_4h.txt.gz"),
                name="expression",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_2 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_4h_2.txt.gz"),
                name="expression 2",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_3 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_8h.txt.gz"),
                name="expression 3",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_4 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_8h_2.txt.gz"),
                name="expression 4",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_5 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_12h.txt.gz"),
                name="expression 5",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_6 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_12h_2.txt.gz"),
                name="expression 6",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            wrong_source = self.prepare_expression(
                f_exp=str(inputs / "discoideum_4h.txt.gz"),
                name="wrong source",
                source="UCSC",
                species="Homo sapiens",
            )
            wrong_species = self.prepare_expression(
                f_exp=str(inputs / "discoideum_4h.txt.gz"),
                name="wrong species",
                source="dictybase",
                species="Dictyostelium purpureum",
            )
            wrong_expression_type = self.prepare_expression(
                f_exp=str(inputs / "discoideum_4h.txt.gz"),
                f_type="FPKM",
                name="wrong expression type",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            wrong_feature_type = self.prepare_expression(
                f_exp=str(inputs / "discoideum_4h.txt.gz"),
                feature_type="transcript",
                name="wrong feature type",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_constant = self.prepare_expression(
                f_exp=str(inputs / "discoideum_constant_8h.txt.gz"),
                name="expression constant",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            different_genes = self.prepare_expression(
                f_exp=str(inputs / "discoideum_different_genes.txt.gz"),
                name="different genes",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_1_copy = self.prepare_expression(
                f_exp=str(inputs / "discoideum_4h.txt.gz"),
                name="same expression",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_2_part = self.prepare_expression(
                f_exp=str(inputs / "discoideum_4h_2_partial_match.txt.gz"),
                name="Partially matching gene names",
                source="dictybase",
                species="Dictyostelium discoideum",
            )

            rel_type_series = RelationType.objects.get(name="series")
            relation = Relation.objects.create(
                contributor=self.contributor,
                collection=self.collection,
                type=rel_type_series,
                category="time-series",
                unit=Relation.UNIT_HOUR,
            )
            assign_perm("view_relation", self.contributor, relation)

            RelationPartition.objects.create(
                relation=relation, entity=exp_1.entity, label="4h", position=4
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_2.entity, label="4h", position=4
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_3.entity, label="8h", position=8
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_4.entity, label="8h", position=8
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_5.entity, label="12h", position=12
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_6.entity, label="12h", position=12
            )
            RelationPartition.objects.create(
                relation=relation, entity=wrong_source.entity, label="4h", position=4
            )
            RelationPartition.objects.create(
                relation=relation, entity=wrong_species.entity, label="4h", position=4
            )
            RelationPartition.objects.create(
                relation=relation,
                entity=wrong_expression_type.entity,
                label="4h",
                position=4,
            )
            RelationPartition.objects.create(
                relation=relation,
                entity=wrong_feature_type.entity,
                label="4h",
                position=4,
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_constant.entity, label="8h", position=8
            )
            RelationPartition.objects.create(
                relation=relation, entity=different_genes.entity, label="8h", position=8
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_1_copy.entity, label="8h", position=8
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_2_part.entity, label="4h", position=4
            )

        inputs = {
            "expressions": [
                exp_1.pk,
                exp_2.pk,
                exp_3.pk,
                exp_4.pk,
                exp_5.pk,
                exp_6.pk,
            ],
        }
        clustering = self.run_process("clustering-hierarchical-etc", inputs)
        self.assertJSON(
            clustering, clustering.output["cluster"], "", outputs / "clustering.json.gz"
        )
        self.assertEqual(len(clustering.process_warning), 0)

        inputs = {
            "expressions": [
                exp_1.pk,
                exp_2.pk,
                exp_3.pk,
                exp_4.pk,
                exp_5.pk,
                exp_6.pk,
            ],
            "genes": [
                "DDB_G0271520",
                "DDB_G0271524",
                "DDB_G0279413",
                "DDB_G0283907",
                "DDB_G0290157",
                "DDB_G0293524",
            ],
            "gene_species": "Dictyostelium discoideum",
            "gene_source": "dictybase",
            "distance": "pearson",
            "linkage": "complete",
        }
        clustering = self.run_process("clustering-hierarchical-etc", inputs)
        self.assertJSON(
            clustering,
            clustering.output["cluster"],
            "",
            outputs / "clustering_genes.json.gz",
        )
        self.assertEqual(len(clustering.process_warning), 0)

        inputs = {
            "expressions": [
                exp_1.pk,
                wrong_source.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-etc", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Input samples are of different Gene ID databases: UCSC and dictybase."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "expressions": [
                exp_1.pk,
                wrong_species.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-etc", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Input samples are of different Species: "
            "Dictyostelium purpureum and Dictyostelium discoideum."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "expressions": [
                exp_1.pk,
                wrong_expression_type.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-etc", inputs, Data.STATUS_ERROR
        )
        error_msg = ["Input samples are of different Expression types: FPKM and TPM."]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "expressions": [
                exp_1.pk,
                wrong_feature_type.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-etc", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Input samples are of different Feature type: transcript and gene."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "expressions": [
                exp_1.pk,
                exp_3.pk,
            ],
            "genes": [
                "DDB_G0271520",
            ],
            "gene_species": "Dictyostelium discoideum",
            "gene_source": "dictybase",
        }
        clustering = self.run_process(
            "clustering-hierarchical-etc", inputs, Data.STATUS_ERROR
        )
        error_msg = ["At least two genes have to be selected."]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "expressions": [
                exp_1.pk,
                exp_3.pk,
            ],
            "genes": [
                "DDB_G0271520",
                "DDB_G0271524",
            ],
            "gene_species": "Wrong",
            "gene_source": "dictybase",
        }
        clustering = self.run_process(
            "clustering-hierarchical-etc", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Selected genes must belong to the same species as "
            "expression files. Instead genes belong to Wrong while "
            "expressions belong to Dictyostelium discoideum."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "expressions": [
                exp_1.pk,
                exp_3.pk,
            ],
            "genes": [
                "DDB_G0271520",
                "DDB_G0271524",
            ],
            "gene_species": "Dictyostelium discoideum",
            "gene_source": "wrong-source",
        }
        clustering = self.run_process(
            "clustering-hierarchical-etc", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Selected genes must be annotated by the same genome "
            "database as expressions. Instead Gene IDs of genes are "
            "from wrong-source and expressions have IDs from dictybase."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "expressions": [
                exp_1.pk,
                exp_constant.pk,
            ],
        }
        clustering = self.run_process("clustering-hierarchical-etc", inputs)
        warn_msg = [
            "Genes not present in all of the selected samples are "
            "excluded from the analysis. Excluded 1 "
            "of them (DDB_G0293524).",
            "3 genes (DDB_G0271520, DDB_G0271524, DDB_G0279413) have "
            "constant expression across time points. Those genes are "
            "excluded from the computation of hierarchical clustering "
            "of genes.",
        ]
        self.assertEqual(clustering.process_warning, warn_msg)

        inputs = {
            "expressions": [
                exp_1.pk,
                exp_1_copy.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-etc", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "All genes have constant expression across time points. "
            "Hierarchical clustering of genes cannot be computed."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "expressions": [exp_1.pk, exp_2.pk],
        }
        clustering = self.run_process(
            "clustering-hierarchical-etc", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Only one time point was provided. At least two time "
            "points are required to run hierarhical clustering."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "expressions": [
                exp_1.pk,
                exp_3.pk,
            ],
            "genes": [
                "DDB_G0000000",
                "DDB_G1111111",
            ],
            "gene_species": "Dictyostelium discoideum",
            "gene_source": "dictybase",
        }
        clustering = self.run_process(
            "clustering-hierarchical-etc", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "At least two of the selected genes have to be present in "
            "all samples to run hierarhical clustering. 0 found in all "
            "samples."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "expressions": [
                exp_1.pk,
                different_genes.pk,
            ],
        }
        clustering = self.run_process(
            "clustering-hierarchical-etc", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "At least two genes shared across all samples are required "
            "to run hierarhical clustering. 0 found in all samples."
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            "expressions": [
                exp_1.pk,
                exp_2_part.pk,
                exp_3.pk,
                exp_4.pk,
            ],
        }
        clustering = self.run_process("clustering-hierarchical-etc", inputs)

    @with_resolwe_host
    @tag_process("find-similar")
    def test_find_similar(self):
        inputs = Path("test_etc_clustering") / "inputs"
        outputs = Path("find_similar") / "outputs"
        with self.preparation_stage():
            exp_1 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_4h.txt.gz"),
                name="expression",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_2 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_4h_2.txt.gz"),
                name="expression 2",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_3 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_8h.txt.gz"),
                name="expression 3",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_4 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_8h_2.txt.gz"),
                name="expression 4",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_5 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_12h.txt.gz"),
                name="expression 5",
                source="dictybase",
                species="Dictyostelium discoideum",
            )
            exp_6 = self.prepare_expression(
                f_exp=str(inputs / "discoideum_12h_2.txt.gz"),
                name="expression 6",
                source="dictybase",
                species="Dictyostelium discoideum",
            )

            rel_type_series = RelationType.objects.get(name="series")
            relation = Relation.objects.create(
                contributor=self.contributor,
                collection=self.collection,
                type=rel_type_series,
                category="time-series",
                unit=Relation.UNIT_HOUR,
            )
            assign_perm("view_relation", self.contributor, relation)

            RelationPartition.objects.create(
                relation=relation, entity=exp_1.entity, label="4h", position=4
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_2.entity, label="4h", position=4
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_3.entity, label="8h", position=8
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_4.entity, label="8h", position=8
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_5.entity, label="12h", position=12
            )
            RelationPartition.objects.create(
                relation=relation, entity=exp_6.entity, label="12h", position=12
            )
        inputs = {
            "expressions": [
                exp_1.pk,
                exp_2.pk,
                exp_3.pk,
                exp_4.pk,
                exp_5.pk,
                exp_6.pk,
            ],
            "gene": "DDB_G0283907",
        }
        similar = self.run_process("find-similar", inputs)
        self.assertJSON(
            similar,
            similar.output["similar_genes"],
            "",
            outputs / "similar_genes.json.gz",
        )

        inputs = {
            "expressions": [
                exp_1.pk,
                exp_2.pk,
                exp_3.pk,
                exp_4.pk,
                exp_5.pk,
                exp_6.pk,
            ],
            "gene": "DDB_G0283907",
            "distance": "pearson",
        }
        similar = self.run_process("find-similar", inputs)
        self.assertJSON(
            similar,
            similar.output["similar_genes"],
            "",
            outputs / "similar_genes_pearson.json.gz",
        )

        inputs = {
            "expressions": [
                exp_1.pk,
                exp_2.pk,
                exp_3.pk,
                exp_4.pk,
                exp_5.pk,
                exp_6.pk,
            ],
            "gene": "WRONG",
            "distance": "pearson",
        }
        similar = self.run_process("find-similar", inputs, Data.STATUS_ERROR)
        error_msg = [
            "Selected query gene was not found. Please make sure the "
            "selected gene name can be found in all expression time "
            "courses."
        ]
        self.assertEqual(similar.process_error, error_msg)
