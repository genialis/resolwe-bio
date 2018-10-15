# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process
from resolwe_bio.utils.test import with_resolwe_host, KBBioProcessTestCase


class ClusteringProcessTestCase(KBBioProcessTestCase):

    @with_resolwe_host
    @tag_process('clustering-hierarchical-samples')
    def test_sample_clustering(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(
                f_exp='clustering_ensembl_1.tab.gz',
                f_type='TPM',
                name='expression',
                source='ENSEMBL',
                species='Homo sapiens')
            wrong_source = self.prepare_expression(
                f_exp='clustering_ensembl_1.tab.gz',
                f_type='TPM',
                name='wrong source',
                source='UCSC',
                species='Homo sapiens')
            wrong_species = self.prepare_expression(
                f_exp='clustering_ensembl_1.tab.gz',
                f_type='TPM',
                name='wrong species',
                source='ENSEMBL',
                species='Mus musculus')
            wrong_expression_type = self.prepare_expression(
                f_exp='clustering_ensembl_1.tab.gz',
                f_type='FPKM',
                name='wrong expression type',
                source='ENSEMBL',
                species='Homo sapiens')
            wrong_feature_type = self.prepare_expression(
                f_exp='clustering_ensembl_1.tab.gz',
                f_type='TPM',
                feature_type='transcript',
                name='wrong feature type',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_2 = self.prepare_expression(
                f_exp='clustering_ensembl_2.tab.gz',
                f_type='TPM',
                name='expression 2',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_3 = self.prepare_expression(
                f_exp='clustering_ensembl_3.tab.gz',
                f_type='TPM',
                name='expression 3',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_5 = self.prepare_expression(
                f_exp='clustering_ensembl_5.tab.gz',
                f_type='TPM',
                name='expression 5',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_6 = self.prepare_expression(
                f_exp='clustering_ensembl_6.tab.gz',
                f_type='TPM',
                name='expression 6',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_7 = self.prepare_expression(
                f_exp='clustering_ensembl_7.tab.gz',
                f_type='TPM',
                name='expression 7',
                source='ENSEMBL',
                species='Homo sapiens')

        inputs = {
            'exps': [
                expression_1.pk,
                wrong_source.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        warning_msg = ['All expression data must be annotated by the same genome database.']
        error_msg = [("Sample 'expression' has 'ENSEMBL' gene IDs, "
                      "while sample 'wrong source' has 'UCSC' gene IDs.")]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                wrong_species.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        warning_msg = ['All expressions must be of the same Species.']
        error_msg = [("Sample 'expression' is 'Homo sapiens', "
                      "while sample 'wrong species' is 'Mus musculus'.")]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                wrong_expression_type.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        warning_msg = ['All expressions must be of the same Expression type.']
        error_msg = [("Expression 'expression' has 'TPM' expression type, "
                      "while sample 'wrong expression type' has 'FPKM' expression type.")]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                wrong_feature_type.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        warning_msg = ['All expressions must be of the same Feature type.']
        error_msg = [("Expression 'expression' has 'gene' feature type, "
                      "while sample 'wrong feature type' has 'transcript' feature type.")]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
            ],
            'preprocessing': {
                'genes': ['gene'],
                'source': 'UCSC',
                'species': 'Homo sapiens',
            },
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        warning_msg = ['Selected genes must be annotated by the same genome database as all expression files.']
        error_msg = [("Gene IDs are from 'UCSC' database, "
                      "while sample 'expression' has gene IDs from 'ENSEMBL' database.")]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
            ],
            'preprocessing': {
                'genes': ['gene'],
                'source': 'ENSEMBL',
                'species': 'Mus musculus',
            },
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        warning_msg = ['Selected genes must be from the same species as all expression files.']
        error_msg = ["Selected genes are 'Mus musculus', while expression 'expression' is 'Homo sapiens'."]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        error_msg = ['Select at least two samples to compute hierarchical clustering of samples.']
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_2.pk,
            ],
            'preprocessing': {
                'genes': [
                    'gene_1',
                ],
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
            },
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        error_msg = [('Select at least two genes to compute hierarchical clustering of samples '
                      'with correlation distance metric or use Euclidean distance metric.')]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_2.pk,
            ],
            'preprocessing': {
                'genes': [
                    'gene_1',
                    'gene_2',
                ],
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
            },
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        error_msg = ['None of the selected genes are present in all samples.']
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_3.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        error_msg = ['The selected samples do not have any common genes.']
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_2.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        error_msg = [('The selected samples contain only one common gene (DEFB125). At '
                      'least two common genes are required to compute hierarchical clustering of '
                      'samples with correlation distance metric. Select a different set of samples '
                      'or use Euclidean distance metric.')]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_2.pk,
            ],
            'preprocessing': {
                'genes': [
                    'ENSG00000178591',
                    'gene_2',
                ],
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
            },
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        error_msg = [('Only one of the selected genes (DEFB125) is present in all samples '
                      'but at least two such genes are required to compute hierarchical clustering '
                      'of samples with correlation distance metric. Select more genes or use '
                      'Euclidean distance metric.')]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_6.pk,
                expression_6.pk,
            ],
            'preprocessing': {
                'z_score': False,
            },
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        error_msg = [('All of the selected samples have constant expression across genes. '
                      'Hierarchical clustering of samples cannot be computed.')]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_7.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs)

        inputs = {
            'exps': [
                expression_3.pk,
                expression_6.pk,
            ],
            'preprocessing': {
                'z_score': False,
            },
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs, Data.STATUS_ERROR)
        error_msg = ['Only one of the selected samples (expression 3) has a non-constant expression '
                     'across genes. However, hierarchical clustering of samples cannot be computed '
                     'with just one sample.']
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_3.pk,
                expression_5.pk,
                expression_6.pk,
                expression_6.pk,
            ],
            'preprocessing': {
                'z_score': False,
            },
        }
        clustering = self.run_process('clustering-hierarchical-samples', inputs)
        warning_msg = [
            ('2 of the selected samples (expression 6, expression 6) have constant '
             'expression across genes. Those samples are excluded from the computation of '
             'hierarchical clustering of samples with correlation distance metric.'),
            ('3 genes (C20orf96, NRSN2-AS1, ZCCHC3) are present in some but not all of the '
             'selected samples. Those genes are excluded from the computation of hierarchical '
             'clustering of samples.'),
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        saved_json, test_json = self.get_json('sample_cluster_data.json.gz', clustering.output['cluster'])
        self.assertEqual(test_json['linkage'], saved_json['linkage'])

    @with_resolwe_host
    @tag_process('clustering-hierarchical-genes')
    def test_gene_clustering(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(
                f_exp='clustering_ensembl_1.tab.gz',
                f_type='TPM',
                name='expression',
                source='ENSEMBL',
                species='Homo sapiens')
            wrong_source = self.prepare_expression(
                f_exp='clustering_ensembl_1.tab.gz',
                f_type='TPM',
                name='wrong source',
                source='UCSC',
                species='Homo sapiens')
            wrong_species = self.prepare_expression(
                f_exp='clustering_ensembl_1.tab.gz',
                f_type='TPM',
                name='wrong species',
                source='ENSEMBL',
                species='Mus musculus')
            wrong_expression_type = self.prepare_expression(
                f_exp='clustering_ensembl_1.tab.gz',
                f_type='FPKM',
                name='wrong expression type',
                source='ENSEMBL',
                species='Homo sapiens')
            wrong_feature_type = self.prepare_expression(
                f_exp='clustering_ensembl_1.tab.gz',
                f_type='TPM',
                feature_type='transcript',
                name='wrong feature type',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_2 = self.prepare_expression(
                f_exp='clustering_ensembl_2.tab.gz',
                f_type='TPM',
                name='expression 2',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_3 = self.prepare_expression(
                f_exp='clustering_ensembl_3.tab.gz',
                f_type='TPM',
                name='expression 3',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_4 = self.prepare_expression(
                f_exp='clustering_ensembl_4.tab.gz',
                f_type='TPM',
                name='expression 4',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_5 = self.prepare_expression(
                f_exp='clustering_ensembl_5.tab.gz',
                f_type='TPM',
                name='expression 5',
                source='ENSEMBL',
                species='Homo sapiens')
            expression_7 = self.prepare_expression(
                f_exp='clustering_ensembl_7.tab.gz',
                f_type='TPM',
                name='expression 7',
                source='ENSEMBL',
                species='Homo sapiens')

        inputs = {
            'exps': [
                expression_1.pk,
                wrong_source.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        warning_msg = ['All expression data must be annotated by the same genome database.']
        error_msg = [("Sample 'expression' has 'ENSEMBL' gene IDs, "
                      "while sample 'wrong source' has 'UCSC' gene IDs.")]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                wrong_species.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        warning_msg = ['All expressions must be of the same Species.']
        error_msg = [("Sample 'expression' is 'Homo sapiens', "
                      "while sample 'wrong species' is 'Mus musculus'.")]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                wrong_expression_type.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        warning_msg = ['All expressions must be of the same Expression type.']
        error_msg = [("Expression 'expression' has 'TPM' expression type, "
                      "while sample 'wrong expression type' has 'FPKM' expression type.")]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                wrong_feature_type.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        warning_msg = ['All expressions must be of the same Feature type.']
        error_msg = [("Expression 'expression' has 'gene' feature type, "
                      "while sample 'wrong feature type' has 'transcript' feature type.")]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
            ],
            'preprocessing': {
                'genes': ['gene'],
                'source': 'UCSC',
                'species': 'Homo sapiens',
            },
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        warning_msg = ['Selected genes must be annotated by the same genome database as all expression files.']
        error_msg = [("Gene IDs are from 'UCSC' database, "
                      "while sample 'expression' has gene IDs from 'ENSEMBL' database.")]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
            ],
            'preprocessing': {
                'genes': ['gene'],
                'source': 'ENSEMBL',
                'species': 'Mus musculus',
            },
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        warning_msg = ['Selected genes must be from the same species as all expression files.']
        error_msg = ["Selected genes are 'Mus musculus', while expression 'expression' is 'Homo sapiens'."]
        self.assertEqual(clustering.process_warning, warning_msg)
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_2.pk,
            ],
            'preprocessing': {
                'genes': [
                    'gene label',
                ],
                'species': 'Homo sapiens',
                'source': 'ENSEMBL',
            },
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        error_msg = ['Select at least two genes to compute hierarchical clustering of genes.']
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        error_msg = [
            ('Select at least two samples to compute hierarchical clustering of genes with '
             'correlation distance metric or use Euclidean distance metric.'),
        ]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_3.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        error_msg = ['The selected samples do not have any common genes.']
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_3.pk,
            ],
            'preprocessing': {
                'genes': [
                    'gene_1',
                    'gene_2',
                ],
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
            },
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        error_msg = ['None of the selected genes are present in all samples.']
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_2.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        error_msg = [('The selected samples contain only one common gene (DEFB125). At '
                      'least two common genes are required to compute hierarchical clustering of '
                      'genes with correlation distance metric. Select a different set of samples '
                      'or use Euclidean distance metric.')]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_2.pk,
            ],
            'preprocessing': {
                'genes': [
                    'ENSG00000178591',
                    'gene_2',
                ],
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
            },
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        error_msg = [('Only one of the selected genes (DEFB125) is present in all samples '
                      'but at least two such genes are required to compute hierarchical clustering '
                      'of genes with correlation distance metric. Select more genes or use '
                      'Euclidean distance metric.')]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_7.pk,
            ],
            'preprocessing': {
                'genes': [
                    'ENSG00000185982',
                    'ENSG00000125903',
                    'ENSG00000186458',  # this gene is missing in both expression files
                ],
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
            },
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs)
        warning_msg = [('1 of the selected genes (DEFB132) is missing in at least one of the '
                        'selected samples. This gene is excluded from the computation of '
                        'hierarchical clustering of genes.')]
        self.assertEqual(clustering.process_warning, warning_msg)
        saved_json, test_json = self.get_json('gene_cluster_filtered.json.gz', clustering.output['cluster'])
        self.assertEqual(test_json['linkage'], saved_json['linkage'])

        inputs = {
            'exps': [
                expression_1.pk,
                expression_1.pk,
            ],
            'preprocessing': {
                'z_score': False,
            }
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        error_msg = [('All of the selected genes have constant expression across samples. '
                      'Hierarchical clustering of genes cannot be computed.')]
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_3.pk,
                expression_4.pk,
            ],
            'preprocessing': {
                'z_score': False,
            },
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs, Data.STATUS_ERROR)
        error_msg = ['Only one of the selected genes (DEFB132) has a non-constant expression '
                     'across samples. However, hierarchical clustering of genes cannot be computed '
                     'with just one gene.']
        self.assertEqual(clustering.process_error, error_msg)

        inputs = {
            'exps': [
                expression_1.pk,
                expression_7.pk,
            ],
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs)

        inputs = {
            'exps': [
                expression_3.pk,
                expression_5.pk,
            ],
            'preprocessing': {
                'z_score': False,
            },
        }
        clustering = self.run_process('clustering-hierarchical-genes', inputs)
        warning_msg = [
            ('3 genes (C20orf96, NRSN2-AS1, ZCCHC3) are present in '
             'some but not all of the selected samples. Those genes are excluded from the '
             'computation of hierarchical clustering of genes.')
        ]
        self.assertEqual(clustering.process_warning, warning_msg)
        saved_json, test_json = self.get_json('gene_cluster_data.json.gz', clustering.output['cluster'])
        self.assertEqual(test_json['linkage'], saved_json['linkage'])
