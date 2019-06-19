# pylint: disable=missing-docstring
import os

from resolwe.flow.models import Data, Collection, Relation
from resolwe.flow.models.entity import RelationPartition, RelationType
from resolwe.test import tag_process

from resolwe_bio.expression_filters.relation import replicate_groups

from resolwe_bio.utils.test import with_resolwe_host, KBBioProcessTestCase


class ExpressionProcessorTestCase(KBBioProcessTestCase):

    fixtures = ['relationtypes.yaml']

    @tag_process('cufflinks', 'cuffmerge')
    def test_cufflinks(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            annotation_gtf = self.prepare_annotation('annotation dicty.gff.gz')
            annotation_gff3 = self.prepare_annotation_gff()

            aligned_reads = self.run_process('alignment-hisat2', {
                'genome': genome.pk,
                'reads': reads.pk,
                'spliced_alignments': {
                    'cufflinks': True
                }
            })

        inputs = {
            'alignment': aligned_reads.pk,
            'annotation': annotation_gtf.pk,
            'genome': genome.pk}
        cuff_exp = self.run_process('cufflinks', inputs)
        self.assertFile(cuff_exp, 'transcripts', 'cufflinks_transcripts.gtf', sort=True)
        self.assertFields(cuff_exp, 'species', 'Dictyostelium discoideum')
        self.assertFields(cuff_exp, 'build', 'dd-05-2009')
        self.assertFields(cuff_exp, 'source', 'DICTYBASE')

        inputs = {
            'alignment': aligned_reads.pk,
            'annotation': annotation_gtf.pk,
            'genome': genome.pk}
        cuff_exp2 = self.run_process('cufflinks', inputs)

        inputs = {
            'expressions': [cuff_exp.pk, cuff_exp2.pk],
            'gff': annotation_gff3.pk,
            'genome': genome.pk}
        cuff_merge_gff3 = self.run_process('cuffmerge', inputs)
        self.assertFile(cuff_merge_gff3, 'annot', 'cuffmerge_transcripts.gtf')
        self.assertFields(cuff_merge_gff3, 'species', 'Dictyostelium discoideum')
        self.assertFields(cuff_merge_gff3, 'build', 'dd-05-2009')
        self.assertFields(cuff_exp, 'source', 'DICTYBASE')

        inputs['gff'] = annotation_gtf.pk
        cuff_merge_gtf = self.run_process('cuffmerge', inputs)
        self.assertFile(cuff_merge_gtf, 'annot', 'cuffmerge_transcripts.gtf')
        self.assertFields(cuff_merge_gtf, 'species', 'Dictyostelium discoideum')
        self.assertFields(cuff_merge_gtf, 'build', 'dd-05-2009')
        self.assertFields(cuff_exp, 'source', 'DICTYBASE')

    @tag_process('cuffquant')
    def test_cuffquant(self):
        with self.preparation_stage():
            inputs = {
                'src': 'cuffquant_mapping.bam',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            bam = self.run_process('upload-bam', inputs)

            annotation = self.prepare_annotation(
                fn='hg19_chr20_small.gtf.gz',
                source='UCSC',
                species='Homo sapiens',
                build='hg19'
            )

        inputs = {
            'alignment': bam.id,
            'annotation': annotation.id}
        cuffquant = self.run_process('cuffquant', inputs)
        self.assertFields(cuffquant, 'species', 'Homo sapiens')
        self.assertFields(cuffquant, 'build', 'hg19')
        self.assertFields(cuffquant, 'source', 'UCSC')

    @with_resolwe_host
    @tag_process('cuffnorm')
    def test_cuffnorm(self):
        with self.preparation_stage():
            collection = Collection.objects.create(
                name='Test collection',
                contributor=self.contributor
            )

            rel_type_group = RelationType.objects.get(name='group')

            replicate_group = Relation.objects.create(
                contributor=self.contributor,
                collection=collection,
                type=rel_type_group,
                category='Replicate'
            )

            inputs = {
                'src': 'cuffquant 1.cxb',
                'source': 'UCSC',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            sample_1 = self.run_process("upload-cxb", inputs)

            inputs = {
                'src': 'cuffquant_2.cxb',
                'source': 'UCSC',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            sample_2 = self.run_process("upload-cxb", inputs)

            inputs = {
                'src': '3-cuffquant.cxb',
                'source': 'UCSC',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            sample_3 = self.run_process("upload-cxb", inputs)

            inputs = {
                'src': '4-cuffquant.cxb',
                'source': 'UCSC',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            sample_4 = self.run_process("upload-cxb", inputs)

            inputs = {
                'src': '5-cuffquant.cxb',
                'source': 'UCSC',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            sample_5 = self.run_process("upload-cxb", inputs)

            inputs = {
                'src': '6-cuffquant.cxb',
                'source': 'UCSC',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            sample_6 = self.run_process("upload-cxb", inputs)

            RelationPartition.objects.create(relation=replicate_group, entity=sample_1.entity_set.first(), label='1')
            RelationPartition.objects.create(relation=replicate_group, entity=sample_2.entity_set.first(), label='1')
            RelationPartition.objects.create(relation=replicate_group, entity=sample_3.entity_set.first(), label='2')
            RelationPartition.objects.create(relation=replicate_group, entity=sample_4.entity_set.first(), label='2')
            RelationPartition.objects.create(relation=replicate_group, entity=sample_5.entity_set.first(), label='2')
            RelationPartition.objects.create(relation=replicate_group, entity=sample_6.entity_set.first(), label='3')

            annotation = self.prepare_annotation(fn='hg19_chr20_small.gtf.gz', source='UCSC',
                                                 species='Homo sapiens', build='hg19')

            self.assertEqual(
                replicate_groups([
                    {'__id': sample_1.id},
                    {'__id': sample_2.id},
                    {'__id': sample_3.id},
                    {'__id': sample_4.id},
                    {'__id': sample_5.id},
                    {'__id': sample_6.id}
                ]),
                [1, 1, 2, 2, 2, 3]
            )

        inputs = {
            'cuffquant': [sample_1.pk, sample_2.pk, sample_3.pk, sample_4.pk, sample_5.pk, sample_6.pk],
            'annotation': annotation.id,
        }
        cuffnorm = self.run_process('cuffnorm', inputs)
        self.assertFile(cuffnorm, 'fpkm_means', 'cuffnorm_all_fpkm_means.txt')
        self.assertFile(cuffnorm, 'genes_fpkm', 'cuffnorm_genes.fpkm_table')
        self.assertFileExists(cuffnorm, 'raw_scatter')
        self.assertFields(cuffnorm, 'source', 'UCSC')
        self.assertFields(cuffnorm, 'species', 'Homo sapiens')
        self.assertFields(cuffnorm, 'build', 'hg19')

        exp = Data.objects.last()
        self.assertFile(exp, 'exp', 'cuffnorm_expression.tab.gz', compression='gzip')
        self.assertFile(exp, 'exp_set', 'cuffnorm_out_exp_set.txt.gz', compression='gzip')
        self.assertJSON(exp, exp.output['exp_set_json'], '', 'cuffnorm_exp_set.json.gz')

    @tag_process('mappability-bcm')
    def test_mappability(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            annotation = self.prepare_annotation_gff()

        mappability = self.run_process('mappability-bcm', {
            'genome': genome.id,
            'gff': annotation.id,
            'length': 50,
        })

        self.assertFileExists(mappability, 'mappability')

    @tag_process('expression-dicty', 'etc-bcm')
    def test_expression_dicty(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            annotation = self.prepare_annotation_gff()

            aligned_reads = self.run_process('alignment-hisat2', {
                'genome': genome.pk,
                'reads': reads.pk
            })

            mappa = self.run_process("upload-mappability", {"src": "purpureum_mappability_50.tab.gz"})

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'mappable': mappa.pk}
        expression = self.run_process('expression-dicty', inputs)
        self.assertFile(expression, 'rpkm', 'expression_bcm_rpkm.tab.gz', compression='gzip')
        self.assertFields(expression, "source", "DICTYBASE")
        self.assertFields(expression, 'species', 'Dictyostelium discoideum')
        self.assertFields(expression, 'build', 'dd-05-2009')
        self.assertFields(expression, 'feature_type', 'gene')

        inputs = {'expressions': [expression.pk, expression.pk]}
        etc = self.run_process('etc-bcm', inputs)
        self.assertJSON(etc, etc.output['etc'], '', 'etc.json.gz')

    @with_resolwe_host
    @tag_process('htseq-count')
    def test_expression_htseq(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            inputs = {
                'src': 'annotation dicty.gtf.gz',
                'source': 'DICTYBASE',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009',
            }
            annotation_correct = self.run_process('upload-gtf', inputs)

            inputs = {
                'src': 'annotation dicty.gtf.gz',
                'source': 'DICTYBASE',
                'species': 'Homo sapiens',
                'build': 'dd-05-2009',
            }
            annotation_wrong_species = self.run_process('upload-gtf', inputs)

            inputs = {
                'src': 'annotation dicty.gtf.gz',
                'source': 'DICTYBASE',
                'species': 'Dictyostelium discoideum',
                'build': 'wrong build',
            }
            annotation_wrong_build = self.run_process('upload-gtf', inputs)

            aligned_reads = self.run_process('alignment-hisat2', {
                'genome': genome.pk,
                'reads': reads.pk,
            })

        inputs = {
            'alignments': aligned_reads.pk,
            'gff': annotation_correct.pk,
            'stranded': 'no',
            'id_attribute': 'transcript_id',
        }
        expression = self.run_process('htseq-count', inputs)
        self.assertFile(expression, 'rc', 'reads_rc.tab.gz', compression='gzip')
        self.assertFile(expression, 'fpkm', 'reads_fpkm.tab.gz', compression='gzip')
        self.assertFile(expression, 'exp', 'reads_tpm.tab.gz', compression='gzip')
        self.assertJSON(expression, expression.output['exp_json'], '', 'expression_htseq.json.gz')
        self.assertFile(expression, 'exp_set', 'htseq_count_out_exp_set.txt.gz', compression='gzip')
        self.assertJSON(expression, expression.output['exp_set_json'], '', 'htseq_count_exp_set.json.gz')
        self.assertFields(expression, 'species', 'Dictyostelium discoideum')
        self.assertFields(expression, 'build', 'dd-05-2009')
        self.assertFields(expression, 'feature_type', 'gene')

        inputs['gff'] = annotation_wrong_species.pk
        expression = self.run_process('htseq-count', inputs, Data.STATUS_ERROR)

        inputs['gff'] = annotation_wrong_build.pk
        expression = self.run_process('htseq-count', inputs, Data.STATUS_ERROR)

    @with_resolwe_host
    @tag_process('htseq-count-raw')
    def test_expression_htseq_cpm(self):
        with self.preparation_stage():
            inputs = {
                'src': 'annotation dicty.gtf.gz',
                'source': 'DICTYBASE',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009',
            }
            annotation = self.run_process('upload-gtf', inputs)

            inputs = {
                'src': 'feature_counts hs.gtf.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'GRCh38_ens90',
            }
            annotation_hs = self.run_process('upload-gtf', inputs)

            bam = {
                'src': 'reads.bam',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009',
            }
            bam = self.run_process('upload-bam', bam)

            inputs = {
                'src': 'feature_counts hs_paired.bam',
                'species': 'Homo sapiens',
                'build': 'GRCh38_ens90',
            }
            bam_paired = self.run_process('upload-bam', inputs)

            inputs = {
                'alignments': bam.pk,
                'gtf': annotation.pk,
                'stranded': 'no',
                'id_attribute': 'transcript_id',
            }
        expression = self.run_process('htseq-count-raw', inputs)
        self.assertFile(expression, 'rc', 'reads_rc.tab.gz', compression='gzip')
        self.assertFile(expression, 'exp', 'reads_cpm.tab.gz', compression='gzip')
        self.assertFields(expression, 'species', 'Dictyostelium discoideum')
        self.assertFields(expression, 'build', 'dd-05-2009')

        inputs = {
            'alignments': bam_paired.pk,
            'gtf': annotation_hs.pk,
        }
        expression = self.run_process('htseq-count-raw', inputs)
        self.assertFile(expression, 'rc', 'htseq_raw_rc.tab.gz', compression='gzip')
        self.assertFile(expression, 'exp', 'htseq_raw_cpm.tab.gz', compression='gzip')
        self.assertFile(expression, 'exp_set', 'htseq_cpm_exp_set.txt.gz', compression='gzip')
        self.assertJSON(expression, expression.output['exp_set_json'], '', 'htseq_cpm_exp_set.json.gz')
        self.assertFields(expression, 'species', 'Homo sapiens')
        self.assertFields(expression, 'build', 'GRCh38_ens90')
        self.assertFields(expression, 'feature_type', 'gene')

    @tag_process('index-fasta-nucl')
    def test_index_fasta_nucl(self):
        with self.preparation_stage():
            inputs = {'src': 'HS chr21_ensembl.fa.gz'}
            genome = self.run_process('upload-fasta-nucl', inputs)

            inputs = {
                'src': 'HS chr21_short.gtf.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_90'
            }
            annotation = self.run_process('upload-gtf', inputs)

        inputs = {'nucl': genome.pk, 'annotation': annotation.pk}
        index_fasta_nucl = self.run_process('index-fasta-nucl', inputs)

        del index_fasta_nucl.output['rsem_index']['total_size']  # Non-deterministic output.
        self.assertFields(index_fasta_nucl, 'rsem_index', {'dir': 'rsem'})
        self.assertFields(index_fasta_nucl, 'source', 'ENSEMBL')
        self.assertFields(index_fasta_nucl, 'species', 'Homo sapiens')
        self.assertFields(index_fasta_nucl, 'build', 'ens_90')

    @with_resolwe_host
    @tag_process('mergeexpressions')
    def test_mergeexpression(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM")
            expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")
            expression_3 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="RC")

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067098', 'DPU_G0067102']
        }

        mergeexpression_1 = self.run_process('mergeexpressions', inputs)
        self.assertFile(mergeexpression_1, "expset", "merged_expset_subset.tab")

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': []
        }

        mergeexpression_2 = self.run_process('mergeexpressions', inputs)
        self.assertFile(mergeexpression_2, "expset", "merged_expset_all.tab")

        inputs = {
            'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067098', 'DPU_G0067102']
        }
        self.run_process('mergeexpressions', inputs, Data.STATUS_ERROR)

    @tag_process('mergeetc')
    def test_etcmerge(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            annotation = self.prepare_annotation_gff()

            aligned_reads = self.run_process('alignment-hisat2', {
                'genome': genome.pk,
                'reads': reads.pk,
            })

            mappa = self.run_process("upload-mappability", {"src": "purpureum_mappability_50.tab.gz"})

            inputs = {
                'alignment': aligned_reads.pk,
                'gff': annotation.pk,
                'mappable': mappa.pk}

            expression = self.run_process('expression-dicty', inputs)

            inputs = {'expressions': [expression.pk, expression.pk]}
            etc = self.run_process('etc-bcm', inputs)

        inputs = {
            'exps': [etc.pk],
            'genes': ['DPU_G0067110', 'DPU_G0067098', 'DPU_G0067102']
        }

        etcmerge = self.run_process('mergeetc', inputs)
        self.assertFile(etcmerge, "expset", "merged_etc.tab.gz", compression='gzip')

    @with_resolwe_host
    @tag_process('feature_counts')
    def test_feature_counts(self):
        with self.preparation_stage():
            inputs = {
                'src': 'feature_counts hs.gtf.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'GRCh38_ens90',
            }
            annotation_gtf = self.run_process('upload-gtf', inputs)
            annotation_gff3 = self.prepare_annotation_gff()

            bam_single_inputs = {
                'src': 'reads.bam',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009'
            }
            bam_single = self.run_process('upload-bam', bam_single_inputs)

            inputs = {
                'src': 'feature_counts hs_paired.bam',
                'species': 'Homo sapiens',
                'build': 'GRCh38_ens90',
            }
            bam_paired = self.run_process('upload-bam', inputs)

        inputs = {
            'alignment': {
                'aligned_reads': bam_paired.id,
            },
            'annotation': {
                'annotation': annotation_gtf.id,
            },
        }

        expression = self.run_process('feature_counts', inputs)
        self.assertFile(expression, 'rc', 'feature_counts_out_rc.tab.gz', compression='gzip')
        self.assertFile(expression, 'fpkm', 'feature_counts_out_fpkm.tab.gz', compression='gzip')
        self.assertFile(expression, 'cpm', 'feature_counts_out_cpm.tab.gz', compression='gzip')
        self.assertFile(expression, 'exp', 'feature_counts_out_tpm.tab.gz', compression='gzip')
        self.assertFile(expression, 'exp_set', 'feature_counts_out_exp_set.txt.gz', compression='gzip')
        self.assertJSON(expression, expression.output['exp_set_json'], '', 'feature_counts_exp_set.json.gz')
        self.assertFields(expression, 'species', 'Homo sapiens')
        self.assertFields(expression, 'build', 'GRCh38_ens90')
        self.assertFields(expression, 'feature_type', 'gene')

        inputs = {
            'alignment': {
                'aligned_reads': bam_single.id,
            },
            'annotation': {
                'annotation': annotation_gff3.id,
                'id_attribute': 'Parent',
            },
        }
        expression = self.run_process('feature_counts', inputs)
        self.assertFile(expression, 'rc', 'reads_rc.tab.gz', compression='gzip')
        self.assertFile(expression, 'fpkm', 'reads_fpkm.tab.gz', compression='gzip')
        self.assertFile(expression, 'exp', 'reads_tpm.tab.gz', compression='gzip')
        self.assertFields(expression, 'feature_type', 'gene')

    @with_resolwe_host
    @tag_process('feature_counts')
    def test_feature_counts_rpkum(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            annotation = self.prepare_annotation(fn='annotation dicty.gtf.gz')
            annotation_gff = self.prepare_annotation_gff()

            aligned_reads = self.run_process('alignment-hisat2', {
                'genome': genome.pk,
                'reads': reads.pk
            })

            mappability = self.run_process("mappability-bcm", {
                "genome": genome.id,
                "gff": annotation_gff.id,
                "length": 50,
            })

        feature_counts = self.run_process('feature_counts', {
            'alignment': {
                'aligned_reads': aligned_reads.id,
            },
            'annotation': {
                'annotation': annotation.id,
                'id_attribute': 'transcript_id',
            },
            'normalization_type': 'RPKUM',
            'mappability': mappability.id,
        })
        self.assertFile(feature_counts, 'exp', 'expression_fc_rpkum.tab.gz', compression='gzip')
        self.assertFields(feature_counts, "source", "DICTYBASE")
        self.assertFields(feature_counts, 'species', 'Dictyostelium discoideum')
        self.assertFields(feature_counts, 'build', 'dd-05-2009')
        self.assertFields(feature_counts, 'feature_type', 'gene')

    @tag_process('salmon-index')
    def test_salmon_index(self):
        with self.preparation_stage():
            cds = self.run_process('upload-fasta-nucl', {'src': 'salmon_cds.fa.gz'})

        inputs = {
            'nucl': cds.id,
            'perfect_hash': True,
            'gencode': False,
            'keep_duplicates': True,
            'source': 'ENSEMBL',
            'species': 'Homo sapiens',
            'build': 'ens_90',
        }
        salmon_index = self.run_process('salmon-index', inputs)

        del salmon_index.output['index']['total_size']  # Non-deterministic output.
        self.assertFields(salmon_index, 'index', {'dir': 'salmon_index'})
        self.assertFields(salmon_index, 'source', 'ENSEMBL')
        self.assertFields(salmon_index, 'species', 'Homo sapiens')
        self.assertFields(salmon_index, 'build', 'ens_90')

    @with_resolwe_host
    @tag_process('salmon-quant')
    def test_salmon_quant(self):
        with self.preparation_stage():
            reads = self.prepare_reads([os.path.join('salmon_quant', 'input', 'hs sim_reads_single.fastq.gz')])
            annotation = self.prepare_annotation(
                os.path.join('salmon_quant', 'input', 'hs annotation.gtf.gz'),
                source='ENSEMBL',
                species='Homo sapiens',
                build='ens_92',
            )
            transcripts = self.run_process('upload-fasta-nucl', {
                'src': os.path.join('salmon_quant', 'input', 'hs cdna.fasta.gz'),
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_92',
            })
            salmon_index = self.run_process('salmon-index', {
                'nucl': transcripts.id,
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_92',
            })

        inputs = {
            'reads': reads.id,
            'salmon_index': salmon_index.id,
            'annotation': annotation.id,
            'options': {
                'min_assigned_frag': 5,
                'gc_bias': True,
                'seq_bias': True,
                'validate_mappings': True,
                'range_factorization_bins': 4,
                'incompat_prior': 0.05,
                'min_score_fraction': 0.7,
                'consensus_slack': 0.25,
                'no_length_correction': False,
                'discard_orphans_quasi': True,
            }
        }
        salmon_quant = self.run_process('salmon-quant', inputs)
        self.assertFile(
            salmon_quant,
            'exp_set',
            os.path.join('salmon_quant', 'output', 'salmon_quant_tpm.tab.gz'),
            compression='gzip',
        )

    @with_resolwe_host
    @tag_process('feature_counts')
    def test_featurecounts_strandedness(self):
        with self.preparation_stage():
            cds = self.run_process('upload-fasta-nucl', {'src': 'salmon_cds.fa.gz'})

            salmon_index = self.run_process('salmon-index', {
                'nucl': cds.id,
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_90',
            })

            annotation = self.run_process('upload-gtf', {
                'src': 'annotation_rsem.gtf.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_90',
            })

            aligned_reads = self.run_process('upload-bam', {
                'src': 'feature counts_detect_strandedness.bam',
                'species': 'Homo sapiens',
                'build': 'ens_90',
            })

        inputs = {
            'alignment': {
                'aligned_reads': aligned_reads.id,
                'assay_type': 'auto',
                'cdna_index': salmon_index.id,
            },
            'annotation': {
                'annotation': annotation.id,
            },
        }

        expression = self.run_process('feature_counts', inputs)
        self.assertFile(expression, 'exp', 'auto_detect_strand_tpm.tab.gz', compression='gzip')

    @tag_process('shrna-quant')
    def test_shrna_quant(self):
        with self.preparation_stage():
            pf_in = './shrna_diffexp/input/'
            pf_out = './shrna_diffexp/output/'

            species = 'Homo sapiens'
            build = 'custom-from-file'
            bam_single_inputs = {
                'src': pf_in + 'SM18_ss.bam',
                'species': species,
                'build': build
            }
            bam = self.run_process('upload-bam', bam_single_inputs)

        inputs = {
            'alignment': bam.id,
            'readlengths': 26,
            'alignscores': -6
        }

        quant = self.run_process('shrna-quant', inputs)
        self.assertFile(quant, 'rc', pf_out + 'SM18_ss_count_matrix.txt.gz', compression='gzip')
        self.assertFile(quant, 'exp', pf_out + 'SM18_ss_count_matrix.txt.gz', compression='gzip')
        self.assertFields(quant, 'exp_type', 'RC')
        self.assertJSON(quant, quant.output['exp_json'], '', pf_out + 'SM18_ss_json.txt.gz')
        self.assertFields(quant, 'source', 'shRNA-gene-sequences')
        self.assertFields(quant, 'species', species)
        self.assertFields(quant, 'build', build)
        self.assertFields(quant, 'feature_type', 'shRNA')
        self.assertFile(quant, 'mapped_species', pf_out + 'SM18_ss_mapped_species.txt.gz', compression='gzip')
