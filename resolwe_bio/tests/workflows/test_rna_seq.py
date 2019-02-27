# pylint: disable=missing-docstring,invalid-name
from resolwe.flow.models import Data
from resolwe.test import with_resolwe_host, tag_process

from resolwe_bio.utils.test import KBBioProcessTestCase


class RNASeqWorkflowTestCase(KBBioProcessTestCase):
    @tag_process('workflow-rnaseq-cuffquant')
    def test_cuffquant_workflow(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            annotation = self.prepare_annotation_gff()

        self.run_process(
            'workflow-rnaseq-cuffquant', {
                'reads': reads.id,
                'genome': genome.id,
                'annotation': annotation.id,
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

    @with_resolwe_host
    @tag_process('workflow-bbduk-star-htseq')
    def test_bbduk_star_htseq_single_workflow(self):
        with self.preparation_stage():
            inputs = {'src': ['hs_single bbduk_star_htseq_reads_single.fastq.gz']}
            reads = self.run_processor('upload-fastq-single', inputs)

            inputs = {'src': 'hs genome.fasta.gz'}
            star_index_fasta = self.run_process('upload-fasta-nucl', inputs)
            adapters = self.prepare_adapters()

            inputs = {
                'src': 'hs annotation.gtf.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_90'
            }
            annotation = self.run_process('upload-gtf', inputs)

            inputs = {'annotation': annotation.id, 'genome2': star_index_fasta.id}
            star_index = self.run_process('alignment-star-index', inputs)

        self.run_process(
            'workflow-bbduk-star-htseq', {
                'reads': reads.id,
                'star_index': star_index.id,
                'adapters': [adapters.id],
                'annotation': annotation.id,
                'stranded': 'yes'
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_bbduk_star_htseq_single_rc.tab.gz', compression='gzip')
        self.assertFile(workflow, 'exp', 'workflow_bbduk_star_htseq_single_cpm.tab.gz', compression='gzip')
        self.assertFields(workflow, 'exp_type', 'CPM')
        self.assertFields(workflow, 'source', 'ENSEMBL')
        self.assertFields(workflow, 'species', 'Homo sapiens')

    @with_resolwe_host
    @tag_process('workflow-bbduk-star-htseq-paired')
    def test_bbduk_star_htseq_paired_workflow(self):
        with self.preparation_stage():
            paired_reads = self.prepare_paired_reads(['hs_paired_R1 workflow_bbduk_star_htseq.fastq.gz'],
                                                     ['hs_paired_R2 workflow_bbduk_star_htseq.fastq.gz'])
            inputs = {
                'src': 'hs annotation.gtf.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_90'
            }
            annotation = self.run_process('upload-gtf', inputs)

            star_index_fasta = self.prepare_adapters('hs genome.fasta.gz')
            inputs = {'annotation': annotation.id, 'genome2': star_index_fasta.id}
            star_index = self.run_process('alignment-star-index', inputs)
            adapters = self.prepare_adapters()

        inputs = {
            'reads': paired_reads.id,
            'adapters': [adapters.id],
            'star_index': star_index.id,
            'annotation': annotation.id,
            'stranded': 'reverse',
        }
        self.run_process('workflow-bbduk-star-htseq-paired', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_bbduk_star_htseq_paired_rc.tab.gz', compression='gzip')
        self.assertFile(workflow, 'exp', 'workflow_bbduk_star_htseq_paired_cpm.tab.gz', compression='gzip')
        self.assertFields(workflow, 'source', 'ENSEMBL')
        self.assertFields(workflow, 'species', 'Homo sapiens')

    @with_resolwe_host
    @tag_process('workflow-bbduk-star-fc-quant-single', 'workflow-bbduk-star-fc-quant-paired')
    def test_bbduk_star_fc_quant_workflow(self):
        with self.preparation_stage():
            inputs = {'src': ['hs_single bbduk_star_htseq_reads_single.fastq.gz']}
            reads = self.run_processor('upload-fastq-single', inputs)

            paired_reads = self.prepare_paired_reads(['hs_paired_R1 workflow_bbduk_star_htseq.fastq.gz'],
                                                     ['hs_paired_R2 workflow_bbduk_star_htseq.fastq.gz'])

            inputs = {'src': 'hs genome.fasta.gz'}
            star_index_fasta = self.run_process('upload-fasta-nucl', inputs)
            adapters = self.prepare_adapters()

            inputs = {
                'src': 'hs annotation.gtf.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_90'
            }
            annotation = self.run_process('upload-gtf', inputs)

            inputs = {'annotation': annotation.id, 'genome2': star_index_fasta.id}
            star_index = self.run_process('alignment-star-index', inputs)

            rrna_reference = self.run_process('upload-fasta-nucl', {
                'src': 'Homo_sapiens_rRNA.fasta.gz',
                'source': 'NCBI',
                'species': 'Homo sapiens',
                'build': 'rRNA',
            })
            rrna_star_index = self.run_process('alignment-star-index', {
                'genome2': rrna_reference.id,
                'advanced': {
                    'genomeSAindexNbases': 2,
                }
            })

            globin_reference = self.run_process('upload-fasta-nucl', {
                'src': 'Homo_sapiens_globin_reference.fasta.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'globin',
            })
            globin_star_index = self.run_process('alignment-star-index', {
                'genome2': globin_reference.id,
                'advanced': {
                    'genomeSAindexNbases': 2,
                }
            })

        self.run_process(
            'workflow-bbduk-star-fc-quant-single', {
                'reads': reads.id,
                'star_index': star_index.id,
                'adapters': [adapters.id],
                'annotation': annotation.id,
                'stranded': 'forward',
                'qc': {
                    'rrna_reference': rrna_star_index.id,
                    'globin_reference': globin_star_index.id,
                },
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        feature_counts = Data.objects.filter(process__slug='feature_counts').last()
        self.assertFile(feature_counts, 'rc', 'workflow_bbduk_star_fc_quant_single_rc.tab.gz', compression='gzip')
        self.assertFile(feature_counts, 'exp', 'workflow_bbduk_star_fc_quant_single_cpm.tab.gz', compression='gzip')
        self.assertFields(feature_counts, 'exp_type', 'CPM')
        self.assertFields(feature_counts, 'source', 'ENSEMBL')
        self.assertFields(feature_counts, 'species', 'Homo sapiens')

        self.run_process(
            'workflow-bbduk-star-fc-quant-paired', {
                'reads': paired_reads.id,
                'adapters': [adapters.id],
                'star_index': star_index.id,
                'annotation': annotation.id,
                'stranded': 'reverse',
                'qc': {
                    'rrna_reference': rrna_star_index.id,
                    'globin_reference': globin_star_index.id,
                },
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        fc_paired = Data.objects.filter(process__slug='feature_counts').last()
        self.assertFile(fc_paired, 'rc', 'workflow_bbduk_star_fc_quant_paired_rc.tab.gz', compression='gzip')
        self.assertFile(fc_paired, 'exp', 'workflow_bbduk_star_fc_quant_paired_cpm.tab.g', compression='gzip')
        self.assertFields(fc_paired, 'exp_type', 'CPM')
        self.assertFields(fc_paired, 'source', 'ENSEMBL')
        self.assertFields(fc_paired, 'species', 'Homo sapiens')

    @with_resolwe_host
    @tag_process('workflow-bbduk-star-featurecounts-single', 'workflow-bbduk-star-featurecounts-paired',
                 'workflow-bbduk-star-featurecounts-qc-single', 'workflow-bbduk-star-featurecounts-qc-paired')
    def test_bbduk_star_featurecounts_workflow(self):
        with self.preparation_stage():
            reads = self.prepare_reads(['hs sim_reads_single.fastq.gz'])
            paired_reads = self.prepare_paired_reads(['hs sim_reads1.fastq.gz'], ['hs sim_reads2.fastq.gz'])
            annotation = self.prepare_annotation(
                fn='hs annotation.gtf.gz',
                source='ENSEMBL',
                species='Homo sapiens',
                build='ens90',
            )
            star_index_fasta = self.run_process('upload-fasta-nucl', {
                'src': 'hs genome.fasta.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens90',
            })
            inputs = {
                'annotation': annotation.id,
                'genome2': star_index_fasta.id,
            }
            star_index = self.run_process('alignment-star-index', inputs)
            adapters = self.prepare_adapters()

            rrna_reference = self.run_process('upload-fasta-nucl', {
                'src': 'Homo_sapiens_rRNA.fasta.gz',
                'source': 'NCBI',
                'species': 'Homo sapiens',
                'build': 'rRNA',
            })
            rrna_star_index = self.run_process('alignment-star-index', {
                'genome2': rrna_reference.id,
                'advanced': {
                    'genomeSAindexNbases': 2,
                }
            })

            globin_reference = self.run_process('upload-fasta-nucl', {
                'src': 'Homo_sapiens_globin_reference.fasta.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'globin',
            })
            globin_star_index = self.run_process('alignment-star-index', {
                'genome2': globin_reference.id,
                'advanced': {
                    'genomeSAindexNbases': 2,
                }
            })

        inputs = {
            'preprocessing': {
                'reads': reads.id,
                'adapters': [adapters.id],
                'custom_adapter_sequences': ['ACTGACTGACTG', 'AAACCCTTT'],
            },
            'alignment': {
                'genome': star_index.id,
            },
            'quantification': {
                'annotation': annotation.id,
            },
        }
        self.run_process('workflow-bbduk-star-featurecounts-single', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        feature_counts = Data.objects.get(process__slug='feature_counts')
        self.assertFile(feature_counts, 'rc', 'feature_counts_rc_single.tab.gz', compression='gzip')
        multiqc = Data.objects.get(process__slug='multiqc')
        self.assertFileExists(multiqc, 'report')

        inputs['preprocessing']['reads'] = paired_reads.id
        self.run_process('workflow-bbduk-star-featurecounts-paired', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        feature_counts = Data.objects.filter(process__slug='feature_counts').last()
        self.assertFile(feature_counts, 'rc', 'feature_counts_rc_paired.tab.gz', compression='gzip')
        multiqc = Data.objects.filter(process__slug='multiqc').last()
        self.assertFileExists(multiqc, 'report')

        inputs = {
            'preprocessing': {
                'reads': reads.id,
                'adapters': [adapters.id],
                'custom_adapter_sequences': ['ACTGACTGACTG', 'AAACCCTTT'],
            },
            'alignment': {
                'genome': star_index.id,
            },
            'quantification': {
                'annotation': annotation.id,
            },
            'qc': {
                'rrna_reference': rrna_star_index.id,
                'globin_reference': globin_star_index.id,
            }
        }

        self.run_process('workflow-bbduk-star-featurecounts-qc-single', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        feature_counts = Data.objects.filter(process__slug='feature_counts').last()
        self.assertFile(feature_counts, 'rc', 'feature_counts_rc_single.tab.gz', compression='gzip')
        globin = Data.objects.filter(process__slug='alignment-star').last()
        self.assertFields(globin, 'build', 'globin')
        multiqc = Data.objects.filter(process__slug='multiqc').last()
        self.assertFileExists(multiqc, 'report')

        inputs['preprocessing']['reads'] = paired_reads.id
        self.run_process('workflow-bbduk-star-featurecounts-qc-paired', inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        feature_counts = Data.objects.filter(process__slug='feature_counts').last()
        self.assertFile(feature_counts, 'rc', 'feature_counts_rc_paired.tab.gz', compression='gzip')
        globin = Data.objects.filter(process__slug='alignment-star').last()
        self.assertFields(globin, 'build', 'globin')
        multiqc = Data.objects.filter(process__slug='multiqc').last()
        self.assertFileExists(multiqc, 'report')

    @with_resolwe_host
    @tag_process('workflow-custom-cutadapt-star-htseq-single', 'workflow-custom-cutadapt-star-htseq-paired')
    def test_custom_cutadapt_star_htseq_workflow(self):
        with self.preparation_stage():
            reads = self.prepare_reads(['SRR2124780_1 1k.fastq.gz'])
            paired_reads = self.prepare_paired_reads(mate1=['SRR2124780_1 1k.fastq.gz'],
                                                     mate2=['SRR2124780_2 1k.fastq.gz'])
            annotation = self.prepare_annotation(
                fn='HS chr21_short.gtf.gz',
                source='UCSC',
                species='Homo sapiens',
                build='hg19'
            )
            star_index_fasta = self.prepare_adapters(fn='HS chr21_ensembl.fa.gz')
            inputs = {'annotation': annotation.id, 'genome2': star_index_fasta.id}

            star_index = self.run_process('alignment-star-index', inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.run_process(
            'workflow-custom-cutadapt-star-htseq-single', {
                'reads': reads.id,
                'genome': star_index.id,
                'gff': annotation.id
            }
        )
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_ccshs.tab.gz', compression='gzip')

        self.run_process(
            'workflow-custom-cutadapt-star-htseq-paired', {
                'reads': paired_reads.id,
                'genome': star_index.id,
                'gff': annotation.id
            }
        )
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_ccshp.tab.gz', compression='gzip')

    @with_resolwe_host
    @tag_process('workflow-custom-cutadapt-star-rsem-single', 'workflow-custom-cutadapt-star-rsem-paired')
    def test_custom_cutadapt_star_rsem_workflow(self):
        with self.preparation_stage():
            single_reads = self.prepare_reads(['reads rsem.fq.gz'])
            paired_reads = self.prepare_paired_reads(mate1=['reads rsem.fq.gz'], mate2=['reads rsem2.fq.gz'])

            inputs = {'src': 'genome_rsem.fa.gz'}
            genome = self.run_process('upload-fasta-nucl', inputs)

            inputs = {
                'src': 'annotation_rsem.gtf.gz',
                'source': 'ENSEMBL',
                'species': 'Homo sapiens',
                'build': 'ens_90'
            }
            annotation = self.run_process('upload-gtf', inputs)

            inputs = {'genome2': genome.pk, 'annotation': annotation.pk}
            star_index = self.run_process('alignment-star-index', inputs)

            inputs = {'nucl': genome.pk, 'annotation': annotation.pk}
            index_fasta_nucl = self.run_process('index-fasta-nucl', inputs)

        inputs = {
            'reads': single_reads.pk,
            'star_index': star_index.pk,
            'expression_index': index_fasta_nucl.pk,
            'stranded': 'yes'
        }
        self.run_process('workflow-custom-cutadapt-star-rsem-single', inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        workflow = Data.objects.last()

        self.assertFile(workflow, 'rc', 'workflow_ccsrs.tab.gz', compression='gzip')
        self.assertFile(workflow, 'genes', 'rsem_genes_single.tab.gz', compression='gzip')
        self.assertFile(workflow, 'transcripts', 'rsem_isoforms_single.tab.gz', compression='gzip')

        inputs = {
            'reads': paired_reads.pk,
            'star_index': star_index.pk,
            'expression_index': index_fasta_nucl.pk,
            'stranded': 'yes'
        }
        self.run_process('workflow-custom-cutadapt-star-rsem-paired', inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        workflow = Data.objects.last()

        self.assertFile(workflow, 'rc', 'workflow_ccsrp.tab.gz', compression='gzip')
        self.assertFile(workflow, 'genes', 'rsem_genes_paired.tab.gz', compression='gzip')
        self.assertFile(workflow, 'transcripts', 'rsem_isoforms_paired.tab.gz', compression='gzip')
        self.assertFile(workflow, 'exp_set', 'rsem_paired_exp_set.txt.gz', compression='gzip')
        self.assertJSON(workflow, workflow.output['exp_set_json'], '', 'rsem_paired_exp_set.json.gz')

    @with_resolwe_host
    @tag_process('workflow-rnaseq-single')
    def test_rnaseq_single_workflow(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            single_reads = self.prepare_reads()
            annotation = self.prepare_annotation('annotation dicty.gtf.gz')
            adapters = self.prepare_adapters()

        self.run_process('workflow-rnaseq-single', {
            'genome': genome.id,
            'reads': single_reads.id,
            'annotation': annotation.id,
            'minlen': 10,
            'stranded': 'no',
            'id_attribute': 'transcript_id'
        })

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_rnaseq_single_rc.tab.gz', compression='gzip')
        self.assertFields(workflow, 'exp_type', 'TPM')
        self.assertFields(workflow, 'source', 'DICTYBASE')

        self.run_process('workflow-rnaseq-single', {
            'genome': genome.id,
            'reads': single_reads.id,
            'annotation': annotation.id,
            'adapters': adapters.id,
            'minlen': 10,
            'stranded': 'no',
            'id_attribute': 'transcript_id'
        })

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_rnaseq_single_rc.tab.gz', compression='gzip')
        self.assertFields(workflow, 'exp_type', 'TPM')
        self.assertFields(workflow, 'source', 'DICTYBASE')

    @with_resolwe_host
    @tag_process('workflow-rnaseq-paired')
    def test_rnaseq_paired_workflow(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            paired_reads = self.prepare_paired_reads()
            annotation = self.prepare_annotation('annotation dicty.gtf.gz')
            adapters = self.prepare_adapters()

        self.run_process('workflow-rnaseq-paired', {
            'genome': genome.id,
            'reads': paired_reads.id,
            'annotation': annotation.id,
            'minlen': 10,
            'trailing': 1,
            'stranded': 'no',
            'id_attribute': 'transcript_id'
        })

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_rnaseq_paired_rc.tab.gz', compression='gzip')
        self.assertFields(workflow, 'exp_type', 'TPM')
        self.assertFields(workflow, 'source', 'DICTYBASE')

        self.run_process('workflow-rnaseq-paired', {
            'genome': genome.id,
            'reads': paired_reads.id,
            'annotation': annotation.id,
            'adapters': adapters.id,
            'minlen': 10,
            'trailing': 1,
            'stranded': 'no',
            'id_attribute': 'transcript_id'
        })

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_rnaseq_paired_rc.tab.gz', compression='gzip')
        self.assertFields(workflow, 'exp_type', 'TPM')
        self.assertFields(workflow, 'source', 'DICTYBASE')

    @with_resolwe_host
    @tag_process('workflow-cutadapt-star-featurecounts-single', 'workflow-cutadapt-star-featurecounts-paired')
    def test_cutadapt_star_featurecounts_workflow(self):
        with self.preparation_stage():
            reads = self.prepare_reads(['SRR7455803_1_short.fastq.gz'])
            paired_reads = self.prepare_paired_reads(
                mate1=['SRR7455803_1_short.fastq.gz'],
                mate2=['SRR7455803_2_short.fastq.gz'])

            inputs = {
                'src': 'PGSC_ITAG_test.gff',
                'build': 'stu_dm_itag-pgsc_4.04_nib',
                'source': 'PGSC/ITAG',
                'species': 'Solanum tuberosum',
            }
            gff3 = self.run_process('upload-gff3', inputs)
            star_index_fasta = self.prepare_adapters(fn='StPGSC4.04n_Chr09.fasta.gz')
            inputs = {'annotation': gff3.id, 'genome2': star_index_fasta.id}
            star_index = self.run_process('alignment-star-index', inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.run_process(
            'workflow-cutadapt-star-featurecounts-single', {
                'preprocessing': {
                    'reads': reads.id,
                },
                'alignment': {
                    'genome': star_index.id,
                },
                'quantification': {
                    'annotation': gff3.id
                }
            }
        )
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        featurecounts = Data.objects.last()
        self.assertFile(featurecounts, 'rc', 'workflow_cutadapt_star_fc_single.tab.gz', compression='gzip')

        self.run_process(
            'workflow-cutadapt-star-featurecounts-paired', {
                'preprocessing': {
                    'reads': paired_reads.id,
                },
                'alignment': {
                    'genome': star_index.id,
                },
                'quantification': {
                    'annotation': gff3.id
                }
            }
        )
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        featurecounts = Data.objects.last()
        self.assertFile(featurecounts, 'rc', 'workflow_cutadapt_star_fc_paired.tab.gz', compression='gzip')
