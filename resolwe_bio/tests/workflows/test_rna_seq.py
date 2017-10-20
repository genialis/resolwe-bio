# pylint: disable=missing-docstring,invalid-name
from django.contrib.auth.models import AnonymousUser
from django.test import LiveServerTestCase

from guardian.shortcuts import assign_perm

from resolwe_bio.utils.test import BioProcessTestCase
from resolwe.flow.models import Data
from resolwe.test.utils import with_resolwe_host


class RNASeqWorkflowTestCase(BioProcessTestCase):
    def test_cuffquant_workflow(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff.gz', 'source': 'DICTYBASE'}
        annotation = self.run_process('upload-gff3', inputs)

        self.run_process(
            'workflow-rnaseq-cuffquant', {
                'reads': reads.id,
                'genome': genome.id,
                'annotation': annotation.id,
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

    def test_bbduk_star_htseq_workflow(self):
        inputs = {'src': ['workflow_bbduk_star test.fastq.gz']}
        reads = self.run_processor('upload-fastq-single', inputs)

        inputs = {'src': 'HS_chr21_ensemble.fa.gz'}
        star_index_fasta = self.run_process('upload-fasta-nucl', inputs)

        inputs = {'src': 'poly_A.fa.gz'}
        adapters1 = self.run_process('upload-fasta-nucl', inputs)

        inputs = {'src': 'HS_chr21_short.gtf.gz', 'source': 'ENSEMBL'}
        annotation = self.run_process('upload-gtf', inputs)

        inputs = {'annotation': annotation.id, 'genome2': star_index_fasta.id}
        star_index = self.run_process('alignment-star-index', inputs)

        self.run_process(
            'workflow-bbduk-star-htseq', {
                'reads': reads.id,
                'star_index': star_index.id,
                'bbduk_adapters': [adapters1.id],
                'annotation': annotation.id,
                'stranded': 'yes'
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_reads_rc.tab.gz', compression='gzip')

    def test_custom_cutadapt_star_htseq_workflow(self):
        reads = self.prepare_reads(['SRR2124780_1 1k.fastq.gz'])
        paired_reads = self.prepare_paired_reads(mate1=['SRR2124780_1 1k.fastq.gz'],
                                                 mate2=['SRR2124780_2 1k.fastq.gz'])
        annotation = self.prepare_annotation(fn='HS_chr21_short.gtf.gz', source='DICTYBASE')
        star_index_fasta = self.prepare_adapters(fn='HS_chr21_ensemble.fa.gz')
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

    def test_custom_cutadapt_star_rsem_workflow(self):
        inputs = {'src': 'genome_rsem.fa.gz'}
        genome = self.run_process('upload-fasta-nucl', inputs)

        inputs = {'src': 'annotation_rsem.gtf.gz', 'source': 'ENSEMBL'}
        annotation = self.run_process('upload-gtf', inputs)

        inputs = {'genome2': genome.pk, 'annotation': annotation.pk}
        star_index = self.run_process('alignment-star-index', inputs)

        inputs = {'nucl': genome.pk, 'annotation': annotation.pk}
        index_fasta_nucl = self.run_process('index-fasta-nucl', inputs)

        inputs = {'src': ['reads_rsem.fq.gz']}
        single_reads = self.run_process('upload-fastq-single', inputs)

        inputs = {'src1': ['reads_rsem.fq.gz'], 'src2': ['reads_rsem2.fq.gz']}
        paired_reads = self.run_process('upload-fastq-paired', inputs)

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

    def test_rnaseq_single_workflow(self):
        genome = self.prepare_genome()
        single_reads = self.prepare_reads()
        annotation = self.prepare_annotation('annotation.gtf.gz')
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

    def test_rnaseq_paired_workflow(self):
        genome = self.prepare_genome()
        paired_reads = self.prepare_paired_reads()
        annotation = self.prepare_annotation('annotation.gtf.gz')
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


class RNASeqDSSTestCase(BioProcessTestCase, LiveServerTestCase):
    @with_resolwe_host
    def test_rnaseq_dss(self):
        single_reads = self.prepare_reads()
        paired_reads = self.prepare_paired_reads()
        genome = self.prepare_genome()
        genome.slug = 'genome-mm10'
        genome.save()
        assign_perm('view_data', AnonymousUser(), genome)
        annotation = self.prepare_annotation('annotation.gtf.gz')
        annotation.slug = 'annotation-mm10'
        annotation.save()
        assign_perm('view_data', AnonymousUser(), annotation)
        adapters = self.prepare_adapters()
        adapters.slug = 'adapters-illumina'
        adapters.save()
        assign_perm('view_data', AnonymousUser(), adapters)

        self.run_process('dss-rna-seq', {
            'genome_and_annotation': 'mm',
            'reads': single_reads.id,
            'id_attribute': 'transcript_id',
            'adapters': 'no'
        })

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_rnaseq_single_rc.tab.gz', compression='gzip')
        self.assertFields(workflow, 'exp_type', 'TPM')
        self.assertFields(workflow, 'source', 'DICTYBASE')

        self.run_process('dss-rna-seq', {
            'genome_and_annotation': 'mm',
            'reads': paired_reads.id,
            'trailing': 1,
            'id_attribute': 'transcript_id',
            'adapters': 'no'
        })

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_rnaseq_paired_rc.tab.gz', compression='gzip')
        self.assertFields(workflow, 'exp_type', 'TPM')
        self.assertFields(workflow, 'source', 'DICTYBASE')
