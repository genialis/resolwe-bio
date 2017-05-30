# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase
from resolwe.flow.models import Data


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

        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_reads_rc.tab.gz', compression='gzip')

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
        workflow = Data.objects.last()
        self.assertFile(workflow, 'rc', 'workflow_rnaseq_paired_rc.tab.gz', compression='gzip')
        self.assertFields(workflow, 'exp_type', 'TPM')
        self.assertFields(workflow, 'source', 'DICTYBASE')
