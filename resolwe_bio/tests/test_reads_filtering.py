from .base import BaseProcessorTestCase
from .utils import PreparedData


class ReadsFilteringProcessorTestCase(BaseProcessorTestCase, PreparedData):

    def test_sormerna_single(self):
        reads = self.prepare_reads('rRNA_forw.fastq.gz')

        inputs = {
            'reads': reads.pk,
            'database_selection': ['rfam-5s-database-id98.fasta'],
            'options': {'threads': 2}}
        filtered_reads = self.run_processor('filtering:sortmerna-2.0-single-end', inputs)
        self.assertFiles(filtered_reads, 'fastq', 'reads_wo_rRNA_single.fastq.gz')
        self.assertFiles(filtered_reads, 'fastq_rRNA', 'reads_rRNA_single.fastq.gz')

    def test_sortmerna_paired(self):
        inputs = {
            'src1': 'rRNA_forw.fastq.gz',
            'src2': 'rRNA_rew.fastq.gz'}
        reads = self.run_processor('import:upload:reads-fastq-paired-end', inputs)
        self.assertFields(reads, 'bases', u' 101, 101')

        inputs = {
            'reads': reads.pk,
            'database_selection': ['rfam-5s-database-id98.fasta'],
            'options': {'threads': 2}}
        filtered_reads = self.run_processor('filtering:sortmerna-2.0-paired-end', inputs)
        self.assertFiles(filtered_reads, 'fastq', 'reads_wo_rRNA_paired_forw.fastq.gz')
        self.assertFiles(filtered_reads, 'fastq2', 'reads_wo_rRNA_paired_rew.fastq.gz')
        self.assertFiles(filtered_reads, 'fastq_rRNA', 'reads_rRNA_paired.fastq.gz')

    def test_prinseq_single(self):
        reads = self.prepare_reads()

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_processor('filtering:prinseq-lite-20.4:single-end', inputs)
        self.assertFiles(filtered_reads, 'fastq', 'filtered_reads_prinseq_single.fastq.gz')

    def test_prinseq_paired(self):
        inputs = {
            'src1': 'rRNA_forw.fastq.gz',
            'src2': 'rRNA_rew.fastq.gz'}
        reads = self.run_processor('import:upload:reads-fastq-paired-end', inputs)

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_processor('filtering:prinseq-lite-20.4:paired-end', inputs)
        self.assertFiles(filtered_reads, 'fastq', 'filtered_reads_prinseq_paired_fw.fastq.gz')
        self.assertFiles(filtered_reads, 'fastq2', 'filtered_reads_prinseq_paired_rw.fastq.gz')

    def test_fastqmcf_single(self):
        reads = self.prepare_reads()

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_processor('filtering:fastq-mcf-1.1.2.537:single-end', inputs)
        self.assertFiles(filtered_reads, 'fastq', 'filtered_reads_fastqmcf_single.fastq.gz')

    def test_fastqmcf_paired(self):
        inputs = {
            'src1': 'rRNA_forw.fastq.gz',
            'src2': 'rRNA_rew.fastq.gz'}
        reads = self.run_processor('import:upload:reads-fastq-paired-end', inputs)

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_processor('filtering:fastq-mcf-1.1.2.537:paired-end', inputs)
        self.assertFiles(filtered_reads, 'fastq', 'filtered_reads_fastqmcf_paired_fw.fastq.gz')
        self.assertFiles(filtered_reads, 'fastq2', 'filtered_reads_fastqmcf_paired_rw.fastq.gz')
