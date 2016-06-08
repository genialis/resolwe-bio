# pylint: disable=missing-docstring
from .utils import skipDockerFailure, BioProcessTestCase


class ReadsFilteringProcessorTestCase(BioProcessTestCase):

    @skipDockerFailure("Fails with: prinseq-lite.pl: command not found")
    def test_prinseq_single(self):
        reads = self.prepare_reads()
        inputs = {'reads': reads.pk}

        filtered_reads = self.run_processor('prinseq-lite-single', inputs)
        self.assertFile(filtered_reads, 'fastq', 'filtered_reads_prinseq_single.fastq.gz', compression='gzip')
        self.assertFields(filtered_reads, 'number', 18)
        self.assertFields(filtered_reads, 'bases', '35')
        self.assertFields(filtered_reads, 'fastqc_url.url', u'fastqc/reads_fastqc/fastqc_report.html')
        self.assertFields(filtered_reads, 'fastqc_archive.file', u'reads_fastqc.zip')

    @skipDockerFailure("Fails with: prinseq-lite.pl: command not found")
    def test_prinseq_paired(self):
        inputs = {
            'src1': 'rRNA_forw.fastq.gz',
            'src2': 'rRNA_rew.fastq.gz'}
        reads = self.run_processor('upload-fastq-paired', inputs)
        inputs = {'reads': reads.pk}

        filtered_reads = self.run_processor('prinseq-lite-paired', inputs)
        self.assertFile(filtered_reads, 'fastq', 'filtered_reads_prinseq_paired_fw.fastq.gz', compression='gzip')
        self.assertFile(filtered_reads, 'fastq2', 'filtered_reads_prinseq_paired_rw.fastq.gz', compression='gzip')
        self.assertFields(filtered_reads, 'number', 13)
        self.assertFields(filtered_reads, 'bases', u' 34-101, 37-101')
        self.assertFields(filtered_reads, 'fastqc_url.url', u'fastqc/rRNA_forw_fastqc/fastqc_report.html')
        self.assertFields(filtered_reads, 'fastqc_url2.url', u'fastqc/rRNA_rew_fastqc/fastqc_report.html')
        self.assertFields(filtered_reads, 'fastqc_archive.file', u'rRNA_forw_fastqc.zip')
        self.assertFields(filtered_reads, 'fastqc_archive2.file', u'rRNA_rew_fastqc.zip')

    @skipDockerFailure("Errors with: KeyError: u'adapters' at "
        "filtered_reads = self.run_processor('fastq-mcf-single', inputs)")
    def test_fastqmcf_single(self):
        reads = self.prepare_reads()
        inputs = {'reads': reads.pk}
        filtered_reads = self.run_processor('fastq-mcf-single', inputs)

        self.assertFile(filtered_reads, 'fastq', 'filtered_reads_fastqmcf_single.fastq.gz', compression='gzip')
        self.assertFields(filtered_reads, 'number', 20)
        self.assertFields(filtered_reads, 'bases', u'33-35')
        self.assertFields(filtered_reads, 'fastqc_url.url', u'fastqc/reads_fastqc/fastqc_report.html')
        self.assertFields(filtered_reads, 'fastqc_archive.file', u'reads_fastqc.zip')

    @skipDockerFailure("Errors with: KeyError: u'adapters' at "
        "filtered_reads = self.run_processor('fastq-mcf-paired', inputs)")
    def test_fastqmcf_paired(self):
        inputs = {
            'src1': 'rRNA_forw.fastq.gz',
            'src2': 'rRNA_rew.fastq.gz'}
        reads = self.run_processor('upload-fastq-paired', inputs)

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_processor('fastq-mcf-paired', inputs)
        self.assertFile(filtered_reads, 'fastq', 'filtered_reads_fastqmcf_paired_fw.fastq.gz', compression='gzip')
        self.assertFile(filtered_reads, 'fastq2', 'filtered_reads_fastqmcf_paired_rw.fastq.gz', compression='gzip')
        self.assertFields(filtered_reads, 'number', 13)
        self.assertFields(filtered_reads, 'bases', u' 34-101, 37-101')
        self.assertFields(filtered_reads, 'fastqc_url.url', u'fastqc/rRNA_forw_fastqc/fastqc_report.html')
        self.assertFields(filtered_reads, 'fastqc_url2.url', u'fastqc/rRNA_rew_fastqc/fastqc_report.html')
        self.assertFields(filtered_reads, 'fastqc_archive.file', u'rRNA_forw_fastqc.zip')
        self.assertFields(filtered_reads, 'fastqc_archive2.file', u'rRNA_rew_fastqc.zip')

    @skipDockerFailure("Fails with: sortmerna_database.py: error: too few "
        "arguments (slugs_path not given)")
    def test_sortmerna_single(self):
        reads = self.prepare_reads('rRNA_forw.fastq.gz')
        inputs = {
            'reads': reads.pk,
            'database_selection': ['rfam-5.8s-database-id98.fasta'],
            'options': {'threads': 2, 'sam': True}}
        filtered_reads = self.run_processor('sortmerna-single', inputs)
        self.assertFile(filtered_reads, 'fastq', 'reads_wo_rRNA_single.fastq.gz', compression='gzip')
        self.assertFile(filtered_reads, 'fastq_rRNA', 'reads_rRNA_single.fastq.gz', compression='gzip')
        self.assertFields(filtered_reads, 'fastq_rRNA_sam.file', 'rRNA_forw_rRNA.sam')
        self.assertFields(filtered_reads, 'stats.file', 'stats.log')
        self.assertFields(filtered_reads, 'number', 13)
        self.assertFields(filtered_reads, 'bases', '101')
        self.assertFields(filtered_reads, 'fastqc_url.url', u'fastqc/rRNA_forw_fastqc/fastqc_report.html')
        self.assertFields(filtered_reads, 'fastqc_archive.file', u'rRNA_forw_fastqc.zip')

    @skipDockerFailure("Fails with: sortmerna_database.py: error: too few "
        "arguments (slugs_path not given)")
    def test_sortmerna_paired(self):
        inputs = {
            'src1': 'rRNA_forw.fastq.gz',
            'src2': 'rRNA_rew.fastq.gz'}
        reads = self.run_processor('upload-fastq-paired', inputs)
        inputs = {
            'reads': reads.pk,
            'database_selection': ['rfam-5.8s-database-id98.fasta'],
            'options': {'threads': 2, 'sam': True}}

        filtered_reads = self.run_processor('sortmerna-paired', inputs)
        self.assertFile(filtered_reads, 'fastq', 'reads_wo_rRNA_paired_forw.fastq.gz', compression='gzip')
        self.assertFile(filtered_reads, 'fastq2', 'reads_wo_rRNA_paired_rew.fastq.gz', compression='gzip')
        self.assertFile(filtered_reads, 'fastq_rRNA', 'reads_rRNA_paired.fastq.gz')
        self.assertFields(filtered_reads, 'fastq_rRNA_sam.file', 'rRNA_forw_rRNA.sam')
        self.assertFields(filtered_reads, 'stats.file', 'stats.log')
        self.assertFields(filtered_reads, 'number', 13)
        self.assertFields(filtered_reads, 'bases', u' 101, 101')
        self.assertFields(filtered_reads, 'fastqc_url.url', "fastqc/rRNA_forw_fastqc/fastqc_report.html")
        self.assertFields(filtered_reads, 'fastqc_url2.url', "fastqc/rRNA_rew_fastqc/fastqc_report.html")
        self.assertFields(filtered_reads, 'fastqc_archive.file', u'rRNA_forw_fastqc.zip')
        self.assertFields(filtered_reads, 'fastqc_archive2.file', u'rRNA_rew_fastqc.zip')
