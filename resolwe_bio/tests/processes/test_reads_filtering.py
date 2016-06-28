# pylint: disable=missing-docstring
from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase


class ReadsFilteringProcessorTestCase(BioProcessTestCase):

    def test_prinseq_single(self):
        reads = self.prepare_reads()
        inputs = {'reads': reads.pk}

        filtered_reads = self.run_processor('prinseq-lite-single', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_prinseq_single.fastq.gz'], compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'url': 'fastqc/reads_fastqc/fastqc_report.html', 'refs': ['fastqc/reads_fastqc'], 'name': 'View'}])

    def test_prinseq_paired(self):
        inputs = {
            'src1': ['rRNA_forw.fastq.gz'],
            'src2': ['rRNA_rew.fastq.gz']}
        reads = self.run_processor('upload-fastq-paired', inputs)

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_processor('prinseq-lite-paired', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_prinseq_paired_fw.fastq.gz'], compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['filtered_reads_prinseq_paired_rw.fastq.gz'], compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'url': 'fastqc/rRNA_forw_fastqc/fastqc_report.html', 'refs': ['fastqc/rRNA_forw_fastqc'], 'name': 'View'}])
        self.assertFields(filtered_reads, "fastqc_url2", [{'url': 'fastqc/rRNA_rew_fastqc/fastqc_report.html', 'refs': ['fastqc/rRNA_rew_fastqc'], 'name': 'View'}])

    def test_fastqmcf_single(self):
        reads = self.prepare_reads()
        inputs = {'reads': reads.pk}
        filtered_reads = self.run_processor('fastq-mcf-single', inputs)

        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_fastqmcf_single.fastq.gz'], compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'url': 'fastqc/reads_fastqc/fastqc_report.html', 'refs': ['fastqc/reads_fastqc'], 'name': 'View'}])

    def test_fastqmcf_paired(self):
        inputs = {
            'src1': ['rRNA_forw.fastq.gz'],
            'src2': ['rRNA_rew.fastq.gz']}
        reads = self.run_processor('upload-fastq-paired', inputs)

        inputs = {'reads': reads.pk}
        filtered_reads = self.run_processor('fastq-mcf-paired', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['filtered_reads_fastqmcf_paired_fw.fastq.gz'], compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['filtered_reads_fastqmcf_paired_rw.fastq.gz'], compression='gzip')
        self.assertFields(filtered_reads, "fastqc_url", [{'url': 'fastqc/rRNA_forw_fastqc/fastqc_report.html', 'refs': ['fastqc/rRNA_forw_fastqc'], 'name': 'View'}])
        self.assertFields(filtered_reads, "fastqc_url2", [{'url': 'fastqc/rRNA_rew_fastqc/fastqc_report.html', 'refs': ['fastqc/rRNA_rew_fastqc'], 'name': 'View'}])

    def test_sortmerna_single(self):
        reads = self.prepare_reads(['rRNA_forw.fastq.gz'])
        rRNAdb_1 = self.run_processor('upload-fasta-nucl', {'src': 'silva-arc-16s-id95.fasta.gz'})
        rRNAdb_2 = self.run_processor('upload-fasta-nucl', {'src': 'silva-arc-23s-id98.fasta.gz'})

        inputs = {
            'reads': reads.id,
            'database_selection': [rRNAdb_1.id, rRNAdb_2.id],
            'options': {'threads': 2, 'sam': True}}
        filtered_reads = self.run_processor('sortmerna-single', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['reads_wo_rRNA_single.fastq.gz'], compression='gzip')
        self.assertFile(filtered_reads, 'fastq_rRNA', 'reads_rRNA_single.fastq.gz', compression='gzip')
        self.assertFields(filtered_reads, 'fastq_rRNA_sam.file', 'rRNA_forw_rRNA.sam')
        self.assertFields(filtered_reads, 'stats.file', 'stats.log')
        self.assertFields(filtered_reads, "fastqc_url", [{'url': 'fastqc/rRNA_forw_filtered_fastqc/fastqc_report.html', 'refs': ['fastqc/rRNA_forw_filtered_fastqc'], 'name': 'View'}])

    def test_sortmerna_paired(self):
        inputs = {
            'src1': ['rRNA_forw.fastq.gz'],
            'src2': ['rRNA_rew.fastq.gz']}
        reads = self.run_processor('upload-fastq-paired', inputs)

        rRNAdb_1 = self.run_processor('upload-fasta-nucl', {'src': 'silva-arc-16s-id95.fasta.gz'})
        rRNAdb_2 = self.run_processor('upload-fasta-nucl', {'src': 'silva-arc-23s-id98.fasta.gz'})

        inputs = {
            'reads': reads.id,
            'database_selection': [rRNAdb_1.id, rRNAdb_2.id],
            'options': {'threads': 2, 'sam': True}}

        filtered_reads = self.run_processor('sortmerna-paired', inputs)
        self.assertFiles(filtered_reads, 'fastq', ['reads_wo_rRNA_paired_forw.fastq.gz'], compression='gzip')
        self.assertFiles(filtered_reads, 'fastq2', ['reads_wo_rRNA_paired_rew.fastq.gz'], compression='gzip')
        self.assertFile(filtered_reads, 'fastq_rRNA', 'reads_rRNA_paired.fastq.gz', compression='gzip')
        self.assertFields(filtered_reads, 'fastq_rRNA_sam.file', 'rRNA_forw_rRNA.sam')
        self.assertFields(filtered_reads, 'stats.file', 'stats.log')
        self.assertFields(filtered_reads, "fastqc_url", [{'url': 'fastqc/rRNA_forw_filtered_fastqc/fastqc_report.html', 'refs': ['fastqc/rRNA_forw_filtered_fastqc'], 'name': 'View'}])
        self.assertFields(filtered_reads, "fastqc_url2", [{'url': 'fastqc/rRNA_rew_filtered_fastqc/fastqc_report.html', 'refs': ['fastqc/rRNA_rew_filtered_fastqc'], 'name': 'View'}])
