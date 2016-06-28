# pylint: disable=missing-docstring
import unittest

from resolwe.flow.models import Data

from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase


class UploadProcessorTestCase(BioProcessTestCase):

    def test_bam_upload(self):
        inputs = {"src": "alignment_name_sorted.bam"}
        upload_bam = self.run_processor("upload-bam", inputs)
        self.assertFile(upload_bam, 'bam', 'alignment_position_sorted.bam')
        self.assertFile(upload_bam, 'bai', 'alignment_bam_upload_index.bai')

        inputs = {"src": "alignment_position_sorted.bam", "src2": "alignment_bam_upload_index.bam.bai"}
        upload_bam = self.run_processor("upload-bam-indexed", inputs, Data.STATUS_ERROR)
        self.assertEqual(upload_bam.process_error[0], 'BAI should have the same name as BAM with .bai extension')

        inputs = {"src": "alignment_position_sorted.bam", "src2": "alignment_position_sorted.bam.bai"}
        upload_bam = self.run_processor("upload-bam-indexed", inputs)
        self.assertFile(upload_bam, 'bam', 'alignment_position_sorted.bam')
        self.assertFile(upload_bam, 'bai', 'alignment_position_sorted.bam.bai')

    def test_upload_expression(self):
        inputs = {"exp_type": "TPM", 'exp_name': 'Expression'}
        self.run_processor("upload-expression", inputs, Data.STATUS_ERROR)

        inputs = {"exp": "exp_1_tpm.tab.gz", "rc": "exp_1_rc.tab.gz", 'exp_name': 'Expression'}
        self.run_processor("upload-expression", inputs, Data.STATUS_ERROR)

        inputs = {"rc": "exp_1_rc.tab.gz", 'exp_name': 'Expression'}
        exp_3 = self.run_processor("upload-expression", inputs)
        self.assertFile(exp_3, "rc", "exp_1_rc.tab.gz")
        self.assertFile(exp_3, 'exp', 'exp_1_rc.tab.gz')
        self.assertJSON(exp_3, exp_3.output['exp_json'], '', 'exp_1.json.gz')

        inputs = {"exp": "exp_1_tpm.tab.gz", "exp_type": "TPM", 'exp_name': 'Expression'}
        exp_4 = self.run_processor("upload-expression", inputs)
        self.assertFile(exp_4, 'exp', 'exp_1_tpm.tab.gz')

        inputs = {"rc": "exp_1_rc.tab.gz", "exp": "exp_1_tpm.tab.gz", "exp_type": "TPM", 'exp_name': 'Expression'}
        exp_5 = self.run_processor("upload-expression", inputs)
        self.assertFields(exp_5, 'exp_type', 'TPM')
        self.assertFile(exp_5, 'exp', 'exp_1_tpm.tab.gz')
        self.assertFile(exp_5, 'rc', 'exp_1_rc.tab.gz')
        self.assertJSON(exp_5, exp_5.output['exp_json'], '', 'exp_1_norm.json.gz')

        inputs = {"rc": "exp_mac_line_ending.txt.gz", 'exp_name': 'Expression'}
        exp_6 = self.run_processor("upload-expression", inputs)
        self.assertJSON(exp_6, exp_6.output['exp_json'], '', 'exp.json.gz')

        inputs = {"rc": "exp_unix_line_ending.txt.gz", 'exp_name': 'Expression'}
        exp_7 = self.run_processor("upload-expression", inputs)
        self.assertJSON(exp_7, exp_7.output['exp_json'], '', 'exp.json.gz')

        inputs = {"rc": "exp_windows_line_ending.txt.gz", 'exp_name': 'Expression'}
        exp_8 = self.run_processor("upload-expression", inputs)
        self.assertJSON(exp_8, exp_8.output['exp_json'], '', 'exp.json.gz')

    def test_upload_cuffquant_expression(self):
        inputs = {"src": "cuffquant_1.cxb"}
        cxb = self.run_processor("upload-cxb", inputs)

        inputs = {
            'exp': 'cuffquant_exp.tab',
            'cxb': cxb.id}
        self.run_processor('upload-expression-cuffnorm', inputs)

    def test_upload_paired_end_reads(self):
        inputs = {"src1": ["mate1.fastq.gz"], "src2": ["mate2.fastq.gz"]}
        self.run_processor("upload-fastq-paired", inputs, Data.STATUS_ERROR)

        inputs = {"src1": ["rRNA_forw.fastq.gz", "rRNA_rew.fastq.gz"],
                  "src2": ["00Hr.fastq.gz", "20Hr.fastq.gz"]
        }

        reads = self.run_processor("upload-fastq-paired", inputs)
        self.assertFiles(reads, "fastq", ["rRNA_forw.fastq.gz", "rRNA_rew.fastq.gz"], compression='gzip')
        self.assertFiles(reads, "fastq2", ["00Hr.fastq.gz", "20Hr.fastq.gz"], compression='gzip')
        self.assertFields(reads, "fastqc_url", [{'url': 'fastqc/rRNA_forw_fastqc/fastqc_report.html', 'refs': ['fastqc/rRNA_forw_fastqc'], 'name': 'View'}, {'url': 'fastqc/rRNA_rew_fastqc/fastqc_report.html', 'refs': ['fastqc/rRNA_rew_fastqc'], 'name': 'View'}])
        self.assertFields(reads, "fastqc_url2", [{'url': 'fastqc/00Hr_fastqc/fastqc_report.html', 'refs': ['fastqc/00Hr_fastqc'], 'name': 'View'}, {'url': 'fastqc/20Hr_fastqc/fastqc_report.html', 'refs': ['fastqc/20Hr_fastqc'], 'name': 'View'}])

    def test_upload_single_end_reads(self):
        inputs = {"src": ["mate1.fastq.gz"]}
        self.run_processor("upload-fastq-single", inputs, Data.STATUS_ERROR)

        inputs = {"src": ["rRNA_forw.fastq.gz", "rRNA_rew.fastq.gz"]}
        reads = self.run_processor("upload-fastq-single", inputs)

        self.assertFiles(reads, "fastq", ["rRNA_forw_single.fastq.gz", "rRNA_rew.fastq.gz"], compression='gzip')
        self.assertFields(reads, "fastqc_url", [{'url': 'fastqc/rRNA_forw_fastqc/fastqc_report.html', 'refs': ['fastqc/rRNA_forw_fastqc'], 'name': 'View'}, {'url': 'fastqc/rRNA_rew_fastqc/fastqc_report.html', 'refs': ['fastqc/rRNA_rew_fastqc'], 'name': 'View'}])

    def test_upload_reads_old_encoding(self):
        inputs = {"src": ["old_encoding.fastq.gz"]}
        reads = self.run_processor("upload-fastq-single", inputs, Data.STATUS_ERROR)

    def test_upload_de(self):
        inputs = {'src': 'deseq2_output.tab.gz'}
        diff_exp = self.run_processor("upload-diffexp", inputs)

        self.assertFile(diff_exp, 'diffexp', 'deseq2_output.tab.gz')
        self.assertJSON(diff_exp, diff_exp.output['volcano_plot'], '', 'deseq2_volcano_plot.json.gz')

    def test_upload_genome(self):
        inputs = {"src": "genome.fasta.gz"}
        genome = self.run_processor('upload-genome', inputs)

        inputs = {"src": "genome.fasta.gz",
                  "bowtie_index": "bt_index.tar.gz",
                  "bowtie2_index": "bt2_index.tar.gz",
                  "bwa_index": "bwa_index.tar.gz",
                  "hisat2_index": "hisat2_index.tar.gz"}
        genome = self.run_processor('upload-genome', inputs)
        self.assertFields(genome, "index_bt.dir", 'bowtie_index')
        self.assertFields(genome, "index_bt2.dir", 'bowtie2_index')
        self.assertFields(genome, "index_bwa.dir", 'BWA_index')
        self.assertFields(genome, "index_hisat2.dir", 'hisat2_index')

    def test_upload_bed(self):
        inputs = {"src": "bad.bed"}
        bed = self.run_processor('upload-bed', inputs, Data.STATUS_ERROR)

        inputs = {"src": "good.bed"}
        bed = self.run_processor('upload-bed', inputs)
        self.assertFile(bed, 'BED', 'good.bed')
