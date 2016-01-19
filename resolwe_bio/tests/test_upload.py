# pylint: disable=missing-docstring
import unittest

from .utils import BioProcessTestCase


class UploadProcessorTestCase(BioProcessTestCase):

    def test_bam_upload(self):
        inputs = {"src": "alignment_name_sorted.bam"}
        upload_bam = self.run_processor("import:upload:mapping-bam", inputs)
        self.assertFiles(upload_bam, 'bam', 'alignment_position_sorted.bam')
        self.assertFiles(upload_bam, 'bai', 'alignment_bam_upload_index.bai')

        inputs = {"src": "alignment_position_sorted.bam", "src2": "alignment_bam_upload_index.bai"}
        try:
            import mongoengine
            self.assertRaises(
                mongoengine.ValidationError, self.run_processor, "import:upload:mapping-bam-indexed", inputs)
        except ImportError:
            # mongoengine is not used on resolwe
            # TODO: update test
            pass

        inputs = {"src": "alignment_position_sorted.bam", "src2": "alignment_bam_upload_index.bam.bai"}
        upload_bam = self.run_processor("import:upload:mapping-bam-indexed", inputs, 'error')
        self.assertFields(upload_bam, 'proc.error', 'BAI should have the same name as BAM with .bai extension')

        inputs = {"src": "alignment_position_sorted.bam", "src2": "alignment_position_sorted.bam.bai"}
        upload_bam = self.run_processor("import:upload:mapping-bam-indexed", inputs)
        self.assertFiles(upload_bam, 'bam', 'alignment_position_sorted.bam')
        self.assertFiles(upload_bam, 'bai', 'alignment_position_sorted.bam.bai')

    def test_upload_expression(self):
        inputs = {"exp_type": "TPM"}
        self.run_processor("import:upload:expression", inputs, 'error')

        inputs = {"exp": "exp_1_tpm.tab.gz", "rc": "exp_1_rc.tab.gz"}
        self.run_processor("import:upload:expression", inputs, 'error')

        inputs = {"rc": "exp_1_rc.tab.gz"}
        exp_3 = self.run_processor("import:upload:expression", inputs)
        self.assertFiles(exp_3, "rc", "exp_1_rc.tab.gz")
        self.assertFiles(exp_3, 'exp', 'exp_1_rc.tab.gz')
        self.assertJSON(exp_3, exp_3.output['exp_json'], '', 'exp_1.json.gz')

        inputs = {"exp": "exp_1_tpm.tab.gz", "exp_type": "TPM"}
        exp_4 = self.run_processor("import:upload:expression", inputs)
        self.assertFiles(exp_4, 'exp', 'exp_1_tpm.tab.gz')

        inputs = {"rc": "exp_1_rc.tab.gz", "exp": "exp_1_tpm.tab.gz", "exp_type": "TPM"}
        exp_5 = self.run_processor("import:upload:expression", inputs)
        self.assertFields(exp_5, 'exp_type', 'TPM')
        self.assertFiles(exp_5, 'exp', 'exp_1_tpm.tab.gz')
        self.assertFiles(exp_5, 'rc', 'exp_1_rc.tab.gz')
        self.assertJSON(exp_5, exp_5.output['exp_json'], '', 'exp_1_norm.json.gz')

        inputs = {"rc": "exp_mac_line_ending.txt.gz"}
        exp_6 = self.run_processor("import:upload:expression", inputs)
        self.assertJSON(exp_6, exp_6.output['exp_json'], '', 'exp.json.gz')

        inputs = {"rc": "exp_unix_line_ending.txt.gz"}
        exp_7 = self.run_processor("import:upload:expression", inputs)
        self.assertJSON(exp_7, exp_7.output['exp_json'], '', 'exp.json.gz')

        inputs = {"rc": "exp_windows_line_ending.txt.gz"}
        exp_8 = self.run_processor("import:upload:expression", inputs)
        self.assertJSON(exp_8, exp_8.output['exp_json'], '', 'exp.json.gz')

    def test_upload_paired_end_reads(self):
        inputs = {"src1": "mate1.fastq.gz", "src2": "mate2.fastq.gz"}
        self.run_processor("import:upload:reads-fastq-paired-end", inputs, 'error')

        inputs = {"src1": "rRNA_forw.fastq.gz", "src2": "rRNA_rew.fastq.gz"}
        reads = self.run_processor("import:upload:reads-fastq-paired-end", inputs)
        self.assertFiles(reads, "fastq", "rRNA_forw.fastq.gz", compression='gzip')
        self.assertFiles(reads, "fastq2", "rRNA_rew.fastq.gz", compression='gzip')
        self.assertFields(reads, "fastqc_archive.file", "rRNA_forw_fastqc.zip")
        self.assertFields(reads, "fastqc_archive2.file", "rRNA_rew_fastqc.zip")
        self.assertFields(reads, "number", 13)
        self.assertFields(reads, "bases", " 101, 101")
        self.assertFields(reads, "fastqc_url.url", "fastqc/rRNA_forw_fastqc/fastqc_report.html")
        self.assertFields(reads, "fastqc_url2.url", "fastqc/rRNA_rew_fastqc/fastqc_report.html")

    def test_upload_single_end_reads(self):
        inputs = {"src": "mate1.fastq.gz"}
        self.run_processor("import:upload:reads-fastq", inputs, 'error')

        inputs = {"src": "rRNA_forw.fastq.gz"}
        reads = self.run_processor("import:upload:reads-fastq", inputs)
        self.assertFiles(reads, "fastq", "rRNA_forw_single.fastq.gz", compression='gzip')
        self.assertFields(reads, "fastqc_archive.file", "rRNA_forw_fastqc.zip")
        self.assertFields(reads, "number", 13)
        self.assertFields(reads, "bases", "101")
        self.assertFields(reads, "fastqc_url.url", "fastqc/rRNA_forw_fastqc/fastqc_report.html")

    def test_upload_reads_old_encoding(self):
        inputs = {"src": "old_encoding.fastq.gz"}
        reads = self.run_processor("import:upload:reads-fastq", inputs)
        self.assertFiles(reads, "fastq", "old_encoding_transformed.fastq.gz", compression='gzip')
        self.assertFields(reads, "fastqc_archive.file", "old_encoding_fastqc.zip")
        self.assertFields(reads, "number", 25)
        self.assertFields(reads, "bases", "40")
        self.assertFields(reads, "fastqc_url.url", "fastqc/old_encoding_fastqc/fastqc_report.html")

    def test_upload_de(self):
        inputs = {'src': 'deseq2_output.tab.gz'}
        diff_exp = self.run_processor("import:upload:diffexp", inputs)

        self.assertFiles(diff_exp, 'diffexp', 'deseq2_output.tab.gz')
        self.assertJSON(diff_exp, diff_exp.output['volcano_plot'], '', 'deseq2_volcano_plot.json.gz')

    def test_upload_genome(self):
        inputs = {"src": "genome.fasta.gz"}
        genome = self.run_processor('import:upload:genome-fasta', inputs)

        self.assertFileExists(genome, "index_bt")
        self.assertFileExists(genome, "index_bt2")
        self.assertFileExists(genome, "index_bwa")

    def test_upload_bed(self):
        inputs = {"src": "bad.bed"}
        bed = self.run_processor('import:upload:bed', inputs, 'error')

        inputs = {"src": "good.bed"}
        bed = self.run_processor('import:upload:bed', inputs)
        self.assertFiles(bed, 'BED', 'good.bed')
