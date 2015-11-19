# pylint: disable=missing-docstring
import mongoengine

from .utils import ProcessTestCase


class UploadProcessorTestCase(ProcessTestCase):

    def test_bam_upload(self):
        inputs = {"src": "alignment_name_sorted.bam"}
        upload_bam = self.run_processor("import:upload:mapping-bam", inputs)
        self.assertFiles(upload_bam, 'bam', 'alignment_position_sorted.bam')
        self.assertFiles(upload_bam, 'bai', 'alignment_bam_upload_index.bai')

        inputs = {"src": "alignment_position_sorted.bam", "src2": "alignment_bam_upload_index.bai"}
        self.assertRaises(mongoengine.ValidationError, self.run_processor, "import:upload:mapping-bam-indexed", inputs)

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

        inputs = {"exp": "00Hr_tpm.tab.gz", "rc": "00Hr_rc.tab.gz"}
        self.run_processor("import:upload:expression", inputs, 'error')

        inputs = {"rc": "00Hr_rc.tab.gz"}
        expressions_3 = self.run_processor("import:upload:expression", inputs)
        self.assertFiles(expressions_3, "rc", "00Hr_rc.tab.gz")
        self.assertFiles(expressions_3, 'exp', '00Hr_tpm_3.tab.gz')

        inputs = {"exp": "00Hr_tpm.tab.gz", "exp_type": "TPM"}
        expressions_4 = self.run_processor("import:upload:expression", inputs)
        self.assertFiles(expressions_4, 'exp', '00Hr_tpm_2.tab.gz')

        inputs = {"rc": "00Hr_rc.tab.gz", "exp": "00Hr_tpm.tab.gz", "exp_type": "TPM"}
        expressions_5 = self.run_processor("import:upload:expression", inputs)

        self.assertFields(expressions_5, 'exp_type', 'TPM')
        self.assertFiles(expressions_5, 'exp', 'expression_tpm.tab.gz')
        self.assertFiles(expressions_5, 'rc', 'expression_rc.tab.gz')
        self.assertJSON(expressions_5, expressions_5.output['exp_json'], '', 'expression.json.gz')

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

        self.assertFileExist(genome, "index_bt")
        self.assertFileExist(genome, "index_bt2")
        self.assertFileExist(genome, "index_bwa")

    def test_upload_bed(self):
        inputs = {"src": "bad.bed"}
        bed = self.run_processor('import:upload:bed', inputs, 'error')

        inputs = {"src": "good.bed"}
        bed = self.run_processor('import:upload:bed', inputs)
        self.assertFiles(bed, 'BED', 'good.bed')
