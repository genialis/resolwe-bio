# pylint: disable=missing-docstring
import six

from resolwe.flow.models import Data
from resolwe_bio.utils.test import BioProcessTestCase


class UploadProcessorTestCase(BioProcessTestCase):

    def test_bam_upload(self):
        inputs = {'src': 'alignment_name_sorted.bam'}
        upload_bam = self.run_process('upload-bam', inputs)
        self.assertFile(upload_bam, 'bam', 'alignment_position_sorted.bam')
        self.assertFile(upload_bam, 'bai', 'alignment_bam_upload_index.bai')

        inputs = {'src': 'alignment_position_sorted.bam', 'src2': 'alignment_bam_upload_index.bam.bai'}
        upload_bam = self.run_process('upload-bam-indexed', inputs, Data.STATUS_ERROR)
        self.assertEqual(upload_bam.process_error[0], 'BAI should have the same name as BAM with .bai extension')

        inputs = {'src': 'alignment_position_sorted.bam', 'src2': 'alignment_position_sorted.bam.bai'}
        upload_bam = self.run_process('upload-bam-indexed', inputs)
        self.assertFile(upload_bam, 'bam', 'alignment_position_sorted.bam')
        self.assertFile(upload_bam, 'bai', 'alignment_position_sorted.bam.bai')

    def test_upload_expression(self):
        inputs = {'exp_type': 'TPM', 'exp_name': 'Expression', 'source': 'UCSC'}
        self.run_process('upload-expression', inputs, Data.STATUS_ERROR)

        inputs = {'exp': 'exp_1_tpm.tab.gz',
                  'rc': 'exp_1_rc.tab.gz',
                  'exp_name': 'Expression',
                  'source': 'UCSC'}
        self.run_process('upload-expression', inputs, Data.STATUS_ERROR)

        inputs = {'rc': 'exp_1_rc.tab.gz', 'exp_name': 'Expression', 'source': 'UCSC'}
        exp_3 = self.run_process('upload-expression', inputs)
        self.assertFile(exp_3, 'rc', 'exp_1_rc.tab.gz')
        self.assertFile(exp_3, 'exp', 'exp_1_rc.tab.gz')
        self.assertJSON(exp_3, exp_3.output['exp_json'], '', 'exp_1.json.gz')

        inputs = {'exp': 'exp_1_tpm.tab.gz',
                  'exp_type': 'TPM',
                  'exp_name': 'Expression',
                  'source': 'UCSC'}
        exp_4 = self.run_process('upload-expression', inputs)
        self.assertFile(exp_4, 'exp', 'exp_1_tpm.tab.gz')

        inputs = {'rc': 'exp_1_rc.tab.gz',
                  'exp': 'exp_1_tpm.tab.gz',
                  'exp_type': 'TPM',
                  'exp_name': 'Expression',
                  'source': 'UCSC'}
        exp_5 = self.run_process('upload-expression', inputs)
        self.assertFields(exp_5, 'exp_type', 'TPM')
        self.assertFile(exp_5, 'exp', 'exp_1_tpm.tab.gz')
        self.assertFile(exp_5, 'rc', 'exp_1_rc.tab.gz')
        self.assertJSON(exp_5, exp_5.output['exp_json'], '', 'exp_1_norm.json.gz')

        inputs = {'rc': 'exp_mac_line_ending.txt.gz', 'exp_name': 'Expression', 'source': 'UCSC'}
        exp_6 = self.run_process('upload-expression', inputs)
        self.assertJSON(exp_6, exp_6.output['exp_json'], '', 'exp.json.gz')

        inputs = {'rc': 'exp_unix_line_ending.txt.gz', 'exp_name': 'Expression', 'source': 'UCSC'}
        exp_7 = self.run_process('upload-expression', inputs)
        self.assertJSON(exp_7, exp_7.output['exp_json'], '', 'exp.json.gz')

        inputs = {'rc': 'exp_windows_line_ending.txt.gz',
                  'exp_name': 'Expression',
                  'source': 'UCSC'}
        exp_8 = self.run_process('upload-expression', inputs)
        self.assertJSON(exp_8, exp_8.output['exp_json'], '', 'exp.json.gz')

    def test_upload_cuffquant_expr(self):
        inputs = {'src': 'cuffquant_1.cxb', 'source': 'UCSC'}
        cxb = self.run_process('upload-cxb', inputs)

        inputs = {
            'exp': 'cuffquant_exp.tab',
            'cxb': cxb.id}
        self.run_process('upload-expression-cuffnorm', inputs)

    def test_upload_paired_end_reads(self):
        inputs = {'src1': ['mate1.fastq.gz'], 'src2': ['mate2.fastq.gz']}
        self.run_process('upload-fastq-paired', inputs, Data.STATUS_ERROR)

        inputs = {'src1': ['rRNA_forw.fastq.gz', 'rRNA_rew.fastq.gz'],
                  'src2': ['00Hr.fastq.gz', '20Hr.fastq.gz']}

        reads = self.run_process('upload-fastq-paired', inputs)
        self.assertFiles(reads, 'fastq', ['rRNA_forw.fastq.gz', 'rRNA_rew.fastq.gz'], compression='gzip')
        self.assertFiles(reads, 'fastq2', ['00Hr.fastq.gz', '20Hr.fastq.gz'], compression='gzip')
        self.assertFields(reads, 'fastqc_url', [{'file': 'fastqc/rRNA_forw_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/rRNA_forw_fastqc'],
                                                 'size': 343222},
                                                {'file': 'fastqc/rRNA_rew_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/rRNA_rew_fastqc'],
                                                 'size': 323297}])
        self.assertFields(reads, 'fastqc_url2', [{'file': 'fastqc/00Hr_fastqc/fastqc_report.html',
                                                  'refs': ['fastqc/00Hr_fastqc'],
                                                  'size': 327878},
                                                 {'file': 'fastqc/20Hr_fastqc/fastqc_report.html',
                                                  'refs': ['fastqc/20Hr_fastqc'],
                                                  'size': 287245}])

    def test_upload_single_end_reads(self):
        inputs = {'src': ['mate1.fastq.gz']}
        self.run_process('upload-fastq-single', inputs, Data.STATUS_ERROR)

        inputs = {'src': ['rRNA_forw.fastq.gz', 'rRNA_rew.fastq.gz']}
        reads = self.run_process('upload-fastq-single', inputs)

        self.assertFiles(reads, 'fastq', ['rRNA_forw_single.fastq.gz', 'rRNA_rew.fastq.gz'], compression='gzip')
        self.assertFields(reads, 'fastqc_url', [{'file': 'fastqc/rRNA_forw_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/rRNA_forw_fastqc'],
                                                 'size': 343222},
                                                {'file': 'fastqc/rRNA_rew_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/rRNA_rew_fastqc'],
                                                 'size': 323297}])

    def test_upload_reads_old_encoding(self):
        inputs = {'src': ['old_encoding.fastq.gz']}
        self.run_process('upload-fastq-single', inputs, Data.STATUS_ERROR)

    def test_upload_de(self):
        inputs = {
            'src': 'deseq2_output.tab.gz',
            'source': 'DICTYBASE',
            'gene_id': 'gene_id',
            'logfc': 'log2FoldChange',
            'fdr': 'padj',
            'pvalue': 'pvalue',
            'stat': 'stat',
        }
        diff_exp = self.run_process('upload-diffexp', inputs)

        self.assertFile(diff_exp, 'raw', 'deseq2_output.tab.gz')
        self.assertJSON(diff_exp, diff_exp.output['de_json'], '', 'deseq2_volcano_plot.json.gz')

    def test_upload_de_check_field_type(self):
        inputs = {
            'src': 'diff_exp_check_geneid_type.tab.gz',
            'source': 'DICTYBASE',
            'gene_id': 'index',
            'logfc': 'log2FoldChange',
            'fdr': 'padj',
            'pvalue': 'pvalue',
            'stat': 'stat',
        }
        diff_exp = self.run_process('upload-diffexp', inputs)
        saved_json, test_json = self.get_json('diff_exp_check_types.json.gz', diff_exp.output['de_json'])
        self.assertEqual(test_json, saved_json)
        all(self.assertIsInstance(gene, six.text_type) for gene in test_json['gene_id'])

    def test_upload_genome(self):
        inputs = {'src': 'genome.fasta.gz'}
        genome = self.run_process('upload-genome', inputs)
        inputs = {"src": "genome.fasta.gz",
                  "bowtie_index": "bt_index.tar.gz",
                  "bowtie2_index": "bt2_index.tar.gz",
                  "bwa_index": "bwa_index.tar.gz",
                  "hisat2_index": "hisat2_index.tar.gz",
                  "subread_index": "subread_index.tar.gz"}
        genome = self.run_process('upload-genome', inputs)
        self.assertFields(genome, "index_bt", {'dir': 'bowtie_index'})
        self.assertFields(genome, "index_bt2", {'dir': 'bowtie2_index'})
        self.assertFields(genome, "index_bwa", {'dir': 'BWA_index'})
        self.assertFields(genome, "index_hisat2", {'dir': 'hisat2_index'})
        self.assertFields(genome, "index_subread", {'dir': 'subread_index'})

    def test_upload_bed(self):
        inputs = {'src': 'bad.bed'}
        bed = self.run_process('upload-bed', inputs, Data.STATUS_ERROR)

        inputs = {'src': 'good.bed'}
        bed = self.run_process('upload-bed', inputs)
        self.assertFile(bed, 'BED', 'good.bed')

    def test_upload_geneset(self):
        inputs = {'src': 'geneset.tab.gz', 'source': 'UCSC'}
        geneset = self.run_process('upload-geneset', inputs)

        self.assertFile(geneset, 'geneset', 'geneset_out.tab.gz', compression='gzip')
        self.assertFields(geneset, 'source', 'UCSC')
        self.assertJSON(geneset, geneset.output['geneset_json'], '', 'geneset.json.gz')

    def test_create_geneset(self):
        inputs = {'genes': ['ABC', 'DEF', 'GHI'], 'source': 'UCSC'}
        geneset = self.run_process('create-geneset', inputs)

        self.assertFile(geneset, 'geneset', 'geneset_2.tab.gz', compression='gzip')
        self.assertJSON(geneset, geneset.output['geneset_json'], '', 'geneset_2.json.gz')

        inputs = {'genes': ['1', '3', '3', '2'], 'source': 'NCBI'}
        geneset_2 = self.run_process('create-geneset', inputs)

        self.assertFile(geneset_2, 'geneset', 'geneset_3.tab.gz', compression='gzip')
        self.assertJSON(geneset_2, geneset_2.output['geneset_json'], '', 'geneset_3.json.gz')
        self.assertEqual(geneset_2.process_warning[0], 'Removed duplicated genes.')

    def test_create_venn(self):
        inputs = {'genes': ['ABC', 'GHI', 'DEF'], 'source': 'UCSC', 'venn': 'venn.json.gz'}
        venn = self.run_process('create-geneset-venn', inputs)

        self.assertFile(venn, 'geneset', 'geneset_venn.tab.gz', compression='gzip')
        self.assertJSON(venn, venn.output['geneset_json'], '', 'geneset_venn.json.gz')
        self.assertJSON(venn, venn.output['venn'], '', 'venn.json.gz')

    def test_upload_image(self):
        inputs = {'src': 'upload image.1.png'}
        png = self.run_process('upload-image-file', inputs)
        self.assertFile(png, 'file', 'upload image.1.png')

    def test_upload_file(self):
        inputs = {'src': 'upload image.1.png'}
        png = self.run_process('upload-file', inputs)
        self.assertFile(png, 'file', 'upload image.1.png')
        self.assertFields(png, 'file', {'file': 'upload image.1.png', 'size': 4799})

    def test_upload_tabular(self):
        inputs = {'src': 'upload_file.tab.gz'}
        tab = self.run_process('upload-tab-file', inputs)
        self.assertFile(tab, 'file', 'tab_file_tabular.tab.gz', compression='gzip')
        self.assertFile(tab, 'src_file', 'upload_file.tab.gz', compression='gzip')

        inputs = {'src': 'upload_file.csv.gz'}
        csv = self.run_process('upload-tab-file', inputs)
        self.assertFile(csv, 'file', 'csv_file_tabular.tab.gz', compression='gzip')
        self.assertFile(csv, 'src_file', 'upload_file.csv.gz', compression='gzip')

        inputs = {'src': 'upload_file.xls'}
        xls = self.run_process('upload-tab-file', inputs)
        self.assertFile(xls, 'file', 'xls_file_tabular.tab.gz', compression='gzip')
        self.assertFile(xls, 'src_file', 'upload_file.xls')

        inputs = {'src': 'upload file.1.xlsx'}
        xlsx = self.run_process('upload-tab-file', inputs)
        self.assertFile(xlsx, 'file', 'xlsx_file_tabular.tab.gz', compression='gzip')
        self.assertFile(xlsx, 'src_file', 'upload file.1.xlsx')
