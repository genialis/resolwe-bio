# pylint: disable=missing-docstring
import six

from django.core.exceptions import ValidationError

from resolwe.flow.models import Data, Secret
from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase


class UploadProcessorTestCase(BioProcessTestCase):

    @tag_process('upload-bam', 'upload-bam-indexed')
    def test_bam_upload(self):
        inputs = {
            'src': 'alignment_name_sorted.bam',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        upload_bam = self.run_process('upload-bam', inputs)
        self.assertFile(upload_bam, 'bam', 'alignment_position_sorted.bam')
        self.assertFile(upload_bam, 'bai', 'alignment_bam_upload_index.bai')
        self.assertFields(upload_bam, 'species', 'Homo sapiens')
        self.assertFields(upload_bam, 'build', 'hg19')

        inputs = {
            'src': 'alignment_position_sorted.bam',
            'src2': 'alignment_bam_upload_index.bam.bai',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        upload_bam = self.run_process('upload-bam-indexed', inputs, Data.STATUS_ERROR)
        self.assertEqual(upload_bam.process_error[0], 'BAI should have the same name as BAM with .bai extension')

        inputs = {
            'src': 'alignment_position_sorted.bam',
            'src2': 'alignment_position_sorted.bam.bai',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        upload_bam = self.run_process('upload-bam-indexed', inputs)
        self.assertFile(upload_bam, 'bam', 'alignment_position_sorted.bam')
        self.assertFile(upload_bam, 'bai', 'alignment_position_sorted.bam.bai')
        self.assertFields(upload_bam, 'species', 'Homo sapiens')
        self.assertFields(upload_bam, 'build', 'hg19')

    @tag_process('upload-expression')
    def test_upload_expression(self):
        inputs = {
            'exp_type': 'TPM',
            'exp_name': 'Expression',
            'source': 'UCSC',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        self.run_process('upload-expression', inputs, Data.STATUS_ERROR)

        inputs = {
            'exp': 'exp_1_tpm.tab.gz',
            'rc': 'exp_1_rc.tab.gz',
            'exp_name': 'Expression',
            'source': 'UCSC',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        self.run_process('upload-expression', inputs, Data.STATUS_ERROR)

        inputs = {
            'rc': 'exp_1_rc.tab.gz',
            'exp_name': 'Expression',
            'source': 'UCSC',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        exp_3 = self.run_process('upload-expression', inputs)
        self.assertFile(exp_3, 'rc', 'exp_1_rc.tab.gz')
        self.assertFile(exp_3, 'exp', 'exp_1_rc.tab.gz')
        self.assertJSON(exp_3, exp_3.output['exp_json'], '', 'exp_1.json.gz')
        self.assertFields(exp_3, 'species', 'Homo sapiens')
        self.assertFields(exp_3, 'build', 'hg19')
        self.assertFields(exp_3, 'feature_type', 'gene')

        inputs = {
            'exp': 'exp_1_tpm.tab.gz',
            'exp_type': 'TPM',
            'exp_name': 'Expression',
            'source': 'UCSC',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        exp_4 = self.run_process('upload-expression', inputs)
        self.assertFile(exp_4, 'exp', 'exp_1_tpm.tab.gz')

        inputs = {
            'rc': 'exp_1_rc.tab.gz',
            'exp': 'exp_1_tpm.tab.gz',
            'exp_type': 'TPM',
            'exp_name': 'Expression',
            'source': 'UCSC',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        exp_5 = self.run_process('upload-expression', inputs)
        self.assertFields(exp_5, 'exp_type', 'TPM')
        self.assertFile(exp_5, 'exp', 'exp_1_tpm.tab.gz')
        self.assertFile(exp_5, 'rc', 'exp_1_rc.tab.gz')
        self.assertJSON(exp_5, exp_5.output['exp_json'], '', 'exp_1_norm.json.gz')

        inputs = {
            'rc': 'exp_mac_line_ending.txt.gz',
            'exp_name': 'Expression',
            'source': 'UCSC',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        exp_6 = self.run_process('upload-expression', inputs)
        self.assertJSON(exp_6, exp_6.output['exp_json'], '', 'exp.json.gz')

        inputs = {
            'rc': 'exp_unix_line_ending.txt.gz',
            'exp_name': 'Expression',
            'source': 'UCSC',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        exp_7 = self.run_process('upload-expression', inputs)
        self.assertJSON(exp_7, exp_7.output['exp_json'], '', 'exp.json.gz')

        inputs = {
            'rc': 'exp_windows_line_ending.txt.gz',
            'exp_name': 'Expression',
            'source': 'UCSC',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        exp_8 = self.run_process('upload-expression', inputs)
        self.assertJSON(exp_8, exp_8.output['exp_json'], '', 'exp.json.gz')

    @tag_process('upload-cxb', 'upload-expression-cuffnorm')
    def test_upload_cuffquant_expr(self):
        inputs = {
            'src': 'cuffquant_1.cxb',
            'source': 'UCSC',
            'species': 'Homo sapiens',
            'build': 'hg19'}
        cxb = self.run_process('upload-cxb', inputs)

        inputs = {
            'exp': 'cuffquant_exp.tab',
            'cxb': cxb.id}
        exp = self.run_process('upload-expression-cuffnorm', inputs)
        self.assertFields(exp, 'feature_type', 'gene')

    @tag_process('upload-fastq-paired')
    def test_upload_paired_end_reads(self):
        inputs = {'src1': ['mate1.fastq.gz'], 'src2': ['mate2.fastq.gz']}
        self.run_process('upload-fastq-paired', inputs, Data.STATUS_ERROR)

        inputs = {'src1': ['rRNA forw.fastq.gz', 'rRNA_rew.fastq.gz'],
                  'src2': ['00Hr.fastq.gz', '20Hr.fastq.gz']}

        reads = self.run_process('upload-fastq-paired', inputs)
        self.assertFiles(reads, 'fastq', ['rRNA forw.fastq.gz', 'rRNA_rew.fastq.gz'], compression='gzip')
        self.assertFiles(reads, 'fastq2', ['00Hr.fastq.gz', '20Hr.fastq.gz'], compression='gzip')
        del reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        del reads.output['fastqc_url'][1]['total_size']  # Non-deterministic output.
        self.assertFields(reads, 'fastqc_url', [{'file': 'fastqc/rRNA forw_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/rRNA forw_fastqc'],
                                                 'size': 343222},
                                                {'file': 'fastqc/rRNA_rew_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/rRNA_rew_fastqc'],
                                                 'size': 323297}])
        del reads.output['fastqc_url2'][0]['total_size']  # Non-deterministic output.
        del reads.output['fastqc_url2'][1]['total_size']  # Non-deterministic output.
        self.assertFields(reads, 'fastqc_url2', [{'file': 'fastqc/00Hr_fastqc/fastqc_report.html',
                                                  'refs': ['fastqc/00Hr_fastqc'],
                                                  'size': 327878},
                                                 {'file': 'fastqc/20Hr_fastqc/fastqc_report.html',
                                                  'refs': ['fastqc/20Hr_fastqc'],
                                                  'size': 287245}])

    @tag_process('upload-fastq-single')
    def test_upload_single_end_reads(self):
        inputs = {'src': ['mate1.fastq.gz']}
        self.run_process('upload-fastq-single', inputs, Data.STATUS_ERROR)

        inputs = {'src': ['rRNA forw.fastq.gz', 'rRNA_rew.fastq.gz']}
        reads = self.run_process('upload-fastq-single', inputs)

        self.assertFiles(reads, 'fastq', ['rRNA_forw_single.fastq.gz', 'rRNA_rew.fastq.gz'], compression='gzip')
        del reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        del reads.output['fastqc_url'][1]['total_size']  # Non-deterministic output.
        self.assertFields(reads, 'fastqc_url', [{'file': 'fastqc/rRNA forw_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/rRNA forw_fastqc'],
                                                 'size': 343222},
                                                {'file': 'fastqc/rRNA_rew_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/rRNA_rew_fastqc'],
                                                 'size': 323297}])

    @tag_process('upload-diffexp')
    def test_upload_de(self):
        inputs = {
            'src': 'deseq2_output.tab.gz',
            'source': 'DICTYBASE',
            'gene_id': 'gene_id',
            'logfc': 'log2FoldChange',
            'fdr': 'padj',
            'pvalue': 'pvalue',
            'stat': 'stat',
            'species': 'Dictyostelium discoideum',
            'build': 'dd-05-2009',
            'feature_type': 'gene'
        }
        diff_exp = self.run_process('upload-diffexp', inputs)

        self.assertFile(diff_exp, 'raw', 'deseq2_output.tab.gz')
        self.assertJSON(diff_exp, diff_exp.output['de_json'], '', 'deseq2_volcano_plot.json.gz')

    @tag_process('upload-diffexp')
    def test_upload_de_check_field_type(self):
        inputs = {
            'src': 'diff_exp_check_geneid_type.tab.gz',
            'source': 'DICTYBASE',
            'gene_id': 'index',
            'logfc': 'log2FoldChange',
            'fdr': 'padj',
            'pvalue': 'pvalue',
            'stat': 'stat',
            'species': 'Dictyostelium discoideum',
            'build': 'dd-05-2009',
            'feature_type': 'gene'
        }
        diff_exp = self.run_process('upload-diffexp', inputs)
        saved_json, test_json = self.get_json('diff_exp_check_types.json.gz', diff_exp.output['de_json'])
        self.assertEqual(test_json, saved_json)
        all(self.assertIsInstance(gene, six.text_type) for gene in test_json['gene_id'])

    @tag_process('upload-genome')
    def test_upload_genome(self):
        inputs = {
            'src': 'genome.fasta.gz',
            'species': 'Dictyostelium discoideum',
            'build': 'dd-05-2009'
        }
        genome = self.run_process('upload-genome', inputs)

        inputs = {
            'src': 'genome.fasta.gz',
            'species': 'Dictyostelium discoideum',
            'build': 'dd-05-2009',
            'advanced': {
                'bowtie_index': 'bt_index.tar.gz',
                'bowtie2_index': 'bt2_index.tar.gz',
                'bwa_index': 'bwa_index.tar.gz',
                'hisat2_index': 'hisat2_index.tar.gz',
                'subread_index': 'subread_index.tar.gz'
            }
        }
        genome = self.run_process('upload-genome', inputs)
        del genome.output['index_bt']['total_size']  # Non-deterministic output.
        self.assertFields(genome, "index_bt", {'dir': 'bowtie_index'})
        del genome.output['index_bt2']['total_size']  # Non-deterministic output.
        self.assertFields(genome, "index_bt2", {'dir': 'bowtie2_index'})
        del genome.output['index_bwa']['total_size']  # Non-deterministic output.
        self.assertFields(genome, "index_bwa", {'dir': 'BWA_index'})
        del genome.output['index_hisat2']['total_size']  # Non-deterministic output.
        self.assertFields(genome, "index_hisat2", {'dir': 'hisat2_index'})
        del genome.output['index_subread']['total_size']  # Non-deterministic output.
        self.assertFields(genome, "index_subread", {'dir': 'subread_index'})

    @tag_process('upload-bed')
    def test_upload_bed(self):
        inputs = {
            'src': 'bad.bed',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        bed = self.run_process('upload-bed', inputs, Data.STATUS_ERROR)

        inputs = {
            'src': 'good.bed',
            'species': 'Homo sapiens',
            'build': 'hg19'
        }
        bed = self.run_process('upload-bed', inputs)
        self.assertFile(bed, 'bed', 'good.bed')

    @tag_process('upload-geneset')
    def test_upload_geneset(self):
        inputs = {
            'src': 'geneset.tab.gz',
            'source': 'UCSC',
            'species': 'Homo sapiens'
        }
        geneset = self.run_process('upload-geneset', inputs)

        self.assertFile(geneset, 'geneset', 'geneset_out.tab.gz', compression='gzip')
        self.assertFields(geneset, 'source', 'UCSC')
        self.assertFields(geneset, 'species', 'Homo sapiens')
        self.assertJSON(geneset, geneset.output['geneset_json'], '', 'geneset.json.gz')

    @tag_process('create-geneset')
    def test_create_geneset(self):
        inputs = {
            'genes': ['ABC', 'DEF', 'GHI'],
            'source': 'UCSC',
            'species': 'Homo sapiens',
        }
        geneset = self.run_process('create-geneset', inputs)

        self.assertFile(geneset, 'geneset', 'geneset_2.tab.gz', compression='gzip')
        self.assertJSON(geneset, geneset.output['geneset_json'], '', 'geneset_2.json.gz')
        self.assertFields(geneset, 'source', 'UCSC')
        self.assertFields(geneset, 'species', 'Homo sapiens')

        inputs = {
            'genes': ['1', '3', '3', '2'],
            'source': 'NCBI',
            'species': 'Homo sapiens',
        }
        geneset_2 = self.run_process('create-geneset', inputs)

        self.assertFile(geneset_2, 'geneset', 'geneset_3.tab.gz', compression='gzip')
        self.assertJSON(geneset_2, geneset_2.output['geneset_json'], '', 'geneset_3.json.gz')
        self.assertEqual(geneset_2.process_warning[0], 'Removed duplicated genes.')

    @tag_process('create-geneset-venn')
    def test_create_venn(self):
        inputs = {
            'genes': ['ABC', 'GHI', 'DEF'],
            'source': 'UCSC',
            'venn': 'venn.json.gz',
            'species': 'Homo sapiens',
        }
        venn = self.run_process('create-geneset-venn', inputs)

        self.assertFile(venn, 'geneset', 'geneset_venn.tab.gz', compression='gzip')
        self.assertJSON(venn, venn.output['geneset_json'], '', 'geneset_venn.json.gz')
        self.assertJSON(venn, venn.output['venn'], '', 'venn.json.gz')
        self.assertFields(venn, 'source', 'UCSC')
        self.assertFields(venn, 'species', 'Homo sapiens')

    @tag_process('upload-fastq-single')
    def test_upload_reformating_single(self):
        inputs = {'src': ['old_encoding.fastq.gz']}
        reads = self.run_process('upload-fastq-single', inputs)
        self.assertFiles(reads, 'fastq', ['old_encoding_transformed.fastq.gz'], compression='gzip')

    @tag_process('upload-fastq-paired')
    def test_upload_reformating_paired(self):
        inputs = {'src1': ['old_encoding.fastq.gz', 'old_encoding1.fastq.gz'],
                  'src2': ['old_encoding_R2.fastq.gz', 'old_encoding1_R2.fastq.gz']}
        reads = self.run_process('upload-fastq-paired', inputs)
        self.assertFiles(reads, 'fastq', ['old_encoding_transformed.fastq.gz',
                                          'old_encoding1_transformed.fastq.gz'], compression='gzip')
        self.assertFiles(reads, 'fastq2', ['old_encoding_transformed_R2.fastq.gz',
                                           'old_encoding1_transformed_R2.fastq.gz'], compression='gzip')

    @tag_process('upload-master-file')
    def test_upload_master_file(self):
        inputs = {
            'src': '56G_masterfile_corrupted.txt',
            'panel_name': '56G panel, v2'
        }
        master_file = self.run_process('upload-master-file', inputs, Data.STATUS_ERROR)

        inputs = {
            'src': 'amplicon_master_file_merged.bed',
            'panel_name': '56G panel, v2'
        }

        with self.assertRaises(ValidationError):
            self.run_process('upload-master-file', inputs)

        inputs = {
            'src': '56G_masterfile_170113.txt.gz',
            'panel_name': '56G panel, v2'
        }
        master_file = self.run_process('upload-master-file', inputs)

        self.assertFile(master_file, 'bedfile', 'amplicon_master_file_merged.bed')
        self.assertFile(master_file, 'nomergebed', 'amplicon_master_file_nomergebed.bed')
        self.assertFile(master_file, 'olapfreebed', 'amplicon_master_file_olapfreebed.bed')
        self.assertFile(master_file, 'primers', 'amplicon_primers.bed')
        self.assertFields(master_file, 'panel_name', '56G panel, v2')

    @tag_process('upload-etc')
    def test_upload_etc(self):
        inputs = {'src': 'etc_upload_input.xls'}
        etc = self.run_process('upload-etc', inputs)

        self.assertFile(etc, 'etcfile', 'test_etc.json.gz')

    @tag_process('upload-fasta-nucl')
    def test_upload_nucl_seq(self):
        inputs = {
            'src': 'genome.fasta.gz',
            'species': 'Dictyostelium discoideum',
        }
        seq = self.run_process('upload-fasta-nucl', inputs)
        self.assertFile(seq, 'fasta', 'genome.fasta.gz', compression='gzip')
        self.assertFields(seq, 'species', 'Dictyostelium discoideum')

    @tag_process('basespace-fastq-single')
    def test_basespace_fastq_single(self):
        # Token with limited scope preobtained from dedicated BaseSpace testing app.
        handle = Secret.objects.create_secret('9bdf059c759a429f8af52ca084130060', self.admin)

        inputs = {'file_ids': ['9461130722', '9461121664'], 'access_token_secret': {'handle': handle}}
        reads = self.run_process('basespace-fastq-single', inputs)

        self.assertFiles(reads, 'fastq', ['Test_S1_L001_R1_001.fastq.gz',
                                          'Test_S1_L002_R1_001.fastq.gz'], compression='gzip')
        del reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        del reads.output['fastqc_url'][1]['total_size']  # Non-deterministic output.
        self.assertFields(reads, 'fastqc_url', [{'file': 'fastqc/Test_S1_L001_R1_001_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/Test_S1_L001_R1_001_fastqc'],
                                                 'size': 343222},
                                                {'file': 'fastqc/Test_S1_L002_R1_001_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/Test_S1_L002_R1_001_fastqc'],
                                                 'size': 343222}])

    @tag_process('basespace-fastq-paired')
    def test_basespace_fastq_paired(self):
        # Token with limited scope preobtained from dedicated BaseSpace testing app.
        handle = Secret.objects.create_secret('d0728b8cceb7455786665453d28c7ebc', self.admin)

        inputs = {
            'upstream_file_ids': ['9864012319', '9863993106'],
            'downstream_file_ids': ['9863999826', '9863993107'],
            'access_token_secret': {'handle': handle}
        }
        reads = self.run_process('basespace-fastq-paired', inputs)

        self.assertFiles(reads, 'fastq', ['TestPaired_S1_L001_R1_001.fastq.gz',
                                          'TestPaired_S1_L002_R1_001.fastq.gz'], compression='gzip')
        self.assertFiles(reads, 'fastq2', ['TestPaired_S1_L001_R2_001.fastq.gz',
                                           'TestPaired_S1_L002_R2_001.fastq.gz'], compression='gzip')
        del reads.output['fastqc_url'][0]['total_size']  # Non-deterministic output.
        del reads.output['fastqc_url'][1]['total_size']  # Non-deterministic output.
        self.assertFields(reads, 'fastqc_url', [{'file': 'fastqc/TestPaired_S1_L001_R1_001_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/TestPaired_S1_L001_R1_001_fastqc'],
                                                 'size': 347881},
                                                {'file': 'fastqc/TestPaired_S1_L002_R1_001_fastqc/fastqc_report.html',
                                                 'refs': ['fastqc/TestPaired_S1_L002_R1_001_fastqc'],
                                                 'size': 347881}])
        del reads.output['fastqc_url2'][0]['total_size']  # Non-deterministic output.
        del reads.output['fastqc_url2'][1]['total_size']  # Non-deterministic output.
        self.assertFields(reads, 'fastqc_url2', [{'file': 'fastqc/TestPaired_S1_L001_R2_001_fastqc/fastqc_report.html',
                                                  'refs': ['fastqc/TestPaired_S1_L001_R2_001_fastqc'],
                                                  'size': 351553},
                                                 {'file': 'fastqc/TestPaired_S1_L002_R2_001_fastqc/fastqc_report.html',
                                                  'refs': ['fastqc/TestPaired_S1_L002_R2_001_fastqc'],
                                                  'size': 351553}])
