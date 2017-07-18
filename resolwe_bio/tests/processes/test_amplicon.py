# pylint: disable=missing-docstring
from os.path import join

from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles


class AmpliconProcessorTestCase(BioProcessTestCase):

    def test_bwa_trim(self):
        inputs = {
            'src1': ['56GSID_10k_mate1.fastq.gz'],
            'src2': ['56GSID_10k_mate2.fastq.gz']}
        reads = self.run_process('upload-fastq-paired', inputs)

        genome = self.run_process('upload-genome', {'src': 'hs_b37_chr2_small.fasta.gz'})
        master_file = self.run_process('upload-master-file', {'src': '56G_masterfile_test.txt'})

        inputs = {
            'master_file': master_file.id,
            'genome': genome.id,
            'reads': reads.id
        }
        bwa_trim = self.run_process('align-bwa-trim', inputs)

        self.assertFile(bwa_trim, 'stats', 'bwa_trim_stats.txt')

    @skipUnlessLargeFiles('56GSID_10k_mate1_RG.bam')
    def test_amplicon_table(self):
        bam = self.run_process('upload-bam', {'src': join('large', '56GSID_10k_mate1_RG.bam')})
        master_file = self.run_process('upload-master-file', {'src': '56G_masterfile_test.txt'})
        coverage = self.run_process('coveragebed', {'alignment': bam.id, 'master_file': master_file.id})

        inputs = {
            'annotation': '56GSID.lf.finalvars.txt',
            'summary': '56GSID_1k.gatkHC_snpEff_summary.html',
            'snpeff_genes': '56GSID_1k.gatkHC_snpEff_genes.txt'
        }

        annot_variants = self.run_process('upload-snpeff', inputs)

        amplicon_table_inputs = {
            'master_file': master_file.id,
            'coverage': coverage.id,
            'annot_vars': [annot_variants.id]
        }

        table = self.run_process('amplicon-table', amplicon_table_inputs)
        self.assertJSON(table, table.output['variant_table'], '', 'amplicon_table_output.json.gz')

        amplicon_table_inputs = {
            'master_file': master_file.id,
            'coverage': coverage.id,
            'annot_vars': [annot_variants.id],
            'all_amplicons': True,
            'table_name': 'All amplicons'
        }

        table = self.run_process('amplicon-table', amplicon_table_inputs)
        self.assertJSON(table, table.output['variant_table'], '', 'amplicon_table_output_all.json.gz')
