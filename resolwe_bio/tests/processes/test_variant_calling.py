# pylint: disable=missing-docstring
import unittest
from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase


class VariantCallingTestCase(BioProcessTestCase):

    def test_variant_calling_samtools(self):
        inputs = {"src": "variant_calling_genome.fasta.gz"}
        genome = self.run_processor('upload-genome', inputs)

        reads = self.prepare_reads(['variant_calling_reads.fastq.gz'])

        inputs = {'genome': genome.pk, 'reads': reads.pk, 'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment-bowtie2', inputs)

        inputs = {'genome': genome.pk, 'mapping': aligned_reads.pk}
        samtools_variants = self.run_processor('vc-samtools', inputs)
        # create a filtering function that will remove the samtools version from the output files
        def filter_version(line):
            return line.startswith(b"##samtoolsVersion")
        self.assertFile(samtools_variants, 'vcf', 'variant_calling_samtools.vcf', filter=filter_version)

    @skipDockerFailure("Fails with: int() argument must be a string or a "
        "number, not 'list' at "
        "self.run_processor('vc-gatk', inputs)")
    def test_variant_calling_gatk(self):
        genome = self.prepare_genome('variant_calling_genome.fasta.gz')
        reads = self.prepare_reads('variant_calling_reads.fastq.gz')

        inputs = {'genome': genome.pk, 'reads': reads.pk, 'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment-bowtie2', inputs)

        samtools_variants = self.run_processor('upload-variants-vcf', {'src': 'variant_calling_samtools.vcf'})

        inputs = {
            'genome': genome.pk,
            'mapping': aligned_reads.pk,
            'known_sites': samtools_variants.pk,
            'known_indels': [samtools_variants.pk],
            'reads_info': {
                'ID': "x",
                'SM': "x",
                'PL': "Illumina",
                'LB': "x",
                'CN': "def",
                'DT': "2014-08-05"},
            'Varc_param': {'stand_emit_conf': 10, 'stand_call_conf': 30}}
        self.run_processor('vc-gatk', inputs)
        # NOTE: output can not be tested

    @skipDockerFailure("Fails with: picard: command not found")
    def test_variant_calling_gatk_joint(self):
        genome = self.prepare_genome('variant_calling_genome.fasta.gz')
        reads = self.prepare_reads('variant_calling_reads.fastq.gz')

        inputs = {'genome': genome.pk, 'reads': reads.pk, 'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment-bowtie2', inputs)

        inputs = {
            'genome': genome.pk,
            'mapping': [aligned_reads.pk],
            'reads_info': {
                'PL': "Illumina",
                'LB': "x",
                'CN': "def",
                'DT': "2014-08-05"},
            'Varc_param': {'stand_emit_conf': 10, 'stand_call_conf': 30}}
        self.run_processor('vc-gatk-joint', inputs)
        # NOTE: output can not be tested

    def test_vc_filtering(self):
        variants = self.run_processor('upload-variants-vcf', {'src': 'variant_calling_filtering.vcf.gz'})

        inputs = {
            'variants': variants.pk,
            'analysis_type': 'snv',
            'parental_strain': 'AX4',
            'mutant_strain': 'mutant',
            'read_depth': 5}

        filtered_variants = self.run_processor('chemut', inputs)
        self.assertFile(filtered_variants, 'vcf', 'variant_calling_filtered_variants.vcf')
