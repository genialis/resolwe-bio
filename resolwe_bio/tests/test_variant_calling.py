# pylint: disable=missing-docstring
from .utils import BioProcessTestCase


class VariantCallingTestCase(BioProcessTestCase):
    def test_variant_calling_samtools(self):
        genome = self.prepare_genome('variant_calling_genome.fasta.gz')
        reads = self.prepare_reads('variant_calling_reads.fastq.gz')

        inputs = {'genome': genome.pk, 'reads': reads.pk, 'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment:bowtie-2-2-3_trim', inputs)

        inputs = {'genome': genome.pk, 'mapping': aligned_reads.pk}
        samtools_variants = self.run_processor('vc-samtools', inputs)
        self.assertFiles(samtools_variants, 'vcf', 'variant_calling_samtools.vcf')

    def test_variant_calling_gatk(self):
        genome = self.prepare_genome('variant_calling_genome.fasta.gz')
        reads = self.prepare_reads('variant_calling_reads.fastq.gz')

        inputs = {'genome': genome.pk, 'reads': reads.pk, 'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment:bowtie-2-2-3_trim', inputs)

        samtools_variants = self.run_processor('import:upload:variants-vcf', {'src': 'variant_calling_samtools.vcf'})

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

    def test_variant_calling_gatk_joint(self):
        genome = self.prepare_genome('variant_calling_genome.fasta.gz')
        reads = self.prepare_reads('variant_calling_reads.fastq.gz')

        inputs = {'genome': genome.pk, 'reads': reads.pk, 'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment:bowtie-2-2-3_trim', inputs)

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
