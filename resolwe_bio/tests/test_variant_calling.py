# pylint: disable=missing-docstring
import unittest
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

    @unittest.skip("Missing tools in runtime")
    def test_vc_filtering(self):
        variants = self.run_processor('import:upload:variants-vcf', {'src': 'variant_calling_filtering.vcf.gz'})

        inputs = {
            'variants': variants.pk,
            'analysis_type': 'snv',
            'parental_strain': 'AX4',
            'mutant_strain': 'mutant',
            'read_depth': 5}

        filtered_variants = self.run_processor('vc_filtering_chem_mutagenesis', inputs)
        self.assertFiles(filtered_variants, 'vcf', 'variant_calling_filtered_variants.vcf')

    @unittest.skip("Run to prepare VC test project data")
    def test_vc_prepare_project(self):
        # upload BAM files and prepare JBrowse Coverage tracks
        mutant_1 = self.run_processor("import:upload:mapping-bam", {"src": "mutant_1.bam"})
        mutant_1_coverage = self.run_processor('jbrowse:bam:coverage', {'bam': mutant_1.pk})

        mutant_2 = self.run_processor("import:upload:mapping-bam", {"src": "mutant_2.bam"})
        mutant_2_coverage = self.run_processor('jbrowse:bam:coverage', {'bam': mutant_2.pk})

        AX4 = self.run_processor("import:upload:mapping-bam", {"src": "AX4.bam"})
        AX4_coverage = self.run_processor('jbrowse:bam:coverage', {'bam': AX4.pk})

        # upload genome file and prepare JBrowse refseq track
        inputs = {"src": "dd_masked_05-13-2009.fasta.gz"}
        genome = self.run_processor('import:upload:genome-fasta', inputs)
        genome_track = self.run_processor('jbrowse:refseq', {'refseq': genome.pk})

        # upload annotation file and prepare JBrowse annotation track
        inputs = {'src': 'dd_07-15-2014.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        gff_track = self.run_processor('jbrowse:gff3', {'gff': annotation.pk})

        # Variant calling - SNVs
        inputs = {
            'genome': genome.pk,
            'mapping': [AX4.pk, mutant_1.pk, mutant_2.pk],
            'reads_info': {
                'PL': "Illumina",
                'LB': "x",
                'CN': "def",
                'DT': "2016-01-25"},
            'Varc_param': {'stand_emit_conf': 10, 'stand_call_conf': 30, 'ploidy': 1, 'glm': 'SNP'}}
        snv = self.run_processor('vc-gatk-joint', inputs)

        # Variant calling - INDELS
        inputs = {
            'genome': genome.pk,
            'mapping': [AX4.pk, mutant_1.pk, mutant_2.pk],
            'reads_info': {
                'PL': "Illumina",
                'LB': "x",
                'CN': "def",
                'DT': "2016-01-25"},
            'Varc_param': {'stand_emit_conf': 10, 'stand_call_conf': 30, 'ploidy': 1, 'glm': 'INDEL'}}
        indel = self.run_processor('vc-gatk-joint', inputs)

        # Filter VC results files
        inputs = {
            'variants': snv.pk,
            'analysis_type': 'snv',
            'parental_strain': 'AX4',
            'mutant_strain': 'mutant',
            'read_depth': 5}

        snv_filtered = self.run_processor('vc_filtering_chem_mutagenesis', inputs)

        inputs = {
            'variants': indel.pk,
            'analysis_type': 'indel',
            'parental_strain': 'AX4',
            'mutant_strain': 'mutant',
            'read_depth': 5}

        indel_filtered = self.run_processor('vc_filtering_chem_mutagenesis', inputs)
