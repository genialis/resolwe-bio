from .base import BaseProcessorTestCase
from .utils import PreparedData
from server.models import Data


class VariantCallingTestCase(BaseProcessorTestCase, PreparedData):
    def test_variant_calling(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads('vc_reads.fastq')

        inputs = {'genome': genome.pk, 'reads': reads.pk, 'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment:bowtie-2-2-3_trim', inputs, Data.STATUS_DONE)
        self.assertFiles(aligned_reads, 'stats', 'VC_bt2_mem_reads_report.txt')

        # samtools variant calling test
        inputs = {'genome': genome.pk, 'mapping': aligned_reads.pk}
        samtools_variants = self.run_processor('vc-samtools', inputs, Data.STATUS_DONE)
        self.assertFiles(samtools_variants, 'vcf', 'VC_reads_align_samtoolscalls.vcf')

        # GATK variant calling test
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
        self.run_processor('vc-gatk', inputs, Data.STATUS_DONE)
        # NOTE: output can not be tested
