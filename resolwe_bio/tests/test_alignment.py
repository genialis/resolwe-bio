from .base import BaseProcessorTestCase


class AlignmentProcessorTestCase(BaseProcessorTestCase):
    def test_bowtie(self):
        inputs = {'src': 'dd_masked_09.05.13.fasta'}
        genome = self.run_processor('import:upload:genome-fasta', inputs)
        self.assertFiles(genome, 'fasta', 'dd_masked_09.05.13.fasta.gz')

        inputs = {'src': 'AX4_on_ka_00Hr_bio1_small.fq.gz'}
        read = self.run_processor('import:upload:reads-fastq', inputs)
        self.assertFields(read, 'bases', 35)

        inputs = {'genome': genome.pk, 'reads': read.pk, 'reporting': {'r': "-a -m 1 --best --strata"}}
        aligned_read = self.run_processor('alignment:bowtie-1-0-0-trimmx', inputs)
