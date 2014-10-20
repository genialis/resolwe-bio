from .base import BaseProcessorTestCase


class AlignmentProcessorTestCase(BaseProcessorTestCase):
    def test_bowtie(self):
        input_ = {'src': 'dd_masked_09.05.13.fasta'}
        genome = self.create_data('import:upload:genome-fasta', input_)
        self.assertFiles(genome, 'fasta', 'dd_masked_09.05.13.fasta.gz')

        input_ = {'src': 'AX4_on_ka_00Hr_bio1_small.fq.gz'}
        read = self.create_data('import:upload:reads-fastq', input_)
        self.assertFields(read, 'bases', 35)

        input_ = {'genome': genome.pk, 'reads': read.pk, 'reporting': {'r': "-a -m 1 --best --strata"}}
        aligned_read = self.create_data('alignment:bowtie-1-0-0-trimmx', input_)
