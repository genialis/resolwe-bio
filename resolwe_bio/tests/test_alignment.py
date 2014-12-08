from .base import BaseProcessorTestCase


class AlignmentProcessorTestCase(BaseProcessorTestCase):
    def test_bowtie(self):
        inputs = {'src': 'genome.fasta.gz'}
        genome = self.run_processor('import:upload:genome-fasta', inputs)
        self.assertDone(genome)
        self.assertFiles(genome, 'fasta', 'genome.fasta.gz')

        inputs = {'src': 'reads.fastq.gz'}
        read = self.run_processor('import:upload:reads-fastq', inputs)
        self.assertDone(read)
        self.assertFields(read, 'bases', 35)

        inputs = {'genome': genome.pk, 'reads': read.pk, 'reporting': {'r': "-a -m 1 --best --strata"}}
        aligned_read = self.run_processor('alignment:bowtie-1-0-0-trimmx', inputs)
        self.assertDone(aligned_read)
