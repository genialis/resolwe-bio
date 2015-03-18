from .base import BaseProcessorTestCase
from server.models import Data


class CoverageProcessorTestCase(BaseProcessorTestCase):
    def prepair_genome(self):
        inputs = {'src': 'genome.fasta.gz'}
        genome = self.run_processor('import:upload:genome-fasta', inputs, Data.STATUS_DONE)
        self.assertFiles(genome, 'fasta', 'genome.fasta.gz')
        return genome

    def prepair_reads(self):
        inputs = {'src': 'reads.fastq.gz'}
        reads = self.run_processor('import:upload:reads-fastq', inputs, Data.STATUS_DONE)
        self.assertFields(reads, 'bases', 35)
        return reads

    def test_coverage(self):
        genome = self.prepair_genome()
        reads = self.prepair_reads()

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs, Data.STATUS_DONE)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs, Data.STATUS_DONE)

        inputs = {'bam': aligned_reads.pk}
        coverage = self.run_processor('bam:coverage', inputs, Data.STATUS_DONE)
        self.assertFiles(coverage, 'bigwig', 'genome_coverage.bw')
