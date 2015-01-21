from .base import BaseProcessorTestCase


class CoverageProcessorTestCase(BaseProcessorTestCase):
    def prepair_genome(self):
        inputs = {'src': 'genome.fasta.gz'}
        genome = self.run_processor('import:upload:genome-fasta', inputs)
        self.assertDone(genome)
        self.assertFiles(genome, 'fasta', 'genome.fasta.gz')
        return genome

    def prepair_reads(self):
        inputs = {'src': 'reads.fastq.gz'}
        reads = self.run_processor('import:upload:reads-fastq', inputs)
        self.assertDone(reads)
        self.assertFields(reads, 'bases', 35)
        return reads

    def test_coverage(self):
        genome = self.prepair_genome()
        reads = self.prepair_reads()

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertDone(annotation)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertDone(aligned_reads)

        inputs = {'bam': aligned_reads.pk}
        coverage = self.run_processor('bam:coverage', inputs)
        self.assertDone(coverage)
        self.assertFiles(coverage, 'bigwig', 'genome_coverage.bw')
