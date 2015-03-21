# pylint: disable=missing-docstring
from .base import BaseProcessorTestCase
from .utils import PreparedData


class CoverageProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_coverage(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)

        inputs = {'bam': aligned_reads.pk}
        coverage = self.run_processor('bam:coverage', inputs)
        self.assertFiles(coverage, 'bigwig', 'genome_coverage.bw')
