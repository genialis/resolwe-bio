# pylint: disable=missing-docstring
from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase


class CoverageProcessorTestCase(BioProcessTestCase):

    def test_coverage(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff.gz'}
        annotation = self.run_processor('upload-gff3', inputs)

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment-tophat2', inputs)

        inputs = {'bam': aligned_reads.pk}

        coverage = self.run_processor('bam-coverage', inputs)
        self.assertFile(coverage, 'bigwig', 'genome_coverage.bw')
