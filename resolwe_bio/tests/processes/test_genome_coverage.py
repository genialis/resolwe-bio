# pylint: disable=missing-docstring
from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase


class CoverageProcessorTestCase(BioProcessTestCase):

    @tag_process('bam-coverage')
    def test_coverage(self):
        with self.preparation_stage():
            genome = self.prepare_genome()
            reads = self.prepare_reads()
            annotation = self.prepare_annotation_gff()

            inputs = {
                'genome': genome.id,
                'reads': reads.id,
                'annotation': annotation.id,
                'PE_options': {
                    'library_type': "fr-unstranded"}}
            aligned_reads = self.run_process('alignment-tophat2', inputs)

        inputs = {'bam': aligned_reads.id}

        coverage = self.run_process('bam-coverage', inputs)
        self.assertFile(coverage, 'bigwig', 'genome_coverage.bw')
