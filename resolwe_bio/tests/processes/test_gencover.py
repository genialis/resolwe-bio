# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class GencoverProcessorTestCase(BioProcessTestCase):

    def test_coverage(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff.gz', 'source': 'DICTYBASE'}
        annotation = self.run_process('upload-gff3', inputs)

        # GTF inport
        inputs = {'src': 'annotation_ok.gtf.gz', 'source': 'DICTYBASE'}
        annotation_gtf = self.run_process('upload-gtf', inputs)

        # redundant GTF inport
        inputs = {'src': 'annotation_red.gtf.gz', 'source': 'DICTYBASE'}
        annotation_gtf_red = self.run_process('upload-gtf', inputs)

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_process('alignment-tophat2', inputs)

        # samtools mapping
        inputs = {
            'mapping': aligned_reads.pk,
            'genome': genome.pk}
        variants = self.run_process('vc-samtools', inputs)

        # Coverage report
        inputs = {
            'mapping': aligned_reads.pk,
            'gtf': annotation_gtf.pk,
            'variants': variants.pk,
            'filter': 3,
            'genes': ['geneX']}

        self.run_process('coverage-garvan', inputs)

        # Missing gene in BAM file test
        inputs = {
            'mapping': aligned_reads.pk,
            'gtf': annotation_gtf_red.pk,
            'variants': variants.pk,
            'filter': 3,
            'genes': ['geneX']}

        exon_cov = self.run_process('coverage-garvan', inputs)
        self.assertFile(exon_cov, 'exon_coverage', 'exons_coverage.txt.gz', compression='gzip')
