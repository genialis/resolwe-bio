# pylint: disable=missing-docstring
from .base import BaseProcessorTestCase
from .utils import PreparedData


class GencoverProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_coverage(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        # GTF inport
        inputs = {'src': 'annotation_ok.gtf'}
        annotation_gtf = self.run_processor('import:upload:annotation-gtf', inputs)
        self.assertFiles(annotation_gtf, 'gtf', 'annotation_ok.gtf')

        # redundant GTF inport
        inputs = {'src': 'annotation_red.gtf'}
        annotation_gtf_red = self.run_processor('import:upload:annotation-gtf', inputs)
        self.assertFiles(annotation_gtf_red, 'gtf', 'annotation_red.gtf')

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)

        # samtools mapping
        inputs = {
            'mapping': aligned_reads.pk,
            'genome': genome.pk}
        variants = self.run_processor('vc-samtools', inputs)

        # Coverage report
        inputs = {
            'mapping': aligned_reads.pk,
            'gtf': annotation_gtf.pk,
            'variants': variants.pk,
            'filter': 3,
            'genes': ['geneX']}

        self.run_processor('coverage:garvan', inputs)

        # Missing gene in BAM file test
        inputs = {
            'mapping': aligned_reads.pk,
            'gtf': annotation_gtf_red.pk,
            'variants': variants.pk,
            'filter': 3,
            'genes': ['geneX']}

        self.run_processor('coverage:garvan', inputs)
