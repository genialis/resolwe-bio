from .base import BaseProcessorTestCase
from .utils import PreparedData
from server.models import Data


class GencoverProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_coverage(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs, Data.STATUS_DONE)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        # GTF inport
        inputs = {'src': 'annotation_ok.gtf'}
        annotation_gtf = self.run_processor('import:upload:annotation-gtf', inputs, Data.STATUS_DONE)
        self.assertFiles(annotation_gtf, 'gtf', 'annotation_ok.gtf')

        # redundant GTF inport
        inputs = {'src': 'annotation_red.gtf'}
        annotation_gtf_red = self.run_processor('import:upload:annotation-gtf', inputs, Data.STATUS_DONE)
        self.assertFiles(annotation_gtf_red, 'gtf', 'annotation_red.gtf')

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs, Data.STATUS_DONE)

        # samtools mapping
        inputs = {
            'mapping': aligned_reads.pk,
            'genome': genome.pk}
        variants = self.run_processor('vc-samtools', inputs, Data.STATUS_DONE)

        # Coverage report
        inputs = {
            'mapping': aligned_reads.pk,
            'gtf': annotation_gtf.pk,
            'variants': variants.pk,
            'filter': 3,
            'genes': ['geneX']}

        genc_results = self.run_processor('coverage:garvan', inputs, Data.STATUS_DONE)
        # self.assertFiles(genc_results, 'bigwig', 'genome_coverage.bw')

        # Missing gene in BAM file test
        inputs = {
            'mapping': aligned_reads.pk,
            'gtf': annotation_gtf_red.pk,
            'variants': variants.pk,
            'filter': 3,
            'genes': ['geneX']}

        genc_results = self.run_processor('coverage:garvan', inputs, Data.STATUS_DONE)
        self.assertFields(genc_results, 'proc.warning', 'Contig scaffold_fake not found in BAM file (for gene geneX)')
