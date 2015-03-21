# pylint: disable=missing-docstring
from .base import BaseProcessorTestCase
from .utils import PreparedData


class PcaProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_pca(self):
        genome = self.prepare_genome()
        reads1 = self.prepare_reads('00Hr.fastq.gz')
        reads2 = self.prepare_reads('20Hr.fastq.gz')

        inputs = {'src': 'annotation.gtf'}
        annotation = self.run_processor('import:upload:annotation-gtf', inputs)
        self.assertFiles(annotation, 'gtf', 'annotation.gtf')

        inputs = {
            'genome': genome.pk,
            'reads': reads1.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_1 = self.run_processor('alignment:tophat-2-0-13', inputs)

        inputs = {
            'genome': genome.pk,
            'reads': reads2.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_2 = self.run_processor('alignment:tophat-2-0-13', inputs)

        inputs = {
            'alignments': aligned_reads_1.pk,
            'gff': annotation.pk,
            'stranded': "no",
            'id_attribute': 'transcript_id'}
        expression_1 = self.run_processor('htseq-count:-0-6-1p1', inputs)

        inputs = {
            'alignments': aligned_reads_2.pk,
            'gff': annotation.pk,
            'stranded': "no",
            'id_attribute': 'transcript_id'}
        expression_2 = self.run_processor('htseq-count:-0-6-1p1', inputs)

        inputs = {'expressions': [expression_1.pk, expression_2.pk]}
        pca = self.run_processor('pca:1-0-0', inputs)
        self.assertJSON(pca, pca.output['pca'], 'flot.data', 'pca.json')
