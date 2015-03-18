from .base import BaseProcessorTestCase
from server.models import Data


class PcaProcessorTestCase(BaseProcessorTestCase):
    def prepair_genome(self):
        inputs = {'src': 'genome.fasta.gz'}
        genome = self.run_processor('import:upload:genome-fasta', inputs, Data.STATUS_DONE)
        self.assertFiles(genome, 'fasta', 'genome.fasta.gz')
        return genome

    def test_pca(self):
        genome = self.prepair_genome()

        inputs = {'src': '00Hr.fastq.gz'}
        reads1 = self.run_processor('import:upload:reads-fastq', inputs, Data.STATUS_DONE)

        inputs = {'src': '20Hr.fastq.gz'}
        reads2 = self.run_processor('import:upload:reads-fastq', inputs, Data.STATUS_DONE)

        inputs = {'src': 'annotation.gtf'}
        annotation = self.run_processor('import:upload:annotation-gtf', inputs, Data.STATUS_DONE)
        self.assertFiles(annotation, 'gtf', 'annotation.gtf')

        inputs = {
            'genome': genome.pk,
            'reads': reads1.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_1 = self.run_processor('alignment:tophat-2-0-13', inputs, Data.STATUS_DONE)

        inputs = {
            'genome': genome.pk,
            'reads': reads2.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_2 = self.run_processor('alignment:tophat-2-0-13', inputs, Data.STATUS_DONE)

        inputs = {
            'alignments': aligned_reads_1.pk,
            'gff': annotation.pk,
            'stranded': "no",
            'id_attribute': 'transcript_id'}
        expression_1 = self.run_processor('htseq-count:-0-6-1p1', inputs, Data.STATUS_DONE)

        inputs = {
            'alignments': aligned_reads_2.pk,
            'gff': annotation.pk,
            'stranded': "no",
            'id_attribute': 'transcript_id'}
        expression_2 = self.run_processor('htseq-count:-0-6-1p1', inputs, Data.STATUS_DONE)

        inputs = {'expressions': [expression_1.pk, expression_2.pk]}
        pca = self.run_processor('pca:1-0-0', inputs, Data.STATUS_DONE)
        self.assertJSON(pca, pca.output['pca'], 'flot.data', 'pca.json')
