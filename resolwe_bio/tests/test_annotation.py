from .base import BaseProcessorTestCase
from server.models import Data


class AnnotationProcessorTestCase(BaseProcessorTestCase):
    def test_transdecoder(self):
        inputs = {'src': 'genome.fasta.gz'}
        genome = self.run_processor('import:upload:genome-fasta', inputs, Data.STATUS_DONE)

        inputs = {'src': '00Hr.fastq.gz'}
        reads1 = self.run_processor('import:upload:reads-fastq', inputs, Data.STATUS_DONE)

        inputs = {'src': '20Hr.fastq.gz'}
        reads2 = self.run_processor('import:upload:reads-fastq', inputs, Data.STATUS_DONE)

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs, Data.STATUS_DONE)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

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
            'alignment': aligned_reads_1.pk,
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_exp = self.run_processor('cufflinks:-2-2-1', inputs, Data.STATUS_DONE)

        inputs = {
            'alignment': aligned_reads_2.pk,
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_exp2 = self.run_processor('cufflinks:-2-2-1', inputs, Data.STATUS_DONE)

        inputs = {
            'expressions': [cuff_exp.pk, cuff_exp2.pk],
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_merge = self.run_processor('cuffmerge:-2-2-1', inputs, Data.STATUS_DONE)

        inputs = {
            'gff': cuff_merge.pk,
            'genome': genome.pk}
        transdecoder = self.run_processor('transdecoder', inputs, Data.STATUS_DONE)
        self.assertFiles(transdecoder, 'gff', 'transdecoder_transcripts.gff')
        self.assertFiles(transdecoder, 'bed', 'transdecoder_transcripts.bed')
