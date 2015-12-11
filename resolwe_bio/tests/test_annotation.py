# pylint: disable=missing-docstring
from .utils import BioProcessTestCase


class AnnotationProcessorTestCase(BioProcessTestCase):
    def test_transdecoder(self):
        inputs = {'src': 'reads_transdecoder.fastq.gz'}
        reads = self.run_processor("import:upload:reads-fastq", inputs)

        inputs = {"src": "genome_transdecoder.fasta.gz"}
        genome = self.run_processor('import:upload:genome-fasta', inputs)

        inputs = {'src': 'annotation_transdecoder.gff.gz'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_exp = self.run_processor('cufflinks:-2-2-1', inputs)

        inputs = {
            'expressions': [cuff_exp.pk],
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_merge = self.run_processor('cuffmerge:-2-2-1', inputs)

        inputs = {
            'gff': cuff_merge.pk,
            'genome': genome.pk}
        transdecoder = self.run_processor('transdecoder', inputs)
        self.assertFiles(transdecoder, 'gff', 'transdecoder_transcripts.gff')
        self.assertFiles(transdecoder, 'bed', 'transdecoder_transcripts.bed')
