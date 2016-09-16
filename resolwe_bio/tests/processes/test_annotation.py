# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class AnnotationProcessorTestCase(BioProcessTestCase):

    def test_transdecoder(self):
        inputs = {'src': ['reads_transdecoder.fastq.gz']}
        reads = self.run_process("upload-fastq-single", inputs)

        inputs = {"src": "genome_transdecoder.fasta.gz"}
        genome = self.run_process('upload-genome', inputs)

        inputs = {'src': 'annotation_transdecoder.gff.gz', 'source': 'dictyBase'}
        annotation = self.run_process('upload-gff3', inputs)

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_process('alignment-tophat2', inputs)

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_exp = self.run_process('cufflinks', inputs)

        inputs = {
            'expressions': [cuff_exp.pk],
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_merge = self.run_process('cuffmerge', inputs)

        inputs = {
            'gff': cuff_merge.pk,
            'genome': genome.pk}
        transdecoder = self.run_process('transdecoder', inputs)
        self.assertFile(transdecoder, 'gff', 'transdecoder_transcripts.gff')
        self.assertFile(transdecoder, 'bed', 'transdecoder_transcripts.bed')
