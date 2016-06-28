# pylint: disable=missing-docstring
from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase


class AnnotationProcessorTestCase(BioProcessTestCase):

    def test_transdecoder(self):
        inputs = {'src': ['reads_transdecoder.fastq.gz']}
        reads = self.run_processor("upload-fastq-single", inputs)

        inputs = {"src": "genome_transdecoder.fasta.gz"}
        genome = self.run_processor('upload-genome', inputs)

        inputs = {'src': 'annotation_transdecoder.gff.gz'}
        annotation = self.run_processor('upload-gff3', inputs)

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment-tophat2', inputs)

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_exp = self.run_processor('cufflinks', inputs)

        inputs = {
            'expressions': [cuff_exp.pk],
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_merge = self.run_processor('cuffmerge', inputs)

        inputs = {
            'gff': cuff_merge.pk,
            'genome': genome.pk}
        transdecoder = self.run_processor('transdecoder', inputs)
        self.assertFile(transdecoder, 'gff', 'transdecoder_transcripts.gff')
        self.assertFile(transdecoder, 'bed', 'transdecoder_transcripts.bed')
