from .base import BaseProcessorTestCase
from .utils import PreparedData
from server.models import Data


class DiffExpProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_cuffdiff(self):
        genome = self.prepare_genome()
        reads1 = self.prepare_reads('00Hr.fastq.gz')
        reads2 = self.prepare_reads('20Hr.fastq.gz')

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
            'alignment': aligned_reads_1.pk,
            'gff': cuff_merge.pk}
        cuffquant = self.run_processor('cuffquant:-2-2-1', inputs, Data.STATUS_DONE)

        inputs = {
            'alignment': aligned_reads_2.pk,
            'gff': cuff_merge.pk}
        cuffquant2 = self.run_processor('cuffquant:-2-2-1', inputs, Data.STATUS_DONE)

        inputs = {
            'cuffquant': [cuffquant.pk, cuffquant2.pk],
            'replicates': ['1', '2'],
            'labels': ['g1', 'g2'],
            'gff': cuff_merge.pk}
        cuffdiff = self.run_processor('cuffdiff:-2-2-1', inputs, Data.STATUS_DONE)
        self.assertFiles(cuffdiff, 'gene_diff_exp', 'cuffdiff_output.gz', gzipped=True)

    def test_bayseq_bcm(self):
        genome = self.prepare_genome()
        reads1 = self.prepare_reads('00Hr.fastq.gz')
        reads2 = self.prepare_reads('20Hr.fastq.gz')

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
            'genome': genome.pk,
            'gff': annotation.pk}
        mappability = self.run_processor('mappability:bcm-1-0-0', inputs, Data.STATUS_DONE)

        inputs = {
            'alignment': aligned_reads_1.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression_1 = self.run_processor('expression:bcm-1-0-0', inputs, Data.STATUS_DONE)

        inputs = {
            'alignment': aligned_reads_2.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression_2 = self.run_processor('expression:bcm-1-0-0', inputs, Data.STATUS_DONE)

        inputs = {
            'name': "00vs20",
            'case': [expression_1.pk],
            'control': [expression_2.pk],
            'replicates': ['1', '2'],
            'mappability': mappability.pk}
        diff_exp = self.run_processor('differentialexpression:bcm-1-0-0', inputs, Data.STATUS_DONE)
        self.assertJSON(diff_exp, diff_exp.output['volcano_plot'], '', 'bayseq_volcano.json')

    def test_deseq2(self):
        genome = self.prepare_genome()
        reads1 = self.prepare_reads('00Hr.fastq.gz')
        reads2 = self.prepare_reads('20Hr.fastq.gz')

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

        inputs = {
            'name': "00vs20",
            'case': [expression_1.pk],
            'control': [expression_2.pk]}
        diff_exp = self.run_processor('differentialexpression:deseq2', inputs, Data.STATUS_DONE)
        self.assertFiles(diff_exp, "diffexp", 'diffexp_deseq2.tab.gz', gzipped=True)
