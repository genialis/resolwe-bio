from .base import BaseProcessorTestCase
from .utils import PreparedData


class ExpressionProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_cufflinks(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_exp = self.run_processor('cufflinks:-2-2-1', inputs)
        self.assertFiles(cuff_exp, 'transcripts', 'cufflinks_transcripts.gtf')

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_exp2 = self.run_processor('cufflinks:-2-2-1', inputs)

        inputs = {
            'expressions': [cuff_exp.pk, cuff_exp2.pk],
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_merge = self.run_processor('cuffmerge:-2-2-1', inputs)
        self.assertFiles(cuff_merge, 'merged_gtf', 'cuffmerge_transcripts.gtf')

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': cuff_merge.pk}
        cuffquant = self.run_processor('cuffquant:-2-2-1', inputs)

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': cuff_merge.pk}
        cuffquant2 = self.run_processor('cuffquant:-2-2-1', inputs)

        inputs = {
            'cuffquant': [cuffquant.pk, cuffquant2.pk],
            'replicates': ['1', '2'],
            'labels': ['g1', 'g2'],
            'gff': cuff_merge.pk}
        cuffnorm = self.run_processor('cuffnorm:-2-2-1', inputs)
        self.assertFiles(cuffnorm, 'expression_set', 'expression_set.tsv.gz', gzipped=True)

    def test_expression_bcm(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)

        inputs = {
            'genome': genome.pk,
            'gff': annotation.pk}
        mappability = self.run_processor('mappability:bcm-1-0-0', inputs)
        self.assertFiles(mappability, 'mappability', 'mappability.tab')

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression = self.run_processor('expression:bcm-1-0-0', inputs)
        self.assertFiles(expression, 'rpkm', 'expression_bcm_rpkm.tab.gz', gzipped=True)

        inputs = {'expressions': [expression.pk, expression.pk]}
        etc = self.run_processor('etc:bcm-1-0-0', inputs)
        self.assertJSON(etc, etc.output['etc'], '', 'etc.json')

        reads2 = self.prepare_reads('00Hr.fastq.gz')
        inputs = {
            'genome': genome.pk,
            'reads': reads2.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads2 = self.run_processor('alignment:tophat-2-0-13', inputs)

        inputs = {
            'alignment': aligned_reads2.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression2 = self.run_processor('expression:bcm-1-0-0', inputs)

        inputs = {
            'expressions': [expression.pk, expression2.pk],
            'genes': ['DDB_G0267184', 'DDB_G0267188', 'DDB_G0267204']
        }

        self.run_processor('mergeexpressions', inputs)

    def test_expression_htseq(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gtf'}
        annotation = self.run_processor('import:upload:annotation-gtf', inputs)
        self.assertFiles(annotation, 'gtf', 'annotation.gtf')

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)

        inputs = {
            'alignments': aligned_reads.pk,
            'gff': annotation.pk,
            'stranded': "no",
            'id_attribute': 'transcript_id'}
        expression = self.run_processor('htseq-count:-0-6-1p1', inputs)
        self.assertFiles(expression, 'rc', 'reads_rc.tab.gz', gzipped=True)
        self.assertFiles(expression, 'fpkm', 'reads_fpkm.tab.gz', gzipped=True)
        self.assertFiles(expression, 'tpm', 'reads_tpm.tab.gz', gzipped=True)
        self.assertJSON(expression, expression.output['exp'], '', 'expression.json')
