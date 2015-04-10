# pylint: disable=missing-docstring
from server.models import Data

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
        self.assertFiles(cuffnorm, 'expset', 'expression_set.tsv.gz', gzipped=True)

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

        mappa = self.run_processor("import:upload:mappability", {"src": "purpureum_mappability_50.tab.gz"})

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'mappable': mappa.pk}
        expression = self.run_processor('expression:bcm-1-0-0', inputs)
        self.assertFiles(expression, 'rpkm', 'expression_bcm_rpkm.tab.gz', gzipped=True)

        inputs = {'expressions': [expression.pk, expression.pk]}
        etc = self.run_processor('etc:bcm-1-0-0', inputs)
        self.assertJSON(etc, etc.output['etc'], '', 'etc.json.gz')

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
            'mappable': mappa.pk}
        expression2 = self.run_processor('expression:bcm-1-0-0', inputs)

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
        expression2 = self.run_processor('htseq-count:-0-6-1p1', inputs)
        self.assertFiles(expression, 'rc', 'reads_rc.tab.gz', gzipped=True)
        self.assertFiles(expression, 'fpkm', 'reads_fpkm.tab.gz', gzipped=True)
        self.assertFiles(expression, 'exp', 'reads_tpm.tab.gz', gzipped=True)
        self.assertJSON(expression, expression.output['exp_json'], '', 'expression_htseq.json.gz')

    def test_mergeexpression(self):
        expression_1 = self.prepare_expression(f_rc='00Hr_rc.tab.gz', f_exp='00Hr_tpm.tab.gz', f_type="TPM")
        expression_2 = self.prepare_expression(f_rc='20Hr_rc.tab.gz', f_exp='20Hr_tpm.tab.gz', f_type="TPM")
        expression_3 = self.prepare_expression(f_rc='20Hr_rc.tab.gz', f_exp='20Hr_tpm.tab.gz', f_type="RC")

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067098', 'DPU_G0067102']
        }

        mergeexpression_1 = self.run_processor('mergeexpressions', inputs)
        self.assertFiles(mergeexpression_1, "expset", "merged_expset_subset.tab")

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': []
        }

        mergeexpression_2 = self.run_processor('mergeexpressions', inputs)
        self.assertFiles(mergeexpression_2, "expset", "merged_expset_all.tab")

        inputs = {
            'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067098', 'DPU_G0067102']
        }
        mergeexpression_3 = self.run_processor('mergeexpressions', inputs, Data.STATUS_ERROR)
