# pylint: disable=missing-docstring
from resolwe.flow.models import Data

from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase
from resolwe.flow.models import Data


class ExpressionProcessorTestCase(BioProcessTestCase):

    def test_cufflinks(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff.gz'}
        annotation = self.run_processor('upload-gff3', inputs)

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment-tophat2', inputs)

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_exp = self.run_processor('cufflinks', inputs)
        self.assertFile(cuff_exp, 'transcripts', 'cufflinks_transcripts.gtf')

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_exp2 = self.run_processor('cufflinks', inputs)

        inputs = {
            'expressions': [cuff_exp.pk, cuff_exp2.pk],
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_merge = self.run_processor('cuffmerge', inputs)
        self.assertFile(cuff_merge, 'merged_gtf', 'cuffmerge_transcripts.gtf')

    def test_cuffquant(self):
        inputs = {"src": "cuffquant_mapping.bam"}
        bam = self.run_processor("upload-bam", inputs)

        annotation = self.prepare_annotation(fn='hg19_chr20_small.gtf.gz')

        inputs = {
            'alignment': bam.id,
            'gff': annotation.id}
        self.run_processor('cuffquant', inputs)

    def test_cuffnorm(self):
        inputs = {"src": "cuffquant_1.cxb"}
        sample_1 = self.run_processor("upload-cxb", inputs)

        inputs = {"src": "cuffquant_2.cxb"}
        sample_2 = self.run_processor("upload-cxb", inputs)

        annotation = self.prepare_annotation(fn='hg19_chr20_small.gtf.gz')

        inputs = {
            'cuffquant': [sample_1.pk, sample_2.pk],
            'annotation': annotation.id,
            'labels': ['Group1', 'Group2'],
            'replicates': ['1', '2']}
        cuffnorm = self.run_processor('cuffnorm', inputs)
        self.assertFile(cuffnorm, 'genes_fpkm', 'cuffnorm_genes.fpkm_table')
        self.assertFile(cuffnorm, 'raw_scatter', 'cuffnorm_scatter_plot.png')

        exp = Data.objects.last()
        self.assertFile(exp, 'exp', 'cuffnorm_expression.tab.gz', compression='gzip')

    def test_expression_bcm(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff.gz'}
        annotation = self.run_processor('upload-gff3', inputs)

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment-tophat2', inputs)

        mappa = self.run_processor("upload-mappability", {"src": "purpureum_mappability_50.tab.gz"})

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'mappable': mappa.pk}
        expression = self.run_processor('expression:bcm', inputs)
        self.assertFile(expression, 'rpkm', 'expression_bcm_rpkm.tab.gz', compression='gzip')

        inputs = {'expressions': [expression.pk, expression.pk]}
        etc = self.run_processor('etc-bcm', inputs)
        self.assertJSON(etc, etc.output['etc'], '', 'etc.json.gz')

    def test_expression_htseq(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gtf.gz'}
        annotation = self.run_processor('upload-gtf', inputs)

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment-tophat2', inputs)

        inputs = {
            'alignments': aligned_reads.pk,
            'gff': annotation.pk,
            'stranded': "no",
            'id_attribute': 'transcript_id'}
        expression = self.run_processor('htseq-count', inputs)
        self.assertFile(expression, 'rc', 'reads_rc.tab.gz', compression='gzip')
        self.assertFile(expression, 'fpkm', 'reads_fpkm.tab.gz', compression='gzip')
        self.assertFile(expression, 'exp', 'reads_tpm.tab.gz', compression='gzip')
        self.assertJSON(expression, expression.output['exp_json'], '', 'expression_htseq.json.gz')

    def test_mergeexpression(self):
        expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM")
        expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")
        expression_3 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="RC")

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067098', 'DPU_G0067102']
        }

        mergeexpression_1 = self.run_processor('mergeexpressions', inputs)
        self.assertFile(mergeexpression_1, "expset", "merged_expset_subset.tab")

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'genes': []
        }

        mergeexpression_2 = self.run_processor('mergeexpressions', inputs)
        self.assertFile(mergeexpression_2, "expset", "merged_expset_all.tab")

        inputs = {
            'exps': [expression_1.pk, expression_2.pk, expression_3.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067098', 'DPU_G0067102']
        }
        self.run_processor('mergeexpressions', inputs, Data.STATUS_ERROR)

    def test_etcmerge(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff.gz'}
        annotation = self.run_processor('upload-gff3', inputs)

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment-tophat2', inputs)

        mappa = self.run_processor("upload-mappability", {"src": "purpureum_mappability_50.tab.gz"})

        inputs = {
            'alignment': aligned_reads.pk,
            'gff': annotation.pk,
            'mappable': mappa.pk}

        expression = self.run_processor('expression:bcm', inputs)

        inputs = {'expressions': [expression.pk, expression.pk]}
        etc = self.run_processor('etc-bcm', inputs)

        inputs = {
            'exps': [etc.pk],
            'genes': ['DPU_G0067110', 'DPU_G0067098', 'DPU_G0067102']
        }

        etcmerge = self.run_processor('mergeetc', inputs)
        self.assertFile(etcmerge, "expset", "merged_etc.tab.gz", compression='gzip')

    def test_ncrna(self):
        inputs = {"src": "ncRNA_sample1.bam"}
        sample_1 = self.run_processor("upload-bam", inputs)

        inputs = {"src": "ncRNA_sample2.bam"}
        sample_2 = self.run_processor("upload-bam", inputs)

        inputs = {'src': 'ncRNA_annotation.gff.gz'}
        annotation = self.run_processor('upload-gff3', inputs)

        inputs = {"src": "ncRNA_genome.fasta.gz"}
        genome = self.run_processor('upload-genome', inputs)

        inputs = {
            'alignment': sample_1.pk,
            'gff': annotation.pk,
            'library_type': "fr-firststrand"}
        cuff_exp_1 = self.run_processor('cufflinks', inputs)

        inputs = {
            'alignment': sample_2.pk,
            'gff': annotation.pk,
            'library_type': "fr-firststrand"}
        cuff_exp_2 = self.run_processor('cufflinks', inputs)

        inputs = {
            'expressions': [cuff_exp_1.pk, cuff_exp_2.pk],
            'gff': annotation.pk,
            'genome': genome.pk}
        merged_annotation = self.run_processor('cuffmerge', inputs)

        annotation_gff3 = self.run_processor('cuffmerge-gtf-to-gff3', {"cuffmerge": merged_annotation.pk})

        inputs = {"genome": genome.pk, "gff": annotation_gff3.pk, "length": 100}
        mappa = self.run_processor("mappability-bcm", inputs)

        inputs = {
            'alignment': sample_1.pk,
            'gff': annotation_gff3.pk,
            'mappable': mappa.pk,
            'stranded': True}
        expression_1 = self.run_processor('expression-bcm-ncrna', inputs)

        inputs = {
            'alignment': sample_2.pk,
            'gff': annotation_gff3.pk,
            'mappable': mappa.pk,
            'stranded': True}
        expression_2 = self.run_processor('expression-bcm-ncrna', inputs)

        inputs = {
            'exps': [expression_1.pk, expression_2.pk],
            'annotation': annotation_gff3.pk}
        ncrna_expressions = self.run_processor('summarizexpressions-ncrna', inputs)
        self.assertFile(ncrna_expressions, 'expset', 'ncRNA_exp_all.tab.gz', compression='gzip')
        self.assertFile(ncrna_expressions, 'ncrna', 'ncRNA_exp.tab.gz', compression='gzip')
