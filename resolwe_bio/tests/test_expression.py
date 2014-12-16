from .base import BaseProcessorTestCase


class ExpressionProcessorTestCase(BaseProcessorTestCase):
    def test_cufflinks(self):
        inputs = {'src': 'tophat_reads_mapped.bam'}
        mapping = self.run_processor('import:upload:mapping-bam', inputs)
        self.assertDone(mapping)
        self.assertFiles(mapping, 'bam', 'tophat_reads_mapped.bam')

        inputs = {'src': 'tophat_reads_mapped.bam'}
        mapping2 = self.run_processor('import:upload:mapping-bam', inputs)
        self.assertDone(mapping2)
        self.assertFiles(mapping2, 'bam', 'tophat_reads_mapped.bam')

        inputs = {'src': 'genome.fasta.gz'}
        genome = self.run_processor('import:upload:genome-fasta', inputs)
        self.assertDone(genome)
        self.assertFiles(genome, 'fasta', 'genome.fasta.gz')

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertDone(annotation)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {'alignment': mapping.pk, 'gff': annotation.pk, 'genome': genome.pk}
        cuff_exp = self.run_processor('cufflinks:-2-2-1', inputs)
        self.assertDone(cuff_exp)
        self.assertFiles(cuff_exp, 'transcripts', 'cufflinks_transcripts.gtf')

        cuff_exp2 = cuff_exp = self.run_processor('cufflinks:-2-2-1', inputs)
        self.assertDone(cuff_exp2)
        self.assertFiles(cuff_exp, 'transcripts', 'cufflinks_transcripts.gtf')

        inputs = {'expressions': str([cuff_exp.pk, cuff_exp2.pk]), 'gff': annotation.pk, 'genome': genome.pk}
        cuff_merge = self.run_processor('cuffmerge:-2-2-1', inputs)
        self.assertDone(cuff_merge)
        self.assertFiles(cuff_merge, 'merged_gtf', 'cuffmerge_annotation.gtf')

        inputs = {'alignments': str([mapping.pk, mapping2.pk]), 'replicates': '1,2', 'labels': 'g1,g2', 'gff': cuff_merge.pk}
        cuffnorm = self.run_processor('cuffnorm:-2-2-1', inputs)
        self.assertDone(cuffnorm)
        self.assertFiles(cuffnorm, 'isoforms_fpkm_tracking', 'cuffnorm_output')

    def test_expression_bcm(self):
        inputs = {'src': 'tophat_reads_mapped.bam'}
        mapping = self.run_processor('import:upload:mapping-bam', inputs)
        self.assertDone(mapping)
        self.assertFiles(mapping, 'bam', 'tophat_reads_mapped.bam')

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertDone(annotation)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {'src': 'genome.fasta.gz'}
        genome = self.run_processor('import:upload:genome-fasta', inputs)
        self.assertDone(genome)
        self.assertFiles(genome, 'fasta', 'genome.fasta.gz')

        inputs = {'genome': genome.pk, 'gff': annotation.pk}
        mappability = self.run_processor('mappability:bcm-1-0-0', inputs)
        self.assertDone(mappability)
        self.assertFiles(mappability, 'mappability', 'mappability.tab')

        inputs = {'alignment': mapping.pk, 'gff': annotation.pk, 'mappable': mappability.pk}
        expression = self.run_processor('expression:bcm-1-0-0', inputs)
        self.assertDone(expression)
        self.assertFiles(expression, 'rpkm', 'expression_bcm_rpkm.tab')
