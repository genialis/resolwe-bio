from .base import BaseProcessorTestCase


class ExpressionProcessorTestCase(BaseProcessorTestCase):
    def prepair_genome(self):
        inputs = {'src': 'genome.fasta.gz'}
        genome = self.run_processor('import:upload:genome-fasta', inputs)
        self.assertDone(genome)
        self.assertFiles(genome, 'fasta', 'genome.fasta.gz')
        return genome

    def prepair_reads(self):
        inputs = {'src': 'reads.fastq.gz'}
        reads = self.run_processor('import:upload:reads-fastq', inputs)
        self.assertDone(reads)
        self.assertFields(reads, 'bases', 35)
        return reads

    def test_cufflinks(self):
        genome = self.prepair_genome()
        reads = self.prepair_reads()

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertDone(annotation)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertDone(aligned_reads)

        inputs = {'alignment': aligned_reads.pk, 'gff': annotation.pk, 'genome': genome.pk}
        cuff_exp = self.run_processor('cufflinks:-2-2-1', inputs)
        self.assertDone(cuff_exp)
        self.assertFiles(cuff_exp, 'transcripts', 'cufflinks_transcripts.gtf')

        inputs = {'alignment': aligned_reads.pk, 'gff': annotation.pk, 'genome': genome.pk}
        cuff_exp2 = self.run_processor('cufflinks:-2-2-1', inputs)
        self.assertDone(cuff_exp2)

        inputs = {'expressions': [cuff_exp.pk, cuff_exp2.pk], 'gff': annotation.pk, 'genome': genome.pk}
        cuff_merge = self.run_processor('cuffmerge:-2-2-1', inputs)
        self.assertDone(cuff_merge)
        self.assertFiles(cuff_merge, 'merged_gtf', 'cuffmerge_transcripts.gtf')

        inputs = {
            'alignments': [aligned_reads.pk, aligned_reads.pk],
            'replicates': ['1', '2'],
            'labels': ['g1', 'g2'],
            'gff': cuff_merge.pk}
        cuffnorm = self.run_processor('cuffnorm:-2-2-1', inputs)
        self.assertDone(cuffnorm)
        self.assertFiles(cuffnorm, 'isoforms_fpkm_tracking', 'cuffnorm_output')

    def test_expression_bcm(self):
        genome = self.prepair_genome()
        reads = self.prepair_reads()

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertDone(annotation)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertDone(aligned_reads)

        inputs = {'genome': genome.pk, 'gff': annotation.pk}
        mappability = self.run_processor('mappability:bcm-1-0-0', inputs)
        self.assertDone(mappability)
        self.assertFiles(mappability, 'mappability', 'mappability.tab')

        inputs = {'alignment': aligned_reads.pk, 'gff': annotation.pk, 'mappable': mappability.pk}
        expression = self.run_processor('expression:bcm-1-0-0', inputs)
        self.assertDone(expression)
        self.assertFiles(expression, 'rpkm', 'expression_bcm_rpkm.tab')

        inputs = {'expressions': [expression.pk, expression.pk]}
        etc = self.run_processor('etc:bcm-1-0-0', inputs)
        self.assertDone(etc)
