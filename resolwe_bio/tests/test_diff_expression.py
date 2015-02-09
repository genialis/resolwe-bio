from .base import BaseProcessorTestCase


class DiffExpProcessorTestCase(BaseProcessorTestCase):
    def prepair_genome(self):
        inputs = {'src': 'genome.fasta.gz'}
        genome = self.run_processor('import:upload:genome-fasta', inputs)
        self.assertDone(genome)
        self.assertFiles(genome, 'fasta', 'genome.fasta.gz')
        return genome

    def test_cuffdiff(self):
        genome = self.prepair_genome()

        inputs = {'src': '00Hr.fastq.gz'}
        reads1 = self.run_processor('import:upload:reads-fastq', inputs)
        self.assertDone(reads1)

        inputs = {'src': '20Hr.fastq.gz'}
        reads2 = self.run_processor('import:upload:reads-fastq', inputs)
        self.assertDone(reads2)

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertDone(annotation)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {
            'genome': genome.pk,
            'reads': reads1.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_1 = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertDone(aligned_reads_1)

        inputs = {
            'genome': genome.pk,
            'reads': reads2.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_2 = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertDone(aligned_reads_2)

        inputs = {
            'alignment': aligned_reads_1.pk,
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_exp = self.run_processor('cufflinks:-2-2-1', inputs)
        self.assertDone(cuff_exp)

        inputs = {
            'alignment': aligned_reads_2.pk,
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_exp2 = self.run_processor('cufflinks:-2-2-1', inputs)
        self.assertDone(cuff_exp2)

        inputs = {
            'expressions': [cuff_exp.pk, cuff_exp2.pk],
            'gff': annotation.pk,
            'genome': genome.pk}
        cuff_merge = self.run_processor('cuffmerge:-2-2-1', inputs)
        self.assertDone(cuff_merge)

        inputs = {
            'alignment': aligned_reads_1.pk,
            'gff': cuff_merge.pk}
        cuffquant = self.run_processor('cuffquant:-2-2-1', inputs)
        self.assertDone(cuffquant)

        inputs = {
            'alignment': aligned_reads_2.pk,
            'gff': cuff_merge.pk}
        cuffquant2 = self.run_processor('cuffquant:-2-2-1', inputs)
        self.assertDone(cuffquant2)

        inputs = {
            'cuffquant': [cuffquant.pk, cuffquant2.pk],
            'replicates': ['1', '2'],
            'labels': ['g1', 'g2'],
            'gff': cuff_merge.pk}
        cuffdiff = self.run_processor('cuffdiff:-2-2-1', inputs)
        self.assertDone(cuffdiff)
        self.assertFiles(cuffdiff, 'gene_diff_exp', 'cuffdiff_output')

    def test_bayseq_bcm(self):
        genome = self.prepair_genome()

        inputs = {'src': '00Hr.fastq.gz'}
        reads1 = self.run_processor('import:upload:reads-fastq', inputs)
        self.assertDone(reads1)

        inputs = {'src': '20Hr.fastq.gz'}
        reads2 = self.run_processor('import:upload:reads-fastq', inputs)
        self.assertDone(reads2)

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertDone(annotation)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {
            'genome': genome.pk,
            'reads': reads1.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_1 = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertDone(aligned_reads_1)

        inputs = {
            'genome': genome.pk,
            'reads': reads2.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_2 = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertDone(aligned_reads_2)

        inputs = {
            'genome': genome.pk,
            'gff': annotation.pk}
        mappability = self.run_processor('mappability:bcm-1-0-0', inputs)
        self.assertDone(mappability)

        inputs = {
            'alignment': aligned_reads_1.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression_1 = self.run_processor('expression:bcm-1-0-0', inputs)
        self.assertDone(expression_1)

        inputs = {
            'alignment': aligned_reads_2.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression_2 = self.run_processor('expression:bcm-1-0-0', inputs)
        self.assertDone(expression_2)

        inputs = {
            'name': "00vs20",
            'case': [expression_1.pk],
            'control': [expression_2.pk],
            'replicates': ['1', '2'],
            'mappability': mappability.pk}
        diff_exp = self.run_processor('differentialexpression:bcm-1-0-0', inputs)
        self.assertDone(diff_exp)
        # self.assertJSON(diff_exp.volcano_plot, '', 'bayseq_volcano.json')
