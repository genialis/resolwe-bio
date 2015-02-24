from .base import BaseProcessorTestCase


class GencoverProcessorTestCase(BaseProcessorTestCase):
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

    def test_coverage(self):
        genome = self.prepair_genome()
        reads = self.prepair_reads()

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertDone(annotation)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        # GTF inport
        inputs = {'src': 'annotation_ok.gtf'}
        annotation_gtf = self.run_processor('import:upload:annotation-gtf', inputs)
        self.assertDone(annotation_gtf)
        self.assertFiles(annotation_gtf, 'gtf', 'annotation_ok.gtf')

        # redundant GTF inport
        inputs = {'src': 'annotation_red.gtf'}
        annotation_gtf_red = self.run_processor('import:upload:annotation-gtf', inputs)
        self.assertDone(annotation_gtf_red)
        self.assertFiles(annotation_gtf_red, 'gtf', 'annotation_red.gtf')

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertDone(aligned_reads)

        # samtools mapping
        inputs = {
            'mapping': aligned_reads.pk,
            'genome': genome.pk}
        variants = self.run_processor('vc-samtools', inputs)
        self.assertDone(variants)

        # Coverage report
        inputs = {
            'mapping': aligned_reads.pk,
            'gtf': annotation_gtf.pk,
            'variants': variants.pk,
            'filter': 3,
            'genes': ['geneX']}

        genc_results = self.run_processor('coverage:garvan', inputs)
        self.assertDone(genc_results)
        # self.assertFiles(genc_results, 'bigwig', 'genome_coverage.bw')

        # Missing gene in BAM file test
        inputs = {
            'mapping': aligned_reads.pk,
            'gtf': annotation_gtf_red.pk,
            'variants': variants.pk,
            'filter': 3,
            'genes': ['geneX']}

        genc_results = self.run_processor('coverage:garvan', inputs)
        self.assertFields(genc_results, 'proc.warning', 'Contig scaffold_fake not found in BAM file (for gene geneX)')
        self.assertDone(genc_results)
