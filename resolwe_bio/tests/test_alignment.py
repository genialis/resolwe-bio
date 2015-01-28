from .base import BaseProcessorTestCase


class AlignmentProcessorTestCase(BaseProcessorTestCase):
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

    def test_bowtie(self):
        genome = self.prepair_genome()
        reads = self.prepair_reads()

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'reporting': {'r': "-a -m 1 --best --strata"}}
        aligned_reads = self.run_processor('alignment:bowtie-1-0-0-trimmx', inputs)
        self.assertDone(aligned_reads)
        self.assertFiles(aligned_reads, 'stats', 'bowtie_reads_report.tab')

    def test_bowtie2(self):
        genome = self.prepair_genome()
        reads = self.prepair_reads()

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment:bowtie-2-2-3_trim', inputs)
        self.assertDone(aligned_reads)
        self.assertFiles(aligned_reads, 'stats', 'bowtie2_reads_report.txt')

    def test_tophat(self):
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
        self.assertFiles(aligned_reads, 'stats', 'tophat_reads_report.txt')

    def test_bwa_bt(self):
        genome = self.prepair_genome()
        reads = self.prepair_reads()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment:bwa_aln-0.7.5a', inputs)
        self.assertDone(aligned_reads)
        self.assertFiles(aligned_reads, 'bam', 'bwa_bt_reads_mapped.bam')
        self.assertFiles(aligned_reads, 'stats', 'bwa_bt_reads_report.txt')

    def test_bwa_sw(self):
        genome = self.prepair_genome()
        reads = self.prepair_reads()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment:bwa_sw-0.7.5a', inputs)
        self.assertDone(aligned_reads)
        self.assertFiles(aligned_reads, 'bam', 'bwa_sw_reads_mapped.bam')
        self.assertFiles(aligned_reads, 'stats', 'bwa_sw_reads_report.txt')

    def test_bwa_mem(self):
        genome = self.prepair_genome()
        reads = self.prepair_reads()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment:bwa_mem-0.7.5a', inputs)
        self.assertDone(aligned_reads)
        self.assertFiles(aligned_reads, 'bam', 'bwa_mem_reads_mapped.bam')
        self.assertFiles(aligned_reads, 'stats', 'bwa_mem_reads_report.txt')
