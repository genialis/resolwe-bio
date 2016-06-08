# pylint: disable=missing-docstring
from .utils import skipDockerFailure, BioProcessTestCase

class AlignmentProcessorTestCase(BioProcessTestCase):

    def test_bowtie(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'reporting': {'r': "-a -m 1 --best --strata"}}
        aligned_reads = self.run_processor('alignment-bowtie', inputs)

    def test_bowtie2(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment-bowtie2', inputs)
        self.assertFile(aligned_reads, 'stats', 'bowtie2_reads_report.txt')

    def test_tophat(self):
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
        self.assertFile(aligned_reads, 'stats', 'tophat_reads_report.txt')

    @skipDockerFailure("Fails with: STAR: command not found")
    def test_star(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff.gz'}
        annotation = self.run_processor('upload-gff3', inputs)

        inputs = {'genome': genome.pk, 'annotation': annotation.pk}
        genome_index = self.run_processor('alignment-star-index', inputs)

        inputs = {
            'genome': genome_index.pk,
            'reads': reads.pk,
            'threads': 1,
            't_coordinates': {
                'quantmode': True,
                'gene_counts': True}}
        aligned_reads = self.run_processor('alignment-star', inputs)
        self.assertFile(aligned_reads, 'gene_counts', 'gene_counts_star.tab.gz', compression='gzip')

    @skipDockerFailure("File contents mismatch for bwa_bt_reads_mapped.bam")
    def test_bwa_bt(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment-bwa-aln', inputs)
        self.assertFile(aligned_reads, 'bam', 'bwa_bt_reads_mapped.bam')
        self.assertFile(aligned_reads, 'stats', 'bwa_bt_reads_report.txt')

    def test_bwa_sw(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment-bwa-sw', inputs)
        self.assertFile(aligned_reads, 'bam', 'bwa_sw_reads_mapped.bam')
        self.assertFile(aligned_reads, 'stats', 'bwa_sw_reads_report.txt')

    @skipDockerFailure("File contents mismatch for bwa_mem_reads_mapped.bam")
    def test_bwa_mem(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment-bwa-mem', inputs)
        self.assertFile(aligned_reads, 'bam', 'bwa_mem_reads_mapped.bam')
        self.assertFile(aligned_reads, 'stats', 'bwa_mem_reads_report.txt')

    def test_hisat2(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk}
        aligned_reads = self.run_processor('alignment-hisat2', inputs)
        self.assertFile(aligned_reads, 'stats', 'hisat2_report.txt')
