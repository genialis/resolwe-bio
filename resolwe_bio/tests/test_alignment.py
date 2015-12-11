# pylint: disable=missing-docstring
from .utils import BioProcessTestCase


class AlignmentProcessorTestCase(BioProcessTestCase):

    def test_bowtie(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'reporting': {'r': "-a -m 1 --best --strata"}}
        aligned_reads = self.run_processor('alignment:bowtie-1-0-0-trimmx', inputs)

    def test_bowtie2(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment:bowtie-2-2-3_trim', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bowtie2_reads_report.txt')

    def test_tophat(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff.gz'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertFiles(aligned_reads, 'stats', 'tophat_reads_report.txt')

    def test_star(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff.gz'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)

        inputs = {'genome': genome.pk, 'annotation': annotation.pk}
        genome_index = self.run_processor('alignment:star:index', inputs)

        inputs = {
            'genome': genome_index.pk,
            'reads': reads.pk,
            'threads': 1,
            't_coordinates': {
                'quantmode': True,
                'gene_counts': True}}
        aligned_reads = self.run_processor('alignment:star', inputs)
        self.assertFiles(aligned_reads, 'gene_counts', 'gene_counts_star.tab.gz', compression='gzip')

    def test_bwa_bt(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment:bwa_aln-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'bam', 'bwa_bt_reads_mapped.bam')
        self.assertFiles(aligned_reads, 'stats', 'bwa_bt_reads_report.txt')

    def test_bwa_sw(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment:bwa_sw-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'bam', 'bwa_sw_reads_mapped.bam')
        self.assertFiles(aligned_reads, 'stats', 'bwa_sw_reads_report.txt')

    def test_bwa_mem(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment:bwa_mem-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'bam', 'bwa_mem_reads_mapped.bam')
        self.assertFiles(aligned_reads, 'stats', 'bwa_mem_reads_report.txt')
