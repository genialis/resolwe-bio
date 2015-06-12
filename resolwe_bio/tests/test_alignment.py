# pylint: disable=missing-docstring
from .base import BaseProcessorTestCase
from .utils import PreparedData


class AlignmentProcessorTestCase(BaseProcessorTestCase, PreparedData):

    def test_bowtie(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()
        filtered_single = self.prepare_filtered_reads_single()
        filtered_paired = self.prepare_filtered_reads_paired()

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'reporting': {'r': "-a -m 1 --best --strata"}}
        aligned_reads = self.run_processor('alignment:bowtie-1-0-0-trimmx', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bowtie_reads_report.tab.gz', compression='gzip')

        inputs = {
            'genome': genome.pk,
            'reads': filtered_single.pk,
            'reporting': {'r': "-a -m 1 --best --strata"}}
        aligned_reads = self.run_processor('alignment:bowtie-1-0-0-trimmx', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bowtie_reads_report_filtered_single.tab.gz', compression='gzip')

        inputs = {
            'genome': genome.pk,
            'reads': filtered_paired.pk,
            'reporting': {'r': "-a -m 1 --best --strata"}}
        aligned_reads = self.run_processor('alignment:bowtie-1-0-0-trimmx', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bowtie_reads_report_filtered_paired.tab.gz', compression='gzip')

    def test_bowtie2(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()
        filtered_single = self.prepare_filtered_reads_single()
        filtered_paired = self.prepare_filtered_reads_paired()

        inputs = {
            'genome': genome.pk,
            'reads': reads.pk,
            'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment:bowtie-2-2-3_trim', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bowtie2_reads_report.txt')

        inputs = {
            'genome': genome.pk,
            'reads': filtered_single.pk,
            'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment:bowtie-2-2-3_trim', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bowtie2_report_filtered_single.txt')

        inputs = {
            'genome': genome.pk,
            'reads': filtered_paired.pk,
            'reporting': {'rep_mode': "def"}}
        aligned_reads = self.run_processor('alignment:bowtie-2-2-3_trim', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bowtie2_report_filtered_paired.txt')

    def test_tophat(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()
        filtered_single = self.prepare_filtered_reads_single()
        filtered_paired = self.prepare_filtered_reads_paired()

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
        self.assertFiles(aligned_reads, 'stats', 'tophat_reads_report.txt')

        inputs = {
            'genome': genome.pk,
            'reads': filtered_single.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertFiles(aligned_reads, 'stats', 'tophat_reads_report_filtered_single.txt')

        inputs = {
            'genome': genome.pk,
            'reads': filtered_paired.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertFiles(aligned_reads, 'stats', 'tophat_reads_report_filtered_paired.txt')

    def test_bwa_bt(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()
        filtered_single = self.prepare_filtered_reads_single()
        filtered_paired = self.prepare_filtered_reads_paired()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment:bwa_aln-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'bam', 'bwa_bt_reads_mapped.bam')
        self.assertFiles(aligned_reads, 'stats', 'bwa_bt_reads_report.txt')

        inputs = {'genome': genome.pk, 'reads': filtered_single.pk}
        aligned_reads = self.run_processor('alignment:bwa_aln-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bwa_bt_reads_report_filtered_single.txt')

        inputs = {'genome': genome.pk, 'reads': filtered_paired.pk}
        aligned_reads = self.run_processor('alignment:bwa_aln-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bwa_bt_reads_report_filtered_paired.txt')

    def test_bwa_sw(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()
        filtered_single = self.prepare_filtered_reads_single()
        filtered_paired = self.prepare_filtered_reads_paired()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment:bwa_sw-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'bam', 'bwa_sw_reads_mapped.bam')
        self.assertFiles(aligned_reads, 'stats', 'bwa_sw_reads_report.txt')

        inputs = {'genome': genome.pk, 'reads': filtered_single.pk}
        aligned_reads = self.run_processor('alignment:bwa_sw-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bwa_sw_reads_report_filtered_single.txt')

        inputs = {'genome': genome.pk, 'reads': filtered_paired.pk}
        aligned_reads = self.run_processor('alignment:bwa_sw-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bwa_sw_reads_report_filtered_paired.txt')

    def test_bwa_mem(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()
        filtered_single = self.prepare_filtered_reads_single()
        filtered_paired = self.prepare_filtered_reads_paired()

        inputs = {'genome': genome.pk, 'reads': reads.pk}
        aligned_reads = self.run_processor('alignment:bwa_mem-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'bam', 'bwa_mem_reads_mapped.bam')
        self.assertFiles(aligned_reads, 'stats', 'bwa_mem_reads_report.txt')

        inputs = {'genome': genome.pk, 'reads': filtered_single.pk}
        aligned_reads = self.run_processor('alignment:bwa_mem-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bwa_mem_reads_report_filtered_single.txt')

        inputs = {'genome': genome.pk, 'reads': filtered_paired.pk}
        aligned_reads = self.run_processor('alignment:bwa_mem-0.7.5a', inputs)
        self.assertFiles(aligned_reads, 'stats', 'bwa_mem_reads_report_filtered_paired.txt')
