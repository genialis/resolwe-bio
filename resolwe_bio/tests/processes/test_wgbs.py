# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class WgbsProcessorTestCase(BioProcessTestCase):

    def test_bsmap(self):
        genome = self.prepare_genome()
        reads_paired = self.prepare_paired_reads(mate1=['fw_reads.fastq.gz'],
                                                 mate2=['rw_reads.fastq.gz'])

        inputs = {
            'genome': genome.id,
            'reads': reads_paired.id,
        }
        bsmap = self.run_process('bsmap', inputs)
        self.assertFile(bsmap, 'stats', 'bsmap_reads_report.txt')

    def test_mcall(self):
        inputs = {'src': 'chr1_part.fasta.gz'}
        genome = self.run_process('upload-genome', inputs)

        inputs = {'src': 'wgbs.bam'}
        bam1 = self.run_process('upload-bam', inputs)

        inputs = {'genome': genome.pk,
                  'genome_identifier': 'hg19',
                  'bam': [bam1.pk],
                  'threads': 3}

        mcall = self.run_process('mcall', inputs)
        self.assertFile(mcall, 'stats', 'wgbs.bam_stat.txt')
        self.assertFile(mcall, 'gbed', 'wgbs.bam.G.bed')
        self.assertFile(mcall, 'hgbed', 'wgbs.bam.HG.bed')
