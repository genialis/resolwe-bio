# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class ScSeqProcessorTestCase(BioProcessTestCase):

    @tag_process('cellranger-mkref')
    def test_cellranger_mkref(self):
        with self.preparation_stage():
            annotation = self.prepare_annotation(fn='HS chr21_short_ensembl.gtf.gz', source='ENSEMBL',
                                                 species='Homo sapiens', build='GRCh38.93')
            inputs = {
                'src': 'HS chr21_short_ensembl.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'GRCh38.93',
            }
            genome = self.run_process('upload-genome', inputs)

        inputs = {'annotation': annotation.id, 'genome': genome.id}
        mkref = self.run_process('cellranger-mkref', inputs)

        self.assertAlmostEqual(mkref.output['genome_index']['size'], 1429926, delta=1)
        self.assertFields(mkref, 'build', 'GRCh38.93')
        self.assertFields(mkref, 'species', 'Homo sapiens')
        self.assertFields(mkref, 'source', 'ENSEMBL')

    @tag_process('cellranger-count')
    def test_cellranger_count(self):
        with self.preparation_stage():
            annotation = self.prepare_annotation(fn='HS chr21_short_ensembl.gtf.gz', source='ENSEMBL',
                                                 species='Homo sapiens', build='GRCh38.93')
            inputs = {
                'src': 'HS chr21_short_ensembl.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'GRCh38.93',
            }
            genome = self.run_process('upload-genome', inputs)

            inputs = {'annotation': annotation.id, 'genome': genome.id}
            genome_index = self.run_process('cellranger-mkref', inputs)

            inputs = {
                'barcodes': ['10x_S1_L001_R1_001.fastq.gz', '10x_S1_L002_R1_001.fastq.gz'],
                'reads': ['10x_S1_L001_R2_001.fastq.gz', '10x_S1_L002_R2_001.fastq.gz'],
            }
            reads = self.run_process('upload-sc-10x', inputs)

        inputs = {
            'reads': reads.id,
            'genome_index': genome_index.id,
        }
        count = self.run_process('cellranger-count', inputs)
        self.assertFile(count, 'matrix_filtered', '10x_scseq_matrix.mtx.gz', compression='gzip')
        self.assertFields(count, 'build', 'GRCh38.93')
        self.assertFields(count, 'species', 'Homo sapiens')
        self.assertFields(count, 'source', 'ENSEMBL')

        # Test 'upload-bam-scseq-indexed' process
        bam = Data.objects.last()
        self.assertFileExists(bam, 'bam')
        self.assertFileExists(bam, 'bai')
        self.assertFile(bam, 'stats', '10x_scseq_stats.txt')
        self.assertFields(bam, 'species', 'Homo sapiens')
        self.assertFields(bam, 'build', 'GRCh38.93')
