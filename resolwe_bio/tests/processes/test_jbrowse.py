# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class JbrowseProcessorTestCase(BioProcessTestCase):

    def test_refseq_track(self):
        genome = self.prepare_genome()
        refseq_track = self.run_process('jbrowse-refseq', {'refseq': genome.pk})
        self.assertFields(refseq_track, 'refseq_track', {'refs': ['seq'],
                                                         'file': 'seq/refSeqs.json'})

    def test_gff3_track(self):
        inputs = {'src': 'annotation.gff.gz', 'source': 'dictyBase'}
        annotation = self.run_process('upload-gff3', inputs)
        gff = self.run_process('jbrowse-gff3', {'gff': annotation.pk})
        self.assertFields(gff, 'annotation_track', {'refs': ['tracks/annotation'],
                                                    'file': 'trackList.json'})

    def test_gtf_track(self):
        annotation = self.prepare_annotation()
        gtf = self.run_process('jbrowse-gtf', {'gtf': annotation.pk})
        self.assertFields(gtf, 'annotation_track', {'refs': ['tracks/annotation'],
                                                    'file': 'trackList.json'})

    def test_bed_track(self):
        inputs = {'src': 'bed_track.bed'}
        bed_file = self.run_process('upload-bed', inputs)
        bed = self.run_process('jbrowse-bed', {'bed': bed_file.pk})
        self.assertFields(bed, 'bed_track', {'refs': ['tracks/bed'],
                                             'file': 'trackList.json'})

    def test_coverage_track(self):
        inputs = {"src": "alignment_coverage.bam"}
        bam = self.run_process("upload-bam", inputs)

        coverage = self.run_process('jbrowse-bam-coverage', {'bam': bam.pk})
        self.assertFile(coverage, 'bigwig_track', 'Jbrowse_genome_coverage.bw')

    def test_norm_coverage_track(self):
        inputs = {"src": "alignment_coverage.bam"}
        bam = self.run_process("upload-bam", inputs)

        inputs = {'bam': bam.pk, 'size': 34000000}
        coverage = self.run_process('jbrowse-bam-coverage-normalized', inputs)
        self.assertFile(coverage, 'bigwig_track', 'Jbrowse_norm_genome_coverage.bw')
