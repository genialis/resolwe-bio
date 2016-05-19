# pylint: disable=missing-docstring
from .utils import skipDockerFailure, BioProcessTestCase


class JbrowseProcessorTestCase(BioProcessTestCase):

    @skipDockerFailure("Fails due to different order of objects in JSON file")
    def test_refseq_track(self):
        genome = self.prepare_genome()
        refseq_track = self.run_processor('jbrowse-refseq', {'refseq': genome.pk})
        self.assertFiles(refseq_track, 'refseq_track', 'refseq.json')

    def test_gff3_track(self):
        inputs = {'src': 'annotation.gff.gz'}
        annotation = self.run_processor('upload-gff3', inputs)
        gff = self.run_processor('jbrowse-gff3', {'gff': annotation.pk})
        self.assertFields(gff, 'annotation_track.refs', ['tracks/annotation'])

    def test_gtf_track(self):
        annotation = self.prepare_annotation()
        gtf = self.run_processor('jbrowse-gtf', {'gtf': annotation.pk})
        self.assertFields(gtf, 'annotation_track.refs', ['tracks/annotation'])

    def test_bed_track(self):
        inputs = {'src': 'bed_track.bed'}
        bed_file = self.run_processor('upload-bed', inputs)
        bed = self.run_processor('jbrowse-bed', {'bed': bed_file.pk})
        self.assertFields(bed, 'bed_track.refs', ['tracks/bed'])

    @skipDockerFailure("Fails with: bedGraphToBigWig: command not found")
    def test_coverage_track(self):
        inputs = {"src": "alignment_coverage.bam"}
        bam = self.run_processor("upload-bam", inputs)

        coverage = self.run_processor('jbrowse-bam-coverage', {'bam': bam.pk})
        self.assertFiles(coverage, 'bigwig_track', 'Jbrowse_genome_coverage.bw')

    @skipDockerFailure("Fails with: bamCoverage: command not found")
    def test_norm_coverage_track(self):
        inputs = {"src": "alignment_coverage.bam"}
        bam = self.run_processor("upload-bam", inputs)

        inputs = {'bam': bam.pk, 'size': 34000000}
        coverage = self.run_processor('jbrowse-bam-coverage-normalized', inputs)
        self.assertFiles(coverage, 'bigwig_track', 'Jbrowse_norm_genome_coverage.bw')
