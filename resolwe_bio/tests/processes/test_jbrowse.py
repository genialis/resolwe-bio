# pylint: disable=missing-docstring
from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase


class JbrowseProcessorTestCase(BioProcessTestCase):

    @tag_process('jbrowse-refseq')
    def test_refseq_track(self):
        with self.preparation_stage():
            genome = self.prepare_genome()

        refseq_track = self.run_process('jbrowse-refseq', {'refseq': genome.pk})
        del refseq_track.output['refseq_track']['total_size']  # Non-deterministic output.
        self.assertFields(refseq_track, 'refseq_track', {'refs': ['seq'],
                                                         'file': 'seq/refSeqs.json'})
        self.assertFields(refseq_track, 'species', 'Dictyostelium discoideum')
        self.assertFields(refseq_track, 'build', 'dd-05-2009')

    @tag_process('jbrowse-gff3')
    def test_gff3_track(self):
        with self.preparation_stage():
            annotation = self.prepare_annotation_gff()

        gff = self.run_process('jbrowse-gff3', {'gff': annotation.pk})
        del gff.output['annotation_track']['total_size']  # Non-deterministic output.

        self.assertFields(gff, 'annotation_track', {'refs': ['tracks/annotation'],
                                                    'file': 'trackList.json'})
        self.assertFields(gff, 'species', 'Dictyostelium discoideum')
        self.assertFields(gff, 'build', 'dd-05-2009')

    @tag_process('jbrowse-gtf')
    def test_gtf_track(self):
        with self.preparation_stage():
            annotation = self.prepare_annotation()

        gtf = self.run_process('jbrowse-gtf', {'gtf': annotation.pk})
        del gtf.output['annotation_track']['total_size']  # Non-deterministic output.
        self.assertFields(gtf, 'annotation_track', {'refs': ['tracks/annotation'],
                                                    'file': 'trackList.json'})
        self.assertFields(gtf, 'species', 'Dictyostelium discoideum')
        self.assertFields(gtf, 'build', 'dd-05-2009')

    @tag_process('jbrowse-bed')
    def test_bed_track(self):
        with self.preparation_stage():
            inputs = {
                'src': 'bed_track.bed',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009'
            }
            bed_file = self.run_process('upload-bed', inputs)

        bed = self.run_process('jbrowse-bed', {'bed': bed_file.pk})
        del bed.output['bed_track']['total_size']  # Non-deterministic output.

        self.assertFields(bed, 'bed_track', {'refs': ['tracks/bed'],
                                             'file': 'trackList.json'})
        self.assertFields(bed, 'species', 'Dictyostelium discoideum')
        self.assertFields(bed, 'build', 'dd-05-2009')

    @tag_process('jbrowse-bam-coverage')
    def test_coverage_track(self):
        with self.preparation_stage():
            inputs = {
                'src': 'alignment_coverage.bam',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            bam = self.run_process('upload-bam', inputs)

        coverage = self.run_process('jbrowse-bam-coverage', {'bam': bam.pk})
        self.assertFile(coverage, 'bigwig_track', 'Jbrowse_genome_coverage.bw')
        self.assertFields(coverage, 'species', 'Homo sapiens')
        self.assertFields(coverage, 'build', 'hg19')

    @tag_process('jbrowse-bam-coverage-normalized')
    def test_norm_coverage_track(self):
        with self.preparation_stage():
            inputs = {
                'src': 'alignment_coverage.bam',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            bam = self.run_process('upload-bam', inputs)

        inputs = {'bam': bam.pk, 'size': 34000000}
        coverage = self.run_process('jbrowse-bam-coverage-normalized', inputs)
        self.assertFile(coverage, 'bigwig_track', 'Jbrowse_norm_genome_coverage.bw')
        self.assertFields(coverage, 'species', 'Homo sapiens')
        self.assertFields(coverage, 'build', 'hg19')
