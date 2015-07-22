# pylint: disable=missing-docstring
from .base import BaseProcessorTestCase
from .utils import PreparedData


class CoverageProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_refseq_track(self):
        genome = self.prepare_genome()
        refseq_track = self.run_processor('jbrowse:refseq', {'refseq': genome.pk})
        self.assertFiles(refseq_track, 'refseq_track', 'refseq.json')

    def test_gff3_track(self):
        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        gff = self.run_processor('jbrowse:gff3', {'gff': annotation.pk})
        self.assertFields(gff, 'annotation_track.refs', ['tracks/annotation'])

    def test_gtf_track(self):
        annotation = self.prepare_annotation()
        gtf = self.run_processor('jbrowse:gtf', {'gtf': annotation.pk})
        self.assertFields(gtf, 'annotation_track.refs', ['tracks/annotation'])

    def test_coverage_track(self):
        inputs = {"src": "alignment_name_sorted.bam"}
        bam = self.run_processor("import:upload:mapping-bam", inputs)

        coverage = self.run_processor('jbrowse:bam:coverage', {'bam': bam.pk})
        self.assertFiles(coverage, 'bigwig_track', 'Jbrowse_genome_coverage.bw')
