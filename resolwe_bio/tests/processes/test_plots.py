# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class PlotsProcessorTestCase(BioProcessTestCase):

    def test_bamplot(self):
        inputs = {'src': 'bamplot_alignment.bam'}
        bam = self.run_process('upload-bam', inputs)

        inputs = {'genome': 'HG19',
                  'input_region': 'chr1:+:41468594-41566948',
                  'bam': [bam.pk],
                  'color': '0,69,134',
                  'names': 'WNT',
                  'yscale': 'uniform',
                  'title': 'SINGLE_REGION',
                  'plot': 'multiple',
                  'rpm': True, }
        self.run_process('bamplot', inputs)

    def test_bamplot_gff(self):
        inputs = {'src': 'bamplot.bed'}
        bed = self.run_process('upload-bed', inputs)

        inputs = {'src': 'bamplot.gff', 'source': 'NCBI'}
        gff = self.run_process('upload-gtf', inputs)

        inputs = {'src': 'bamplot_alignment.bam'}
        bam = self.run_process('upload-bam', inputs)

        inputs = {'genome': 'HG19',
                  'input_gff': gff.pk,
                  'bam': [bam.pk],
                  'color': '255,192,0',
                  'names': 'GROUP3_MB',
                  'yscale': 'uniform',
                  'title': 'SINGLE_REGION',
                  'plot': 'multiple',
                  'rpm': True,
                  'bed': [bed.pk], }
        self.run_process('bamplot', inputs)
