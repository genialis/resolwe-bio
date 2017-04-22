# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class PlotsProcessorTestCase(BioProcessTestCase):

    def test_bamplot(self):
        inputs = {'src': 'bamplot_alignment.bam'}
        bam = self.run_process('upload-bam', inputs)

        inputs = {'src': 'bamplot_alignment.bam'}
        bam1 = self.run_process('upload-bam', inputs)

        inputs = {'genome': 'HG19',
                  'input_region': 'chr1:+:41468594-41566948',
                  'bam': [bam.pk, bam1.pk],
                  'color': '0,69,134',
                  'names': ['WNT', 'bbb'],
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
                  'names': ['GROUP3_MB'],
                  'yscale': 'uniform',
                  'title': 'SINGLE_REGION',
                  'plot': 'multiple',
                  'rpm': True,
                  'bed': [bed.pk], }
        self.run_process('bamplot', inputs)

    def test_bamliquidator(self):
        inputs = {'src': 'bamplot_ alignment1.bam'}
        bam1 = self.run_process('upload-bam', inputs)

        inputs = {'src': 'bamplot_alignment.bam'}
        bam = self.run_process('upload-bam', inputs)

        inputs = {'bam': [bam1.id, bam.id],
                  'cell_type': 'MCD cell',
                  'extension': 200}

        bamliquidator = self.run_process('bamliquidator', inputs)
        self.assertFields(bamliquidator, 'summary', {'file': 'output/summary.html', 'size': 524296})

    def test_bamliquidator_gff(self):
        inputs = {'src': 'bamplot_ alignment1.bam'}
        bam1 = self.run_process('upload-bam', inputs)

        inputs = {'src': 'bamplot_alignment.bam'}
        bam = self.run_process('upload-bam', inputs)

        inputs = {'src': 'bamplot.gff', 'source': 'NCBI'}
        gff = self.run_process('upload-gtf', inputs)

        inputs = {'bam': [bam1.id, bam.id],
                  'analysis_type': 'region',
                  'cell_type': 'MCD cell',
                  'regions_gtf': gff.id,
                  'extension': 200, }

        self.run_process('bamliquidator', inputs)
