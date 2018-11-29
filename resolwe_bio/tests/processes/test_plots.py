# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase


class PlotsProcessorTestCase(BioProcessTestCase):

    @tag_process('bamplot')
    def test_bamplot(self):
        with self.preparation_stage():
            inputs = {
                'src': 'bamplot_alignment.bam',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            bam = self.run_process('upload-bam', inputs)

            inputs = {
                'src': 'bamplot_alignment.bam',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            bam1 = self.run_process('upload-bam', inputs)

            bed_input = {
                'src': 'bamplot.bed',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            bed = self.run_process('upload-bed', bed_input)

            inputs = {
                'src': 'bamplot.gff',
                'source': 'NCBI',
                'species': 'Homo sapiens',
                'build': 'GRCh38'
            }
            gff = self.run_process('upload-gtf', inputs)

        inputs = {
            'genome': 'HG19',
            'input_region': 'chr1:+:41468594-41566948',
            'bam': [bam.pk, bam1.pk],
            'color': '0,69,134',
            'names': ['WNT', 'bbb'],
            'yscale': 'uniform',
            'title': 'SINGLE_REGION',
            'plot': 'multiple',
            'rpm': True,
        }
        self.run_process('bamplot', inputs)

        inputs = {
            'genome': 'HG19',
            'input_gff': gff.pk,
            'bam': [bam.pk],
            'color': '255,192,0',
            'names': ['GROUP3_MB'],
            'yscale': 'uniform',
            'title': 'SINGLE_REGION',
            'plot': 'multiple',
            'rpm': True,
            'bed': [bed.pk],
        }
        self.run_process('bamplot', inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

    @tag_process('bamliquidator')
    def test_bamliquidator(self):
        with self.preparation_stage():
            inputs = {
                'src': 'bamplot_ alignment1.bam',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            bam1 = self.run_process('upload-bam', inputs)

            inputs = {
                'src': 'bamplot_alignment.bam',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            bam = self.run_process('upload-bam', inputs)

            inputs = {
                'src': 'bamplot.gff',
                'source': 'NCBI',
                'species': 'Homo sapiens',
                'build': 'GRCh38'
            }
            gff = self.run_process('upload-gtf', inputs)

        inputs = {'bam': [bam1.id, bam.id],
                  'cell_type': 'MCD cell',
                  'extension': 200}

        bamliquidator = self.run_process('bamliquidator', inputs)
        del bamliquidator.output['summary']['total_size']  # Non-deterministic output.
        self.assertFields(bamliquidator, 'summary', {'file': 'output/summary.html',
                                                     'size': 524296})

        inputs['regions_gtf'] = gff.id
        inputs['analysis_type'] = 'region'
        self.run_process('bamliquidator', inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
