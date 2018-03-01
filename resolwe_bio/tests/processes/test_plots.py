# pylint: disable=missing-docstring
from resolwe.flow.models import Process
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

    @tag_process('bamplot')
    def test_bamplot_gff(self):
        with self.preparation_stage():
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
                'src': 'bamplot_alignment.bam',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
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

        inputs = {'bam': [bam1.id, bam.id],
                  'cell_type': 'MCD cell',
                  'extension': 200}

        bamliquidator = self.run_process('bamliquidator', inputs)
        del bamliquidator.output['summary']['total_size']  # Non-deterministic output.
        self.assertFields(bamliquidator, 'summary', {'file': 'output/summary.html',
                                                     'size': 524296})

    @tag_process('bamliquidator')
    def test_bamliquidator_gff(self):
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
                  'analysis_type': 'region',
                  'cell_type': 'MCD cell',
                  'regions_gtf': gff.id,
                  'extension': 200, }

        self.run_process('bamliquidator', inputs)

    @tag_process('plot-dexeq')
    def test_dexseq_plot(self):
        with self.preparation_stage():
            # Mock upload dexseq process
            process = Process.objects.create(
                name='Upload dexseq mock process',
                requirements={
                    'expression-engine': 'jinja',
                    'resources': {
                        'network': True,
                    },
                    'executor': {
                        'docker': {
                            'image': 'resolwebio/base:ubuntu-18.04',
                        },
                    },
                },
                contributor=self.contributor,
                type='data:differentialexpression:dexseq:',
                input_schema=[
                    {
                        'name': 'dxr',
                        'type': 'basic:file:',
                    },
                ],
                output_schema=[
                    {
                        'name': 'dxr',
                        'type': 'basic:file:',
                    },
                ],
                run={
                    'language': 'bash',
                    'program': r"""
re-import {{ dxr.file_temp|default(dxr.file) }} {{ dxr.file }} "RData" "RData" 0.1
re-save-file dxr "${NAME}".RData
"""
                }
            )

            inputs = {
                'dxr': 'dxr.RData',
            }
            dexseq = self.run_process(process.slug, inputs)

        inputs = {
            'dexseq': dexseq.pk,
            'gene': 'ENSMUSG00000024831',
            'transcript': True,
            'sample': True,
        }

        dexseq_plot = self.run_process('plot-dexeq', inputs)

        self.assertFileExists(dexseq_plot, 'pdf')
