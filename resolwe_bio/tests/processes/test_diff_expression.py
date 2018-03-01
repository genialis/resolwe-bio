# pylint: disable=missing-docstring
from resolwe.flow.models import Data, Process
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class DiffExpProcessorTestCase(BioProcessTestCase):

    @tag_process('cuffdiff')
    def test_cuffdiff(self):
        with self.preparation_stage():
            inputs = {
                'src': 'cuffquant_1.cxb',
                'source': 'UCSC',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            cuffquant = self.run_process("upload-cxb", inputs)

            inputs = {
                'src': 'cuffquant_2.cxb',
                'source': 'UCSC',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            cuffquant2 = self.run_process("upload-cxb", inputs)

            annotation = self.prepare_annotation(fn='hg19_chr20_small.gtf.gz', source='UCSC')

        inputs = {
            'case': [cuffquant.id],
            'control': [cuffquant2.id],
            'annotation': annotation.id}
        cuffdiff = self.run_process('cuffdiff', inputs)
        self.assertFile(cuffdiff, 'raw', 'raw_cuffdiff.tab.gz', compression='gzip')
        self.assertFile(cuffdiff, 'de_file', 'de_file_cuffdiff.tab.gz', compression='gzip')
        self.assertJSON(cuffdiff, cuffdiff.output['de_json'], '', 'cuffdiff.json.gz')
        self.assertFields(cuffdiff, 'source', 'UCSC')
        self.assertFields(cuffdiff, 'species', 'Homo sapiens')
        self.assertFields(cuffdiff, 'build', 'hg19')
        self.assertFields(cuffdiff, 'feature_type', 'gene')

    @tag_process('differentialexpression-deseq2')
    def test_deseq2_genes(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', source='DICTYBASE')
            expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', source='DICTYBASE')
            expression_3 = self.prepare_expression(f_rc='exp_3_rc.tab.gz', source='DICTYBASE')
            expression_4 = self.prepare_expression(f_rc='exp_4_rc.tab.gz', source='DICTYBASE')

        inputs = {
            'case': [expression_1.pk, expression_3.pk],
            'control': [expression_2.pk, expression_4.pk],
            'filter': 0,
        }

        diff_exp = self.run_process('differentialexpression-deseq2', inputs)

        self.assertFileExists(diff_exp, 'raw')
        self.assertJSON(diff_exp, diff_exp.output['de_json'], '', 'deseq2.json.gz')
        self.assertFields(diff_exp, 'source', 'DICTYBASE')
        self.assertFields(diff_exp, 'species', 'Dictyostelium discoideum')
        self.assertFields(diff_exp, 'build', 'dd-05-2009')
        self.assertFields(diff_exp, 'feature_type', 'gene')

    @tag_process('differentialexpression-deseq2')
    def test_deseq2_source(self):
        with self.preparation_stage():
            expression_dictybase = self.prepare_expression(source='DICTYBASE')
            expression_ucsc = self.prepare_expression(source='UCSC')

        inputs = {
            'case': [expression_dictybase.pk],
            'control': [expression_ucsc.pk]
        }

        self.run_process('differentialexpression-deseq2', inputs, Data.STATUS_ERROR)

    @tag_process('differentialexpression-edger')
    def test_edger(self):
        with self.preparation_stage():
            inputs = {
                'rc': 'exp_1_rc.tab.gz',
                'exp_name': 'Expression',
                'source': 'DICTYBASE',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009',
            }
            expression_1 = self.run_process('upload-expression', inputs)

            inputs = {
                'rc': 'exp_2_rc.tab.gz',
                'exp_name': 'Expression',
                'source': 'DICTYBASE',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009',
            }
            expression_2 = self.run_process('upload-expression', inputs)

            inputs = {
                'rc': 'exp_3_rc.tab.gz',
                'exp_name': 'Expression',
                'source': 'DICTYBASE',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009',
            }
            expression_3 = self.run_process('upload-expression', inputs)

            inputs = {
                'rc': 'exp_4_rc.tab.gz',
                'exp_name': 'Expression',
                'source': 'DICTYBASE',
                'species': 'Dictyostelium discoideum',
                'build': 'dd-05-2009',
            }
            expression_4 = self.run_process('upload-expression', inputs)

        inputs = {
            'case': [expression_1.id, expression_3.id],
            'control': [expression_2.id, expression_4.id]
        }

        diff_exp = self.run_process('differentialexpression-edger', inputs)
        self.assertFile(diff_exp, 'raw', 'diffexp_edgeR.tab.gz', compression='gzip')
        self.assertJSON(diff_exp, diff_exp.output['de_json'], '', 'edgeR.json.gz')
        self.assertFields(diff_exp, 'source', 'DICTYBASE')
        self.assertFields(diff_exp, 'species', 'Dictyostelium discoideum')
        self.assertFields(diff_exp, 'build', 'dd-05-2009')
        self.assertFields(diff_exp, 'feature_type', 'gene')

    @tag_process('differentialexpression-dexseq')
    def test_dexseq(self):
        with self.preparation_stage():
            # Mock upload featureCounts process
            process = Process.objects.create(
                name='Upload featureCounts mock process',
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
                data_name='{{ fc.file|default("?") }}',
                type='data:expression:featurecounts:',
                flow_collection="sample",
                input_schema=[
                    {
                        'name': 'fc',
                        'type': 'basic:file:',
                    },
                ],
                output_schema=[
                    {
                        'name': 'feature_counts_output',
                        'type': 'basic:file:',
                    },
                    {
                        'name': 'source',
                        'type': 'basic:string:',
                    },
                    {
                        'name': 'species',
                        'type': 'basic:string:',
                    },
                    {
                        'name': 'build',
                        'type': 'basic:string:',
                    },
                    {
                        'name': 'feature_type',
                        'type': 'basic:string:',
                    },
                ],
                run={
                    'language': 'bash',
                    'program': r"""
re-import {{ fc.file_temp|default(fc.file) }} {{ fc.file }} "txt" "txt" 0.1 compress
re-save-file feature_counts_output "${NAME}".txt.gz
re-save source 'ENSEMBL'
re-save species 'Mus musculus'
re-save build 'GRCm38.p5'
re-save feature_type 'exon'
"""
                }
            )

            inputs = {'fc': 'mus_musculus_short_fCount1.txt.gz'}
            featurecounts_1 = self.run_process(process.slug, inputs)

            inputs['fc'] = 'mus_musculus_short_fCount2.txt.gz'
            featurecounts_2 = self.run_process(process.slug, inputs)

            inputs['fc'] = 'mus_musculus_short_fCount3.txt.gz'
            featurecounts_3 = self.run_process(process.slug, inputs)

            inputs['fc'] = 'mus_musculus_short_fCount4.txt.gz'
            featurecounts_4 = self.run_process(process.slug, inputs)

            annotation = self.prepare_annotation(
                fn='annotation_ensembl_mus_musculus_short_dexseq.gtf.gz',
                source='ENSEMBL',
                species='Mus musculus',
                build='GRCm38.p5',
            )

        inputs = {
            'case': [featurecounts_1.pk, featurecounts_2.pk],
            'control': [featurecounts_3.pk, featurecounts_4.pk],
            'annotation': annotation.pk,
        }

        diff_exp = self.run_process('differentialexpression-dexseq', inputs)

        self.assertFileExists(diff_exp, 'raw')
        self.assertJSON(diff_exp, diff_exp.output['de_json'], '', 'dexseq.json.gz')
        self.assertFileExists(diff_exp, 'dxr')
        self.assertFields(diff_exp, 'source', 'ENSEMBL')
        self.assertFields(diff_exp, 'species', 'Mus musculus')
        self.assertFields(diff_exp, 'build', 'GRCm38.p5')
        self.assertFields(diff_exp, 'feature_type', 'exon')
