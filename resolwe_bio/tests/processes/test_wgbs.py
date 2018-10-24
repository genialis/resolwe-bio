# pylint: disable=missing-docstring
from resolwe.flow.models import Process
from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase


class WgbsProcessorTestCase(BioProcessTestCase):

    @tag_process('walt')
    def test_walt(self):
        with self.preparation_stage():
            inputs = {
                'src': 'hg19_chr2 17k.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'hg19',
            }
            genome = self.run_process('upload-genome', inputs)
            reads_paired = self.prepare_paired_reads(mate1=['3A_WT_WGBS_chr2_17k R1.fastq.gz'],
                                                     mate2=['3A_WT_WGBS_chr2_17k_R2.fastq.gz'])

        inputs = {
            'genome': genome.id,
            'reads': reads_paired.id,
        }
        walt = self.run_process('walt', inputs)
        self.assertFile(walt, 'stats', 'walt_report.txt')
        self.assertFileExists(walt, 'mr')
        self.assertFileExists(walt, 'unmapped_f')
        self.assertFileExists(walt, 'unmapped_r')
        self.assertFields(walt, 'species', 'Homo sapiens')
        self.assertFields(walt, 'build', 'hg19')

    @tag_process('methcounts')
    def test_methcounts(self):
        with self.preparation_stage():
            inputs = {
                'src': 'hg19_chr2 17k.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'hg19',
            }
            genome = self.run_process('upload-genome', inputs)

            # Mock upload WALT alignment process
            process = Process.objects.create(
                name='Upload WALT alignment file (.mr) mock process',
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
                type='data:alignment:mr:walt',
                input_schema=[
                    {
                        'name': 'fc',
                        'type': 'basic:file:',
                    },
                ],
                output_schema=[
                    {
                        'name': 'mr',
                        'type': 'basic:file:',
                    },
                    {
                        'name': 'species',
                        'type': 'basic:string:',
                    },
                    {
                        'name': 'build',
                        'type': 'basic:string:',
                    },
                ],
                run={
                    'language': 'bash',
                    'program': r"""
re-import {{ fc.file_temp|default(fc.file) }} {{ fc.file }} "mr" "mr" 0.1 compress
re-save-file mr "${NAME}".mr.gz
re-save species 'Homo sapiens'
re-save build 'hg19'
"""
                }
            )

            inputs = {'fc': '3A_WT_WGBS_chr2 17k.mr.gz'}
            walt = self.run_process(process.slug, inputs)

        inputs = {
            'genome': genome.id,
            'alignment': walt.id,
        }
        methcounts = self.run_process('methcounts', inputs)
        self.assertFile(methcounts, 'stats', 'methconts_report.txt')
        self.assertFileExists(methcounts, 'meth')
        self.assertFileExists(methcounts, 'bigwig')
        self.assertFields(methcounts, 'species', 'Homo sapiens')
        self.assertFields(methcounts, 'build', 'hg19')

    @tag_process('hmr')
    def test_hmr(self):
        with self.preparation_stage():
            # Mock upload methcounts process
            process = Process.objects.create(
                name='Upload methcounts file (.meth) mock process',
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
                type='data:wgbs:methcounts:',
                input_schema=[
                    {
                        'name': 'fc',
                        'type': 'basic:file:',
                    },
                ],
                output_schema=[
                    {
                        'name': 'meth',
                        'type': 'basic:file:',
                    },
                    {
                        'name': 'species',
                        'type': 'basic:string:',
                    },
                    {
                        'name': 'build',
                        'type': 'basic:string:',
                    },
                ],
                run={
                    'language': 'bash',
                    'program': r"""
re-import {{ fc.file_temp|default(fc.file) }} {{ fc.file }} "meth" "meth" 0.1 compress
re-save-file meth "${NAME}".meth.gz
re-save species 'Homo sapiens'
re-save build 'hg19'
"""
                }
            )

            inputs = {'fc': '3A_WT_WGBS_chr2 17k.meth.gz'}
            methcounts = self.run_process(process.slug, inputs)

        hmr = self.run_process('hmr', {'methcounts': methcounts.id})
        self.assertFile(hmr, 'hmr', '3A_WT_WGBS_chr2_17k.hmr.gz', compression='gzip')
        self.assertFields(hmr, 'species', 'Homo sapiens')
        self.assertFields(hmr, 'build', 'hg19')
