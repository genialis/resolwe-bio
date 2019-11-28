# pylint: disable=missing-docstring
import os

from resolwe.flow.models import Process
from resolwe.test import tag_process
from resolwe_bio.utils.test import BioProcessTestCase


class SlamdunkProcessorTestCase(BioProcessTestCase):

    @tag_process('slamdunk-all-paired')
    def test_slamdunk_paired(self):
        with self.preparation_stage():
            paired_reads = self.prepare_paired_reads(['hs_slamseq_R1_complemented.fastq.gz'],
                                                     ['hs_slamseq_R2.fastq.gz'])
            transcripts = self.run_process('upload-fasta-nucl', {
                'src': os.path.join('slamseq', 'input', 'hs_transcript.fasta'),
                'species': 'Homo sapiens',
                'build': 'Gencode 32'
            })
            bedfile = self.run_process('upload-bed', {
                'src': os.path.join('slamseq', 'input', 'hs_transcript.bed'),
                'species': 'Homo sapiens',
                'build': 'Gencode 32'
            })
        inputs = {
            'reads': paired_reads.id,
            'ref_seq': transcripts.id,
            'regions': bedfile.id,
            'filter_multimappers': True,
            'max_alignments': 1,
            'read_length': 75
        }
        slamdunk = self.run_process('slamdunk-all-paired', inputs)
        self.assertFile(slamdunk, 'tcount', os.path.join('slamseq', 'output', 'hs_slamseq_tcount.tsv'))

    @tag_process('alleyoop-rates')
    def test_alleyoop_rates(self):
        with self.preparation_stage():
            reference = self.run_process('upload-fasta-nucl', {
                'src': os.path.join('slamseq', 'input', 'hs_transcript.fasta'),
                'species': 'Homo sapiens',
                'build': 'Gencode 32'
            })
            process = Process.objects.create(
                name='Upload Slamdunk data mock process',
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
                type='data:alignment:bam:slamdunk:',
                input_schema=[
                    {
                        'name': 'src',
                        'type': 'basic:file:',
                    },
                    {
                        'name': 'index',
                        'type': 'basic:file:',
                    },
                ],
                output_schema=[
                    {
                        'name': 'bam',
                        'type': 'basic:file:',
                    },
                    {
                        'name': 'bai',
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
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "bam" "bam" 0.1 extract
re-save-file bam "${NAME}".bam

re-import {{ index.file_temp|default(index.file) }} {{ index.file }} "bai" "bai" 0.1 extract
re-save-file bai "${NAME}".bai
re-save species "Homo sapiens"
re-save build "Gencode 32"
"""
                }
            )
            slamdunk = self.run_process(process.slug, {
                'src': os.path.join('slamseq', 'output', 'hs_slamseq_slamdunk_mapped_filtered.bam'),
                'index': os.path.join('slamseq', 'output', 'hs_slamseq_slamdunk_mapped_filtered.bam.bai')
            })

        alleyoop_rates = self.run_process('alleyoop-rates', {
            'ref_seq': reference.id,
            'slamdunk': slamdunk.id
        })
        self.assertFile(alleyoop_rates, 'report', os.path.join('slamseq', 'output', 'hs_alleyoop_overallrates.txt'))

    @tag_process('alleyoop-utr-rates')
    def test_alleyoop_utrrates(self):
        with self.preparation_stage():
            reference = self.run_process('upload-fasta-nucl', {
                'src': os.path.join('slamseq', 'input', 'hs_transcript.fasta'),
                'species': 'Homo sapiens',
                'build': 'Gencode 32'
            })
            bedfile = self.run_process('upload-bed', {
                'src': os.path.join('slamseq', 'input', 'hs_transcript.bed'),
                'species': 'Homo sapiens',
                'build': 'Gencode 32'
            })
            process = Process.objects.create(
                name='Upload Slamdunk data mock process',
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
                type='data:alignment:bam:slamdunk:',
                input_schema=[
                    {
                        'name': 'src',
                        'type': 'basic:file:',
                    },
                    {
                        'name': 'index',
                        'type': 'basic:file:',
                    },
                ],
                output_schema=[
                    {
                        'name': 'bam',
                        'type': 'basic:file:',
                    },
                    {
                        'name': 'bai',
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
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "bam" "bam" 0.1 extract
re-save-file bam "${NAME}".bam

re-import {{ index.file_temp|default(index.file) }} {{ index.file }} "bai" "bai" 0.1 extract
re-save-file bai "${NAME}".bai
re-save species "Homo sapiens"
re-save build "Gencode 32"
"""
                }
            )
            slamdunk = self.run_process(process.slug, {
                'src': os.path.join('slamseq', 'output', 'hs_slamseq_slamdunk_mapped_filtered.bam'),
                'index': os.path.join('slamseq', 'output', 'hs_slamseq_slamdunk_mapped_filtered.bam.bai')
            })
        rates = self.run_process('alleyoop-utr-rates', {
            'ref_seq': reference.id,
            'regions': bedfile.id,
            'slamdunk': slamdunk.id
        })
        self.assertFile(
            rates,
            'report',
            os.path.join('slamseq', 'output', 'hs_alleyoop_mutationrates.txt')
        )
