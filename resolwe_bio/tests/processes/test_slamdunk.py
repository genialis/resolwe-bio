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
