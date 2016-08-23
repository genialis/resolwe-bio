# pylint: disable=missing-docstring
from resolwe.flow.models import Data

from resolwe_bio.utils.test import BioProcessTestCase


class AbyssProcessorTestCase(BioProcessTestCase):

    def test_abyss(self):
        se_reads = self.prepare_reads(['reads.fastq.gz'])

        inputs = {'src1': ['reads_paired_abyss_1.fastq.gz'], 'src2': ['reads_paired_abyss_2.fastq.gz']}
        reads = self.run_process('upload-fastq-paired', inputs)

        inputs = {'paired_end': reads.pk,
                  'se': se_reads.pk,
                  'options': {'k': 25,
                              'name': 'ecoli'}}

        self.run_process('assembler-abyss', inputs)

        inputs = {'paired_end': reads.pk,
                  'se': se_reads.pk,
                  'use_unmapped': True,
                  'options': {'k': 25,
                              'name': 'ecoli'}}

        self.run_process('assembler-abyss', inputs, Data.STATUS_ERROR)
