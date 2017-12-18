# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class DemultiplexProcessorTestCase(BioProcessTestCase):

    @tag_process('upload-multiplexed-paired')
    def test_demultiplex(self):
        inputs = {
            'reads': 'pool24.read1.small.qseq.bz2',
            'reads2': 'pool24.read3.small.qseq.bz2',
            'barcodes': 'pool24.read2.small.qseq.bz2',
            'annotation': 'pool24.tsv'
        }
        obj = self.run_process('upload-multiplexed-paired', inputs)
        self.assertFields(obj, 'matched', '5 reads (7.04 %)')
        self.assertFields(obj, 'notmatched', '51 reads (71.83 %)')
        self.assertFields(obj, 'badquality', '13 reads (18.31 %)')
        self.assertFields(obj, 'skipped', '0 reads (0.00 %)')

        reads = Data.objects.last()
        self.assertFiles(reads, "fastq", ["demultiplexed_reads_fw.fastq.gz"], compression='gzip')
        self.assertFiles(reads, "fastq2", ["demultiplexed_reads_rw.fastq.gz"], compression='gzip')

        inputs = {
            'reads': 'pool24.read3.small.qseq.bz2',
            'reads2': 'pool24.read1.small.qseq.bz2',
            'barcodes': 'pool24.read2.small.qseq.bz2',
            'annotation': 'pool24.tsv'
        }
        obj = self.run_process('upload-multiplexed-paired', inputs)
        self.assertFields(obj, 'matched', '5 reads (7.04 %)')
        self.assertFields(obj, 'notmatched', '51 reads (71.83 %)')
        self.assertFields(obj, 'badquality', '13 reads (18.31 %)')
        self.assertFields(obj, 'skipped', '0 reads (0.00 %)')

        inputs = {
            'reads': 'pool24.read1.small.qseq.bz2',
            'reads2': 'pool24.read3.small.qseq.bz2',
            'barcodes': 'pool24.read2.small.qseq.bz2',
            'annotation': 'pool24.2.tsv'
        }
        obj = self.run_process('upload-multiplexed-paired', inputs)
        self.assertFields(obj, 'matched', '5 reads (7.04 %)')
        self.assertFields(obj, 'notmatched', '51 reads (71.83 %)')
        self.assertFields(obj, 'badquality', '13 reads (18.31 %)')
        self.assertFields(obj, 'skipped', '0 reads (0.00 %)')

        inputs = {
            'reads': 's_1_1.qseq.small.txt.bz2',
            'reads2': 's_1_3.qseq.small.txt.bz2',
            'barcodes': 's_1_2.qseq.small.txt.bz2',
            'annotation': 's_1.tsv'
        }
        obj = self.run_process('upload-multiplexed-paired', inputs)
        self.assertFields(obj, 'matched', '12 reads (11.88 %)')
        self.assertFields(obj, 'notmatched', '80 reads (79.21 %)')
        self.assertFields(obj, 'badquality', '8 reads (7.92 %)')
        self.assertFields(obj, 'skipped', '0 reads (0.00 %)')
