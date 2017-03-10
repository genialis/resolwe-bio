# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class GenomeBrowserProcessorTestCase(BioProcessTestCase):

    def test_igv_bam(self):
        bam = self.prepare_bam()

        inputs = {'src': 'reads.bam'}
        bam1 = self.run_process('upload-bam', inputs)

        inputs = {
            'genomeid': 'hg19',
            'bam': [bam.pk, bam1.pk],
            'locus': 'chr7:79439229-79481604'
        }

        igv_session = self.run_process('igv', inputs)

        # remove changing lines from the output
        def filter_resource(line):
            if b'<Resource path=' in line:
                return True

        self.assertFile(igv_session, 'igv_file', 'igv_session_bam.xml', file_filter=filter_resource)
