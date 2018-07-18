# pylint: disable=missing-docstring
from resolwe.test import tag_process
from django.conf import settings

from resolwe_bio.utils.test import BioProcessTestCase


class GenomeBrowserProcessorTestCase(BioProcessTestCase):

    @tag_process('igv')
    def test_igv_bam(self):
        with self.preparation_stage():
            inputs = {
                'src': 'reads.bam',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }

            bam = self.run_process('upload-bam', inputs)

        inputs = {
            'genomeid': 'hg19',
            'bam': [bam.id],
            'locus': 'chr7:79439229-79481604'
        }

        igv_session = self.run_process('igv', inputs)

        # Remove resource path lines that have changing data object ids.
        def filter_resource(line):
            if '<Resource path="{}/data/'.format(settings.RESOLWE_HOST_URL).encode('utf-8') in line:
                return True

        self.assertFile(igv_session, 'igv_file', 'igv_session_bam.xml', file_filter=filter_resource)
