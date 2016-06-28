# pylint: disable=missing-docstring
from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase


class CompatibilityProcessorTestCase(BioProcessTestCase):

    def test_reference_compatibility(self):
        inputs = {"src": "sp_test.fasta"}
        genome = self.run_processor('upload-genome', inputs)

        mapping = self.prepare_bam()
        annotation = self.prepare_annotation()

        inputs = {'reference': genome.pk, 'bam': mapping.pk, 'annot': annotation.pk}
        compatibility_test = self.run_processor('reference_compatibility', inputs)
        self.assertFile(compatibility_test, 'report_file', 'sp_test_compatibility_report.txt')

    def test_feature_location(self):
        inputs = {'src': 'mm10_small.gtf.gz'}
        annotation = self.run_processor('upload-gtf', inputs)

        inputs = {'annotation': annotation.pk,
                  'feature_type': 'exon',
                  'id_type': 'transcript_id',
                  'summarize_exons': True}
        features = self.run_processor('feature_location', inputs)
        self.assertJSON(features, features.output['feature_location'], '', 'feature_locations.json.gz')
