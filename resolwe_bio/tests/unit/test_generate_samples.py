# pylint: disable=missing-docstring
from django.core.management import call_command

from resolwe.flow.tests import ProcessTestCase
from resolwe_bio.models import Sample


class GenerateSamplesTest(ProcessTestCase):

    def setUp(self):
        super(GenerateSamplesTest, self).setUp()

    def test_generate_samples(self):
        call_command('generate_samples', '-s=1', '-p=0', '--rseed')
        sample = Sample.objects.last()
        if sample:
            sample_data = sample.data.all()
            self.assertEqual(sample_data[0].slug, 'gs-reads')
            self.assertEqual(sample_data[1].slug, 'mapping')
            self.assertEqual(sample_data[2].slug, 'expression')
        else:
            self.fail("Sample not crated")
