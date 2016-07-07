# pylint: disable=missing-docstring
from django.core.management import call_command

from resolwe.flow.models import Data, Storage
from resolwe_bio.utils.test import BioProcessTestCase

class GenerateEtcTest(BioProcessTestCase):

    def setUp(self):
        super(GenerateEtcTest, self).setUp()

    def test_generate_samples(self):
        call_command('generate_etc', '-e=1', '--rseed')
        data = Data.objects.last()
        storage = Storage.objects.last()
        self.assertEqual(data.slug, 'etc')
        self.assertEqual(data.name, 'D. discoideum')
        self.assertEqual(data.descriptor['parental_strain'], 'AX4')
        self.assertEqual(data.output['etcfile']['file'], 'etc.json.gz')
        self.assertEqual(storage.json['etc']['timePoints'], [0, 4, 8, 12, 16, 20, 24])
