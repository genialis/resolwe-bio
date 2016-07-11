# pylint: disable=missing-docstring
from django.core.management import call_command

from resolwe.flow.models import Data, Storage
from resolwe_bio.utils.test import BioProcessTestCase
from resolwe_bio.models import Sample


class GenerateSamplesTest(BioProcessTestCase):

    def test_generate_samples(self):
        call_command('generate_samples', '-s=1', '-p=0', '--rseed')
        sample = Sample.objects.last()
        if sample:
            sample_data = sample.data.all().order_by('id')
            self.assertEqual(sample_data[0].slug, 'gs-reads')
            self.assertEqual(sample_data[1].slug, 'mapping')
            self.assertEqual(sample_data[2].slug, 'expression')
        else:
            self.fail("Sample not created")


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


class GenerateDiffexprTest(BioProcessTestCase):

    def test_generate_diffexpr(self):
        call_command('generate_diffexpr', '-n=1', '-g=2')
        diffexpr = Data.objects.last()
        if diffexpr:
            self.assertEqual(diffexpr.process.type, 'data:differentialexpression:cuffdiff:')
            self.assertEqual(len(diffexpr.input['case']), 2)
            self.assertEqual(len(diffexpr.input['control']), 2)

            expressions = diffexpr.input['case'] + diffexpr.input['control']

            for expr_id in expressions:
                expression = Data.objects.get(id=expr_id)
                self.assertEqual(expression.process.type, 'data:cufflinks:cuffquant:')
        else:
            self.fail("Differential expression not created")
