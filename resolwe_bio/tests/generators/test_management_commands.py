# pylint: disable=missing-docstring
from os.path import join

import six

from django.core.management import call_command

from resolwe.flow.models import Data
from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles
from resolwe_bio.models import Sample


class GenerateSamplesTest(BioProcessTestCase):

    @skipUnlessLargeFiles('expressions-py2.json.gz', 'expressions-py3.json.gz',
                          'expressions-py2.tab.gz', 'expressions-py3.tab.gz')
    def test_generate_samples(self):
        call_command('generate_samples', '-s=1', '-p=0', '--rseed')
        sample = Sample.objects.last()
        if sample:
            sample_data = sample.data.all().order_by('id')
            self.assertEqual(sample_data[0].slug, 'gs-reads')
            self.assertEqual(sample_data[1].slug, 'mapping')
            self.assertEqual(sample_data[2].slug, 'expression')
            expr = sample_data[2]
            # NOTE: Python 2 and 3 produce different results even when setting random.seed() to the same number
            if six.PY2:
                self.assertJSON(expr, expr.output['exp_json'], '', join('large', 'expressions-py2.json.gz'))
                self.assertFile(expr, 'exp', join('large', 'expressions-py2.tab.gz'), compression='gzip')
            else:
                self.assertJSON(expr, expr.output['exp_json'], '', join('large', 'expressions-py3.json.gz'))
                self.assertFile(expr, 'exp', join('large', 'expressions-py3.tab.gz'), compression='gzip')
        else:
            self.fail("Sample not created")


class GenerateEtcTest(BioProcessTestCase):

    @skipUnlessLargeFiles('etc-py2.json.gz', 'etc-py3.json.gz')
    def test_generate_samples(self):
        call_command('generate_etc', '-e=1', '--rseed')
        data = Data.objects.last()
        self.assertEqual(data.slug, 'etc')
        self.assertEqual(data.name, 'D. discoideum')
        self.assertEqual(data.descriptor['parental_strain'], 'AX4')
        # NOTE: Python 2 and 3 produce different results even when setting random.seed() to the same number
        if six.PY2:
            self.assertJSON(data, data.output['etc'], '', join('large', 'etc-py2.json.gz'))
            self.assertFile(data, 'etcfile', join('large', 'etc-py2.json.gz'), compression='gzip')
        else:
            self.assertJSON(data, data.output['etc'], '', join('large', 'etc-py3.json.gz'))
            self.assertFile(data, 'etcfile', join('large', 'etc-py3.json.gz'), compression='gzip')


class GenerateDiffExprTest(BioProcessTestCase):

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
