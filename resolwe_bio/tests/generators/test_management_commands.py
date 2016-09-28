# pylint: disable=missing-docstring
from os.path import join

import six

from django.core.management import call_command

from resolwe.flow.models import Data
from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles
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
            expr = sample_data[2]
            # NOTE: Python 2 and 3 produce different results even when setting random.seed() to the same number
            if six.PY2:
                self.assertJSON(expr, expr.output['exp_json'], '', 'expressions-py2.json.gz')
                self.assertFile(expr, 'exp', 'expressions-py2.tab.gz', compression='gzip')
            else:
                self.assertJSON(expr, expr.output['exp_json'], '', 'expressions-py3.json.gz')
                self.assertFile(expr, 'exp', 'expressions-py3.tab.gz', compression='gzip')
        else:
            self.fail("Sample not created")


class GenerateEtcTest(BioProcessTestCase):

    @skipUnlessLargeFiles('etc-py2.json.gz', 'etc-py3.json.gz')
    def test_generate_etc(self):
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

    @skipUnlessLargeFiles('DE-cuffdiff-py2.json.gz', 'DE-cuffdiff-py3.json.gz',
                          'DE-cuffdiff-py2.tab.gz', 'DE-cuffdiff-py3.tab.gz')
    def test_generate_diffexpr_cuffdiff(self):
        call_command('generate_diffexpr_cuffdiff', '-n=1', '-g=2', '--rseed')
        diffexpr = Data.objects.last()
        if diffexpr:
            self.assertEqual(diffexpr.process.type, 'data:differentialexpression:cuffdiff:')
            self.assertEqual(len(diffexpr.input['case']), 2)
            self.assertEqual(len(diffexpr.input['control']), 2)
            # NOTE: Python 2 and 3 produce different results even when setting random.seed() to the same number
            if six.PY2:
                self.assertJSON(diffexpr, diffexpr.output['de_json'], '',
                                join('large', 'DE-cuffdiff-py2.json.gz'))
                self.assertFile(diffexpr, 'de_file',
                                join('large', 'DE-cuffdiff-py2.tab.gz'), compression='gzip')
            else:
                self.assertJSON(diffexpr, diffexpr.output['de_json'], '',
                                join('large', 'DE-cuffdiff-py3.json.gz'))
                self.assertFile(diffexpr, 'de_file',
                                join('large', 'DE-cuffdiff-py3.tab.gz'), compression='gzip')

            expressions = diffexpr.input['case'] + diffexpr.input['control']
            for expr_id in expressions:
                expression = Data.objects.get(id=expr_id)
                self.assertEqual(expression.process.type, 'data:cufflinks:cuffquant:')
        else:
            self.fail("Differential expression not created")

    @skipUnlessLargeFiles('DE-deseq-py2.json.gz', 'DE-deseq-py3.json.gz',
                          'DE-deseq-py2.tab.gz', 'DE-deseq-py3.tab.gz')
    def test_generate_diffexpr_deseq(self):
        call_command('generate_diffexpr_deseq', '-n=1', '-g=2', '--rseed')
        diffexpr = Data.objects.last()
        if diffexpr:
            self.assertEqual(diffexpr.process.type, 'data:differentialexpression:deseq2:')
            self.assertEqual(len(diffexpr.input['case']), 2)
            self.assertEqual(len(diffexpr.input['control']), 2)
            # NOTE: Python 2 and 3 produce different results even when setting random.seed() to the same number
            if six.PY2:
                self.assertJSON(diffexpr, diffexpr.output['de_json'], '',
                                join('large', 'DE-deseq-py2.json.gz'))
                self.assertFile(diffexpr, 'de_file',
                                join('large', 'DE-deseq-py2.tab.gz'), compression='gzip')
            else:
                self.assertJSON(diffexpr, diffexpr.output['de_json'], '',
                                join('large', 'DE-deseq-py3.json.gz'))
                self.assertFile(diffexpr, 'de_file',
                                join('large', 'DE-deseq-py3.tab.gz'), compression='gzip')

            expressions = diffexpr.input['case'] + diffexpr.input['control']
            for expr_id in expressions:
                expression = Data.objects.get(id=expr_id)
                self.assertEqual(expression.process.type, 'data:expression:')
        else:
            self.fail("Differential expression not created")


class GenerateGeneSetTest(BioProcessTestCase):

    def test_generate_geneset(self):
        call_command('generate_geneset', '-n=1', '--rseed')
        geneset = Data.objects.last()
        if geneset:
            self.assertEqual(geneset.process.type, 'data:geneset:')
            self.assertEqual(geneset.output['source'], 'mm10')
            # NOTE: Python 2 and 3 produce different results even when setting random.seed() to the same number
            if six.PY2:
                self.assertJSON(geneset, geneset.output['geneset_json'], '', 'geneset-py2.json.gz')
                self.assertFile(geneset, 'geneset', 'geneset-py2.tab.gz', compression='gzip')
            else:
                self.assertJSON(geneset, geneset.output['geneset_json'], '', 'geneset-py3.json.gz')
                self.assertFile(geneset, 'geneset', 'geneset-py3.tab.gz', compression='gzip')
        else:
            self.fail("Gene set not created")
