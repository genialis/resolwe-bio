# pylint: disable=missing-docstring
from os.path import join

import mock
import six

from django.core.management import call_command

from resolwe.flow.models import Data
from resolwe_bio.models import Sample
from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles


class GeneratorBioProcessTestCase(BioProcessTestCase):

    def assertEqualIgnoreNumbers(self, first, second, msg=None):  # pylint: disable=invalid-name
        """Check that objects are equal, ignoring number differences.

        NOTE: Types of both objects (even if numbers) must be the same.

        :param first: first object
        :param second: second object
        :param str msg: custom error message on failure

        """
        if isinstance(first, six.integer_types) or isinstance(first, float):
            # only check that types are the same, ignore values
            self.assertIs(type(first), type(second))
        else:
            self.assertEqual(first, second, msg)

    def assertDictEqualIgnoreNumbers(self, dict1, dict2, msg=None):  # pylint: disable=invalid-name
        """Check that dictionaries are equal, ignoring number differences.

        :param dict d1: first dictionary
        :param dict d2: second dictionary
        :param str msg: custom error message on failure

        """
        self.assertIsInstance(dict1, dict, 'First argument is not a dictionary')
        self.assertIsInstance(dict2, dict, 'Second argument is not a dictionary')

        for k in dict1:
            self.assertIn(k, dict2)
            val1, val2 = dict1[k], dict2[k]
            if isinstance(val1, dict):
                self.assertDictEqualIgnoreNumbers(val1, val2)
            elif isinstance(val1, list) or isinstance(val1, tuple):
                self.assertEqual(len(val1), len(val2))
                for elem1, elem2 in zip(val1, val2):
                    self.assertEqualIgnoreNumbers(elem1, elem2)
            else:
                self.assertEqualIgnoreNumbers(val1, val2)


class GenerateSamplesTest(GeneratorBioProcessTestCase):

    def test_generate_samples(self):
        call_command('generate_samples', '-s=1', '-p=0', '--rseed')
        sample = Sample.objects.last()
        if sample:
            sample_data = sample.data.all().order_by('id')
            self.assertEqual(sample_data[0].slug, 'gs-reads')
            self.assertEqual(sample_data[1].slug, 'mapping')
            self.assertEqual(sample_data[2].slug, 'expression')

            expr = sample_data[2]
            # NOTE: The random.gammavariate() that is used for generating expression data produces
            # slightly differently rounded numbers on different platforms (e.g. Linux and MacOS) so
            # the check ignores the generated numbers
            with mock.patch.object(self, 'assertAlmostEqualGeneric', self.assertDictEqualIgnoreNumbers):
                # NOTE: Python 2 and 3 produce different results even when setting random.seed() to
                # the same number due to https://docs.python.org/3/whatsnew/3.2.html#random
                if six.PY2:
                    self.assertJSON(expr, expr.output['exp_json'], '', 'expressions-py2.json.gz')
                else:
                    self.assertJSON(expr, expr.output['exp_json'], '', 'expressions-py3.json.gz')
            self.assertFileExists(expr, 'exp')
        else:
            self.fail("Sample not created")


class GenerateEtcTest(GeneratorBioProcessTestCase):

    def test_generate_etc(self):
        call_command('generate_etc', '-e=1', '--rseed')
        data = Data.objects.last()
        self.assertEqual(data.slug, 'etc')
        self.assertEqual(data.name, 'D. discoideum')
        self.assertEqual(data.descriptor['parental_strain'], 'AX4')

        # NOTE: The random.gammavariate() that is used for generating expression data produces
        # slightly differently rounded numbers on different platforms (e.g. Linux and MacOS) so
        # the check ignores the generated numbers
        with mock.patch.object(self, 'assertAlmostEqualGeneric', self.assertDictEqualIgnoreNumbers):
            # NOTE: Python 2 and 3 produce different results even when setting random.seed() to the
            # same number due to https://docs.python.org/3/whatsnew/3.2.html#random
            if six.PY2:
                self.assertJSON(data, data.output['etc'], '', 'etc-py2.json.gz')
            else:
                self.assertJSON(data, data.output['etc'], '', 'etc-py3.json.gz')
        self.assertFileExists(data, 'etcfile')


class GenerateDiffExprTest(GeneratorBioProcessTestCase):

    @skipUnlessLargeFiles('DE-cuffdiff-py2.json.gz', 'DE-cuffdiff-py3.json.gz',
                          'DE-cuffdiff-py2.tab.gz', 'DE-cuffdiff-py3.tab.gz')
    def test_generate_diffexpr_cuffdiff(self):
        call_command('generate_diffexpr_cuffdiff', '-n=1', '-g=2', '--rseed')
        diffexpr = Data.objects.last()
        if diffexpr:
            self.assertEqual(diffexpr.process.type, 'data:differentialexpression:cuffdiff:')
            self.assertEqual(len(diffexpr.input['case']), 2)
            self.assertEqual(len(diffexpr.input['control']), 2)
            # NOTE: Python 2 and 3 produce different results even when setting random.seed() to the
            # same number due to https://docs.python.org/3/whatsnew/3.2.html#random
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
            # NOTE: Python 2 and 3 produce different results even when setting random.seed() to the
            # same number due to https://docs.python.org/3/whatsnew/3.2.html#random
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


class GenerateGeneSetTest(GeneratorBioProcessTestCase):

    def test_generate_geneset(self):
        call_command('generate_geneset', '-n=1', '--rseed')
        geneset = Data.objects.last()
        if geneset:
            self.assertEqual(geneset.process.type, 'data:geneset:')
            self.assertEqual(geneset.output['source'], 'UCSC')
            # NOTE: Python 2 and 3 produce different results even when setting random.seed() to the
            # same number due to https://docs.python.org/3/whatsnew/3.2.html#random
            if six.PY2:
                self.assertJSON(geneset, geneset.output['geneset_json'], '', 'geneset-py2.json.gz')
                self.assertFile(geneset, 'geneset', 'geneset-py2.tab.gz', compression='gzip')
            else:
                self.assertJSON(geneset, geneset.output['geneset_json'], '', 'geneset-py3.json.gz')
                self.assertFile(geneset, 'geneset', 'geneset-py3.tab.gz', compression='gzip')
        else:
            self.fail("Gene set not created")
