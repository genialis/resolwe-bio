# pylint: disable=missing-docstring
from __future__ import absolute_import, division, print_function, unicode_literals

from django.contrib.auth import get_user_model
from django.test import TestCase

from resolwe.flow.models import Data, DescriptorSchema, Process
from resolwe_bio.models import Sample


class SignalsSetTest(TestCase):
    def setUp(self):
        user_model = get_user_model()

        self.user = user_model.objects.create_user('test_user')
        DescriptorSchema.objects.create(name='Sample', slug='sample', contributor=self.user)
        self.process = Process.objects.create(name='Test process', contributor=self.user, flow_collection='sample')
        # `Sample`is created automatically when `Data` object is created
        self.data = Data.objects.create(name='Test data', contributor=self.user, process=self.process)

    def test_delete_last_data(self):
        self.data.delete()
        self.assertEqual(Sample.objects.count(), 0)

    def test_new_sample(self):
        data = Data.objects.create(name='Test data', contributor=self.user, process=self.process)
        sample = Sample.objects.last()
        self.assertTrue(sample.data.filter(pk=data.pk).exists())
