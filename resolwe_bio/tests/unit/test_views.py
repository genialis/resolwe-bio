# pylint: disable=missing-docstring
from __future__ import absolute_import, division, print_function, unicode_literals

from mock import MagicMock

from django.contrib.auth import get_user_model
from django.test import TestCase

from guardian.shortcuts import assign_perm, remove_perm
from rest_framework import exceptions
from rest_framework.test import APIRequestFactory

from resolwe.flow.models import Collection, Data, Process
from resolwe_bio.models import Sample
from resolwe_bio.views import SampleViewSet

factory = APIRequestFactory()  # pylint: disable=invalid-name


class SampleViewSetTest(TestCase):
    def setUp(self):
        user_model = get_user_model()

        self.user = user_model.objects.create_user('test_user')
        self.collection = Collection.objects.create(name="Test Collection", contributor=self.user)
        self.sample = Sample.objects.create(name="Test sample", contributor=self.user)
        process = Process.objects.create(name="Test process", contributor=self.user)
        self.data = Data.objects.create(name="Test data", contributor=self.user, process=process)
        self.data_2 = Data.objects.create(name="Test data 2", contributor=self.user, process=process)

        # another Data object to make sure that other objects are not processed
        Data.objects.create(name="Dummy data", contributor=self.user, process=process)

        self.sample.data.add(self.data)

        assign_perm('add_collection', self.user, self.collection)
        assign_perm('add_sample', self.user, self.sample)

        self.sampleviewset = SampleViewSet()

    def test_add_to_collection(self):
        request_mock = MagicMock(data={'ids': [self.collection.pk]}, user=self.user)
        self.sampleviewset.get_object = lambda: self.sample

        self.sampleviewset.add_to_collection(request_mock)

        self.assertEqual(self.collection.data.count(), 1)
        self.assertEqual(self.sample.collections.count(), 1)

    def test_remove_from_collection(self):
        # Manually add Sample and it's Data objects to the Collection
        self.sample.collections.add(self.collection.pk)
        self.collection.data.add(self.data)

        request_mock = MagicMock(data={'ids': [self.collection.pk]}, user=self.user)
        self.sampleviewset.get_object = lambda: self.sample

        self.sampleviewset.remove_from_collection(request_mock)

        self.assertEqual(self.collection.data.count(), 0)
        self.assertEqual(self.sample.collections.count(), 0)

    def test_add_remove_permissions(self):
        request_mock = MagicMock(data={'ids': [self.collection.pk]}, user=self.user)
        self.sampleviewset.get_object = lambda: self.sample

        remove_perm('add_collection', self.user, self.collection)

        with self.assertRaises(exceptions.PermissionDenied):
            self.sampleviewset.remove_from_collection(request_mock)

        with self.assertRaises(exceptions.PermissionDenied):
            self.sampleviewset.add_to_collection(request_mock)

    def test_add_data(self):
        self.sample.collections.add(self.collection)

        request_mock = MagicMock(data={'ids': [self.data.pk]}, user=self.user)
        self.sampleviewset.get_object = lambda: self.sample

        self.sampleviewset.add_data(request_mock)

        self.assertEqual(self.sample.data.count(), 1)
        self.assertEqual(self.collection.data.count(), 1)

    def test_remove_data(self):
        self.sample.data.add(self.data_2)
        self.sampleviewset.get_object = lambda: self.sample

        # sample is removed only when last data object is removed
        request_mock = MagicMock(data={'ids': [self.data.pk]}, user=self.user)
        self.sampleviewset.remove_data(request_mock)
        self.assertEqual(Sample.objects.count(), 1)
        request_mock = MagicMock(data={'ids': [self.data_2.pk]}, user=self.user)
        self.sampleviewset.remove_data(request_mock)
        self.assertEqual(Sample.objects.count(), 0)
