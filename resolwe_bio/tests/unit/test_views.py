from __future__ import absolute_import, division, print_function, unicode_literals

from mock import MagicMock

from django.contrib.auth import get_user_model
from django.test import TestCase

from guardian import shortcuts
from rest_framework import exceptions

from resolwe.flow.models import Collection, Data, Process
from resolwe_bio.models import Sample
from resolwe_bio.views import SampleViewSet


class SampleViewSetTest(TestCase):
    def setUp(self):
        User = get_user_model()

        self.user = User.objects.create_user('test_user')
        self.collection = Collection.objects.create(name="Test Collection", contributor=self.user)
        self.sample = Sample.objects.create(name="Test sample", contributor=self.user)
        process = Process.objects.create(name="Test process", contributor=self.user)
        self.data = Data.objects.create(name="Test data", contributor=self.user, process=process)

        # another Data object to make sure that other objects are not processed
        Data.objects.create(name="Dummy data", contributor=self.user, process=process)

        self.sample.data.add(self.data)

        shortcuts.assign_perm('add_collection', self.user, self.collection)

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

        shortcuts.remove_perm('add_collection', self.user, self.collection)

        with self.assertRaises(exceptions.PermissionDenied):
            self.sampleviewset.remove_from_collection(request_mock)

        with self.assertRaises(exceptions.PermissionDenied):
            self.sampleviewset.add_to_collection(request_mock)
