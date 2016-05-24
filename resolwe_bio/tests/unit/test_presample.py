from __future__ import absolute_import, division, print_function, unicode_literals

from django.contrib.auth import get_user_model
from django.core.urlresolvers import reverse

from rest_framework.test import APITestCase, APIRequestFactory, force_authenticate
from rest_framework import status

from resolwe_bio.models import Sample
from resolwe_bio.views import SampleViewSet, PresampleViewSet


class PresampleTestCase(APITestCase):
    def setUp(self):
        self.user = get_user_model().objects.create(username='user', is_superuser=True)
        self.sample = Sample.objects.create(contributor=self.user, name="Test sample")

        sample_viewset = SampleViewSet()
        self.sample_queryset = sample_viewset.get_queryset()
        presample_viewset = PresampleViewSet()
        self.presample_queryset = presample_viewset.get_queryset()

        url_mapping = {
            'get': 'retrieve',
            'put': 'update',
            'patch': 'partial_update',
            'delete': 'destroy',
        }

        self.sample_view = SampleViewSet.as_view(url_mapping)
        self.presample_view = PresampleViewSet.as_view(url_mapping)

        self.factory = APIRequestFactory()

    @staticmethod
    def get_detail_url(endpoint, pk):
        return reverse('resolwebio-api:{}-detail'.format(endpoint), kwargs={'pk': pk})

    def test_querysets(self):
        self.assertEqual(self.sample_queryset.count(), 0)
        self.assertEqual(self.presample_queryset.count(), 1)

        self.sample.presample = False  # mark sample as annotated
        self.sample.save()
        self.assertEqual(self.sample_queryset.count(), 1)
        self.assertEqual(self.presample_queryset.count(), 0)

    def test_upgrade_presample(self):
        """`Presample` can be transformed into `Sample`"""
        url = self.get_detail_url('presample', self.sample.pk)
        request = self.factory.patch(url, {'presample': False}, format='json')
        force_authenticate(request, user=self.user)
        response = self.presample_view(request, pk=self.sample.pk)

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.sample.refresh_from_db()
        self.assertFalse(self.sample.presample)

    def test_revert_sample(self):
        """`Sample` cannot be reverted back in to `Presample`"""
        self.sample.presample = False
        self.sample.save()

        url = self.get_detail_url('sample', self.sample.pk)
        request = self.factory.patch(url, {'presample': False}, format='json')
        force_authenticate(request, user=self.user)
        response = self.sample_view(request, pk=self.sample.pk)

        # `response.content` is "A server error occurred.", but this is OK,
        # because request is treated as request with no data (after `presample`
        # is removed). This is the reason for status code 204.
        self.assertEqual(response.status_code, status.HTTP_204_NO_CONTENT)
        self.sample.refresh_from_db()
        self.assertFalse(self.sample.presample)

    def test_wrong_endpoint(self):
        """`Sample` cannot be changed through `Presample` endpoint"""
        self.sample.presample = False
        self.sample.save()

        url = self.get_detail_url('presample', self.sample.pk)
        request = self.factory.patch(url, {'presample': False}, format='json')
        force_authenticate(request, user=self.user)
        response = self.presample_view(request, pk=self.sample.pk)

        self.assertEqual(response.status_code, status.HTTP_404_NOT_FOUND)
        self.sample.refresh_from_db()
        self.assertFalse(self.sample.presample)
