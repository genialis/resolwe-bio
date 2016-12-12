# pylint: disable=missing-docstring,invalid-name,no-member
from __future__ import absolute_import, division, print_function, unicode_literals

from django.contrib.auth.models import User
from django.core.urlresolvers import reverse
from django.core.management import call_command

from rest_framework import status
from rest_framework.test import APITestCase

from ..models import Feature


class FeatureTestCase(APITestCase):
    def setUp(self):
        self.features = []
        for i in range(10):
            self.features.append(Feature.objects.create(
                source='NCBI',
                feature_id='FT-{}'.format(i),
                species='Lorem ipsum',
                type=Feature.TYPE_GENE,
                sub_type=Feature.SUBTYPE_PROTEIN_CODING,
                name='FOO{}'.format(i),
                full_name='Foobarius machinus',
                aliases=['BAR{}'.format(i), 'BTMK{}'.format(i), 'SHARED']
            ))

        call_command('rebuild_index', interactive=False, verbosity=0)

    def assertFeatureEqual(self, data, feature):
        self.assertEqual(data['source'], feature.source)
        self.assertEqual(data['species'], feature.species)
        self.assertEqual(data['type'], feature.type)
        self.assertEqual(data['sub_type'], feature.sub_type)
        self.assertEqual(data['name'], feature.name)
        self.assertEqual(data['full_name'], feature.full_name)
        self.assertEqual(data['aliases'], feature.aliases)

    def test_feature_search(self):
        FEATURE_SEARCH_URL = reverse('resolwebio-api:kb_feature_search-list')

        # Test without any query.
        response = self.client.get(FEATURE_SEARCH_URL, format='json')
        self.assertEqual(len(response.data), len(self.features))

        # Test with non-matching query.
        response = self.client.get(FEATURE_SEARCH_URL, {'query': 'FO'}, format='json')
        self.assertEqual(len(response.data), 0)

        for feature in self.features:
            # Test query by feature name.
            response = self.client.get(FEATURE_SEARCH_URL, {'query': feature.name}, format='json')
            self.assertEqual(len(response.data), 1)
            self.assertFeatureEqual(response.data[0], feature)

            # Test query by alias.
            for alias in feature.aliases:
                if alias == 'SHARED':
                    continue

                response = self.client.get(FEATURE_SEARCH_URL, {'query': alias}, format='json')
                self.assertEqual(len(response.data), 1)
                self.assertFeatureEqual(response.data[0], feature)

        # Test query by shared alias.
        response = self.client.get(FEATURE_SEARCH_URL, {'query': 'SHARED'}, format='json')
        self.assertEqual(len(response.data), len(self.features))

        # Test query by multiple gene IDs.
        response = self.client.post(FEATURE_SEARCH_URL, {'query': ['FOO1', 'FOO2', 'FT-7']}, format='json')
        self.assertEqual(len(response.data), 3)

        response = self.client.post(FEATURE_SEARCH_URL, {'query': ['FOO1', 'SHARED', 'FT-7']}, format='json')
        self.assertEqual(len(response.data), len(self.features))

        # Test query by source.
        response = self.client.get(FEATURE_SEARCH_URL, {'query': 'FOO1', 'source': 'NCBI'}, format='json')
        self.assertEqual(len(response.data), 1)

        response = self.client.get(FEATURE_SEARCH_URL, {'query': 'FOO1', 'source': 'FOO'}, format='json')
        self.assertEqual(len(response.data), 0)

        response = self.client.post(FEATURE_SEARCH_URL, {'query': ['FOO1', 'FOO2', 'FT-7'], 'source': 'NCBI'},
                                    format='json')
        self.assertEqual(len(response.data), 3)

        response = self.client.post(FEATURE_SEARCH_URL, {'query': ['FOO1', 'FOO2', 'FT-7'], 'source': 'FOO'},
                                    format='json')
        self.assertEqual(len(response.data), 0)

    def test_feature_autocomplete(self):
        FEATURE_AUTOCOMPLETE_URL = reverse('resolwebio-api:kb_feature_autocomplete-list')

        # Test non-matching query.
        response = self.client.get(FEATURE_AUTOCOMPLETE_URL, {'query': 'FOU'}, format='json')
        self.assertEqual(len(response.data), 0)

        # Test partial name query.
        response = self.client.get(FEATURE_AUTOCOMPLETE_URL, {'query': 'FO'}, format='json')
        self.assertEqual(len(response.data), len(self.features))

        # Test partial alias query.
        response = self.client.get(FEATURE_AUTOCOMPLETE_URL, {'query': 'SHAR'}, format='json')
        self.assertEqual(len(response.data), len(self.features))

        for feature in self.features:
            # Test query by full feature name.
            response = self.client.get(FEATURE_AUTOCOMPLETE_URL, {'query': feature.name}, format='json')
            self.assertEqual(len(response.data), 1)
            self.assertFeatureEqual(response.data[0], feature)

            # Test query by alias.
            for alias in feature.aliases:
                if alias == 'SHARED':
                    continue

                response = self.client.get(FEATURE_AUTOCOMPLETE_URL, {'query': alias}, format='json')
                self.assertEqual(len(response.data), 1)
                self.assertFeatureEqual(response.data[0], feature)

    def test_feature_admin(self):
        # Test that only an admin can access the endpoint.
        response = self.client.get(reverse('resolwebio-api:feature-list'), format='json')
        self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)

        # Authenticate as normal user.
        normal_user = User.objects.create_user('tester', 'tester@genialis.com', 'tester')
        self.client.force_authenticate(user=normal_user)
        response = self.client.get(reverse('resolwebio-api:feature-list'), format='json')
        self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)

        # Authenticate as admin.
        admin_user = User.objects.create_superuser('admin', 'admin@genialis.com', 'admin')
        self.client.force_authenticate(user=admin_user)

        # Test listing and detailed access features.
        response = self.client.get(reverse('resolwebio-api:feature-list'), format='json')
        self.assertEqual(len(response.data), len(self.features))
        for data, feature in zip(sorted(response.data, key=lambda x: x['id']), self.features):
            self.assertFeatureEqual(data, feature)

            detail = self.client.get(
                reverse('resolwebio-api:feature-detail', kwargs={'pk': data['id']}),
                format='json'
            )
            self.assertFeatureEqual(detail.data, feature)

        # Test adding new features.
        response = self.client.post(
            reverse('resolwebio-api:feature-list'),
            {
                'source': 'NCBI',
                'feature_id': 'TEST',
                'species': 'Lorem ipsum',
                'type': Feature.TYPE_GENE,
                'sub_type': Feature.SUBTYPE_PROTEIN_CODING,
                'name': 'TEST',
                'full_name': 'Test test',
                'aliases': ['ANOTHER', 'ALIAS'],
            },
        )
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(Feature.objects.get(source='NCBI', feature_id='TEST').full_name, 'Test test')

        # Test duplicate insert (should automatically update).
        response = self.client.post(
            reverse('resolwebio-api:feature-list'),
            {
                'source': 'NCBI',
                'feature_id': 'TEST',
                'species': 'Lorem ipsum',
                'type': Feature.TYPE_GENE,
                'sub_type': Feature.SUBTYPE_PROTEIN_CODING,
                'name': 'TEST',
                'full_name': 'Modified test',
                'aliases': ['ANOTHER', 'ALIAS'],
            },
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(Feature.objects.get(source='NCBI', feature_id='TEST').full_name, 'Modified test')

        # Test missing feature_id insert.
        response = self.client.post(
            reverse('resolwebio-api:feature-list'),
            {
                'source': 'NCBI',
                'species': 'Lorem ipsum',
                'type': Feature.TYPE_GENE,
                'sub_type': Feature.SUBTYPE_PROTEIN_CODING,
                'name': 'TEST',
                'full_name': 'Test test',
                'aliases': ['ANOTHER', 'ALIAS'],
            },
        )
        self.assertEqual(response.status_code, status.HTTP_400_BAD_REQUEST)
