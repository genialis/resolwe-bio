# pylint: disable=missing-docstring,invalid-name,no-member
from __future__ import absolute_import, division, print_function, unicode_literals

from django.core.urlresolvers import reverse

from rest_framework import status
from rest_framework.test import APITestCase

from ..models import Mapping


class MappingTestCase(APITestCase):
    def setUp(self):
        self.mappings = []
        for i in range(10):
            self.mappings.append(Mapping.objects.create(
                relation_type='crossdb',
                source_db='A',
                source_id='FT-{}'.format(i),
                target_db='B',
                target_id='ANOTHER-{}'.format(i),
            ))

    def assertMappingEqual(self, data, mapping):
        self.assertEqual(data['relation_type'], mapping.relation_type)
        self.assertEqual(data['source_db'], mapping.source_db)
        self.assertEqual(data['source_id'], mapping.source_id)
        self.assertEqual(data['target_db'], mapping.target_db)
        self.assertEqual(data['target_id'], mapping.target_id)

    def test_lookup(self):
        MAPPING_URL = reverse('resolwebio-api:mapping-list')

        # Test without any query.
        response = self.client.get(MAPPING_URL, format='json')
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), len(self.mappings))

        # Test lookup by source_db, target_db and a single source feature identifier.
        response = self.client.get(MAPPING_URL, {
            'source_db': 'A',
            'target_db': 'B',
            'source_id': 'FT-0',
        }, format='json')
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), 1)
        self.assertMappingEqual(response.data[0], self.mappings[0])

        # Test lookup by source_db, target_db and a list of source feature identifiers.
        response = self.client.post(MAPPING_URL, {
            'source_db': 'A',
            'target_db': 'B',
            'source_id__in': ['FT-0', 'FT-1', 'FT-5']
        }, format='json')
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), 3)
        self.assertMappingEqual(response.data[0], self.mappings[0])
        self.assertMappingEqual(response.data[1], self.mappings[1])
        self.assertMappingEqual(response.data[2], self.mappings[5])
