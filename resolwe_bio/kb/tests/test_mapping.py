from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from resolwe.test import TestCase

from ..models import Mapping


class MappingTestCase(APITestCase, TestCase):
    def setUp(self):
        super(MappingTestCase, self).setUp()

        self.mappings = []
        for i in range(10):
            self.mappings.append(
                Mapping.objects.create(
                    relation_type="crossdb",
                    source_db="SRC",
                    source_id="FT{}".format(i),
                    source_species="Mus musculus",
                    target_db="TGT",
                    target_id="ANOTHER{}".format(i),
                    target_species="Mus musculus",
                )
            )

    def assertMappingEqual(self, data, mapping):
        self.assertEqual(data["relation_type"], mapping.relation_type)
        self.assertEqual(data["source_db"], mapping.source_db)
        self.assertEqual(data["source_id"], mapping.source_id)
        self.assertEqual(data["source_species"], mapping.source_species)
        self.assertEqual(data["target_db"], mapping.target_db)
        self.assertEqual(data["target_id"], mapping.target_id)
        self.assertEqual(data["target_species"], mapping.target_species)

    def test_lookup(self):
        MAPPING_URL = reverse("resolwebio-api:kb_mapping_search")

        # Test without any query.
        response = self.client.get(MAPPING_URL, format="json")
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), len(self.mappings))

        # Test lookup by source_db, target_db and a single source feature identifier.
        response = self.client.get(
            MAPPING_URL,
            {"source_db": "SRC", "target_db": "TGT", "source_id": "FT0",},
            format="json",
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), 1)
        self.assertMappingEqual(response.data[0], self.mappings[0])

        # Test lookup by source_db, target_db and a list of source feature identifiers.
        response = self.client.get(
            MAPPING_URL,
            {"source_db": "SRC", "target_db": "TGT", "source_id__in": "FT0,FT1,FT5"},
            format="json",
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), 3)
        self.assertMappingEqual(response.data[0], self.mappings[0])
        self.assertMappingEqual(response.data[1], self.mappings[1])
        self.assertMappingEqual(response.data[2], self.mappings[5])

        # Test query with a lot of source ids.
        response = self.client.get(
            MAPPING_URL,
            {
                "source_db": "SRC",
                "target_db": "TGT",
                "source_id__in": ",".join([str(x) for x in range(512)]),
            },
            format="json",
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
