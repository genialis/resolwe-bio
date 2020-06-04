from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase

from resolwe.test import TestCase

from ..models import Feature


class FeatureTestCase(TestCase, APITestCase):
    @staticmethod
    def create_feature(index, source):
        return Feature.objects.create(
            source=source,
            feature_id="FT-{}".format(index),
            species="Lorem ipsum",
            type=Feature.TYPE_GENE,
            sub_type=Feature.SUBTYPE_PROTEIN_CODING,
            name="FOO{}".format(index),
            full_name="Foobarius machinus",
            aliases=["BAR{}".format(index), "BTMK{}".format(index), "SHARED"],
        )

    @classmethod
    def setUpTestData(cls):
        cls.features = []

        for i in range(7):
            cls.features.append(FeatureTestCase.create_feature(i, "NCBI"))

        for i in range(7, 10):
            cls.features.append(FeatureTestCase.create_feature(i, "XSRC"))

    def assertFeatureEqual(self, data, feature):
        self.assertEqual(data["source"], feature.source)
        self.assertEqual(data["species"], feature.species)
        self.assertEqual(data["type"], feature.type)
        self.assertEqual(data["sub_type"], feature.sub_type)
        self.assertEqual(data["name"], feature.name)
        self.assertEqual(data["full_name"], feature.full_name)
        self.assertEqual(data["aliases"], feature.aliases)

    def test_feature_search(self):
        FEATURE_SEARCH_URL = reverse("resolwebio-api:kb_feature_search")

        # Test without any query.
        response = self.client.get(FEATURE_SEARCH_URL, format="json")
        self.assertEqual(len(response.data), len(self.features))

        # Test with empty query.
        response = self.client.get(FEATURE_SEARCH_URL, {"query": ""}, format="json")
        self.assertEqual(len(response.data), 10)

        # Test with non-matching query.
        response = self.client.get(FEATURE_SEARCH_URL, {"query": "F1"}, format="json")
        self.assertEqual(len(response.data), 0)

        for feature in self.features:
            # Test query by feature name.
            response = self.client.get(
                FEATURE_SEARCH_URL, {"query": feature.name}, format="json"
            )
            self.assertEqual(len(response.data), 1)
            self.assertFeatureEqual(response.data[0], feature)

            # Test query by alias.
            for alias in feature.aliases:
                if alias == "SHARED":
                    continue

                response = self.client.get(
                    FEATURE_SEARCH_URL, {"query": alias}, format="json"
                )
                self.assertEqual(len(response.data), 1)
                self.assertFeatureEqual(response.data[0], feature)

        # Test query by shared alias.
        response = self.client.get(
            FEATURE_SEARCH_URL, {"query": "SHARED"}, format="json"
        )
        self.assertEqual(len(response.data), len(self.features))

        # Test query by multiple gene IDs.
        response = self.client.get(
            FEATURE_SEARCH_URL, {"query": "FOO1,FOO2,FT-7"}, format="json"
        )
        self.assertEqual(len(response.data), 3)

        response = self.client.get(
            FEATURE_SEARCH_URL, {"query": "FOO1,SHARED,FT-7"}, format="json"
        )
        self.assertEqual(len(response.data), len(self.features))

        # Test query by source.
        response = self.client.get(
            FEATURE_SEARCH_URL, {"query": "FOO1", "source": "NCBI"}, format="json"
        )
        self.assertEqual(len(response.data), 1)

        response = self.client.get(
            FEATURE_SEARCH_URL, {"query": "FOO1", "source": "FOO"}, format="json"
        )
        self.assertEqual(len(response.data), 0)

        response = self.client.get(
            FEATURE_SEARCH_URL,
            {"query": "FOO1,FOO2,FT-3", "source": "NCBI"},
            format="json",
        )
        self.assertEqual(len(response.data), 3)

        response = self.client.get(
            FEATURE_SEARCH_URL,
            {"query": "FOO7,FOO8,FT-9", "source": "XSRC"},
            format="json",
        )
        self.assertEqual(len(response.data), 3)

        response = self.client.get(
            FEATURE_SEARCH_URL,
            {"query": "FOO1,FOO2,FT-7", "source": "FOO"},
            format="json",
        )
        self.assertEqual(len(response.data), 0)

        # Test query with a lot of features.
        response = self.client.get(
            FEATURE_SEARCH_URL,
            {"query": ",".join([str(x) for x in range(1024)])},
            format="json",
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        # Test search by feature_id.
        response = self.client.get(
            FEATURE_SEARCH_URL, {"feature_id": "FOO", "source": "XSRC"}, format="json"
        )
        self.assertEqual(len(response.data), 0)

        response = self.client.get(
            FEATURE_SEARCH_URL, {"feature_id": "FOO1", "source": "XSRC"}, format="json"
        )
        self.assertEqual(len(response.data), 0)

        response = self.client.get(
            FEATURE_SEARCH_URL, {"feature_id": "FT-9", "source": "XSRC"}, format="json"
        )
        self.assertEqual(len(response.data), 1)

        response = self.client.get(
            FEATURE_SEARCH_URL,
            {"feature_id": "FT-1,FT-2", "source": "NCBI"},
            format="json",
        )
        self.assertEqual(len(response.data), 2)

    def test_feature_autocomplete(self):
        FEATURE_AUTOCOMPLETE_URL = reverse("resolwebio-api:kb_feature_autocomplete")

        # Test empty query.
        response = self.client.get(
            FEATURE_AUTOCOMPLETE_URL, {"query": ""}, format="json"
        )
        self.assertEqual(len(response.data), 0)

        response = self.client.get(
            FEATURE_AUTOCOMPLETE_URL, {"query": "", "source": "NCBI"}, format="json"
        )
        self.assertEqual(len(response.data), 0)

        response = self.client.get(FEATURE_AUTOCOMPLETE_URL, format="json")
        self.assertEqual(len(response.data), 0)

        # Test non-matching query.
        response = self.client.get(
            FEATURE_AUTOCOMPLETE_URL, {"query": "FOU"}, format="json"
        )
        self.assertEqual(len(response.data), 0)

        # Test partial name query.
        response = self.client.get(
            FEATURE_AUTOCOMPLETE_URL, {"query": "FO"}, format="json"
        )
        self.assertEqual(len(response.data), len(self.features))

        # Test partial name query with source.
        response = self.client.get(
            FEATURE_AUTOCOMPLETE_URL, {"query": "FO", "source": "NCBI"}, format="json"
        )
        self.assertEqual(len(response.data), 7)

        # Test partial alias query.
        response = self.client.get(
            FEATURE_AUTOCOMPLETE_URL, {"query": "SHAR"}, format="json"
        )
        self.assertEqual(len(response.data), len(self.features))

        # Test partial alias query with source.
        response = self.client.get(
            FEATURE_AUTOCOMPLETE_URL, {"query": "SHAR", "source": "XSRC"}, format="json"
        )
        self.assertEqual(len(response.data), 3)

        for feature in self.features:
            # Test query by full feature name.
            response = self.client.get(
                FEATURE_AUTOCOMPLETE_URL,
                {"query": feature.name, "source": feature.source},
                format="json",
            )
            self.assertEqual(len(response.data), 1)
            self.assertFeatureEqual(response.data[0], feature)

            # Test query by alias.
            for alias in feature.aliases:
                if alias == "SHARED":
                    continue

                response = self.client.get(
                    FEATURE_AUTOCOMPLETE_URL, {"query": alias}, format="json"
                )
                self.assertEqual(len(response.data), 1)
                self.assertFeatureEqual(response.data[0], feature)

