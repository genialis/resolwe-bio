from rest_framework import status
from rest_framework.test import APIRequestFactory, force_authenticate

from resolwe.flow.models import Data, DescriptorSchema, Entity, Process
from resolwe.flow.views import DataViewSet, EntityViewSet
from resolwe.test import TestCase

factory = APIRequestFactory()


class BaseViewSetFiltersTest(TestCase):
    def _check_filter(
        self, query_args, expected, expected_status_code=status.HTTP_200_OK
    ):
        """Check that query_args filter to expected queryset."""
        request = factory.get("/", query_args, format="json")
        force_authenticate(request, self.admin)
        response = self.viewset(request)

        if status.is_success(response.status_code):
            self.assertEqual(len(response.data), len(expected))
            self.assertCountEqual(
                [item.pk for item in expected],
                [item["id"] for item in response.data],
            )
        else:
            self.assertEqual(response.status_code, expected_status_code)
            response.render()
            return response


class TestDataViewSetFilters(BaseViewSetFiltersTest):
    def setUp(self):
        super().setUp()

        self.viewset = DataViewSet.as_view(
            actions={
                "get": "list",
            }
        )

        self.proc = Process.objects.create(
            type="data:test:process",
            slug="test-process",
            version="1.0.0",
            contributor=self.contributor,
            entity_type="test-schema",
            entity_descriptor_schema="test-schema",
            input_schema=[
                {"name": "input_data", "type": "data:test:", "required": False}
            ],
            output_schema=[
                {"name": "source", "type": "basic:string:"},
                {"name": "species", "type": "basic:string:"},
                {"name": "build", "type": "basic:string:"},
                {"name": "feature_type", "type": "basic:string:"},
            ],
        )

        self.descriptor_schema = DescriptorSchema.objects.create(
            slug="test-schema",
            version="1.0.0",
            contributor=self.contributor,
        )

        self.data = []
        for index in range(10):
            data = Data.objects.create(
                name="Data {}".format(index),
                contributor=self.contributor,
                process=self.proc,
                status=Data.STATUS_DONE,
                output={
                    "source": "NCBI{}".format(index),
                    "species": "Mus musculus" if index < 5 else "Homo sapiens",
                    "build": "B{}".format(index),
                    "feature_type": "gene" if index < 5 else "foo",
                },
            )

            self.data.append(data)

    def test_filter_source(self):
        self._check_filter({"source": "NCBI1"}, [self.data[1]])
        self._check_filter({"source": "NCBI5"}, [self.data[5]])
        self._check_filter({"source": "UCSC"}, [])

    def test_filter_species(self):
        self._check_filter({"species": "Mus musculus"}, self.data[:5])
        self._check_filter({"species": "Mus"}, self.data[:5])
        self._check_filter({"species": "musculus"}, self.data[:5])
        self._check_filter({"species": "Homo sapiens"}, self.data[5:])
        self._check_filter({"species": "Homo"}, self.data[5:])
        self._check_filter({"species": "sapiens"}, self.data[5:])

    def test_filter_build(self):
        self._check_filter({"build": "B1"}, [self.data[1]])
        self._check_filter({"build": "B5"}, [self.data[5]])
        self._check_filter({"build": "XXX"}, [])

    def test_filter_feature_type(self):
        self._check_filter({"feature_type": "gene"}, self.data[:5])
        self._check_filter({"feature_type": "foo"}, self.data[5:])

    def test_filter_text(self):
        self._check_filter({"text": "data"}, self.data)
        self._check_filter({"text": "mus"}, self.data[:5])
        self._check_filter({"text": "sapiens"}, self.data[5:])
        self._check_filter({"text": "contributor"}, self.data)
        self._check_filter({"text": "blablabla"}, [])

        # Check that changes are applied immediately.
        self.data[0].output["species"] = "Rat rattus"
        self.data[0].save()

        self._check_filter({"text": "rat rattus"}, [self.data[0]])


class TestEntityViewSetFilters(BaseViewSetFiltersTest):
    def setUp(self):
        super().setUp()

        self.viewset = EntityViewSet.as_view(
            actions={
                "get": "list",
            }
        )

        self.descriptor_schema = DescriptorSchema.objects.create(
            slug="sample",
            contributor=self.contributor,
            schema=[
                {
                    "label": "General",
                    "name": "general",
                    "group": [
                        {
                            "name": "species",
                            "type": "basic:string:",
                            "label": "Species",
                        },
                    ],
                }
            ],
        )

        self.entities = [
            Entity.objects.create(
                name="Test entity 1",
                contributor=self.contributor,
                descriptor_schema=self.descriptor_schema,
                descriptor={"general": {"species": "Homo Sapiens"}},
            ),
            Entity.objects.create(
                name="Test entity 2",
                contributor=self.contributor,
                descriptor_schema=self.descriptor_schema,
                descriptor={"general": {"species": "Mus musculus"}},
            ),
        ]

    def test_filter_species(self):
        self._check_filter({"species": "Mus musculus"}, [self.entities[1]])
        self._check_filter({"species": "Mus"}, [self.entities[1]])
        self._check_filter({"species": "musculus"}, [self.entities[1]])
        self._check_filter({"species": "Homo sapiens"}, [self.entities[0]])
        self._check_filter({"species": "Homo"}, [self.entities[0]])
        self._check_filter({"species": "sapiens"}, [self.entities[0]])

    def test_filter_text(self):
        # By species.
        self._check_filter({"text": "Mus musculus"}, [self.entities[1]])
        self._check_filter({"text": "Mus"}, [self.entities[1]])
        self._check_filter({"text": "musculus"}, [self.entities[1]])
        self._check_filter({"text": "Homo sapiens"}, [self.entities[0]])
        self._check_filter({"text": "Homo"}, [self.entities[0]])
        self._check_filter({"text": "sapiens"}, [self.entities[0]])

        # Check that changes are applied immediately.
        self.entities[0].descriptor["general"]["species"] = "Rat rattus"
        self.entities[0].save()

        self._check_filter({"text": "rat rattus"}, [self.entities[0]])
