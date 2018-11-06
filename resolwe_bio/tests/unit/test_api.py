# pylint: disable=missing-docstring
from rest_framework import status
from rest_framework.test import APIRequestFactory, force_authenticate

from resolwe.flow.models import Data, DescriptorSchema, Process
from resolwe.flow.views import DataViewSet
from resolwe.test import TestCase

factory = APIRequestFactory()  # pylint: disable=invalid-name


class TestDataViewSetFilters(TestCase):
    def setUp(self):
        super().setUp()

        self.data_viewset = DataViewSet.as_view(actions={
            'get': 'list',
        })

        self.proc = Process.objects.create(
            type='data:test:process',
            slug='test-process',
            version='1.0.0',
            contributor=self.contributor,
            entity_type='test-schema',
            entity_descriptor_schema='test-schema',
            input_schema=[{'name': 'input_data', 'type': 'data:test:', 'required': False}],
            output_schema=[
                {'name': 'source', 'type': 'basic:string:'},
                {'name': 'species', 'type': 'basic:string:'},
                {'name': 'build', 'type': 'basic:string:'},
                {'name': 'feature_type', 'type': 'basic:string:'},
            ]
        )

        self.descriptor_schema = DescriptorSchema.objects.create(
            slug='test-schema',
            version='1.0.0',
            contributor=self.contributor,
        )

        self.data = []
        for index in range(10):
            data = Data.objects.create(
                name='Data {}'.format(index),
                contributor=self.contributor,
                process=self.proc,
                status=Data.STATUS_DONE,
                output={
                    'source': 'NCBI{}'.format(index),
                    'species': 'Mus musculus' if index < 5 else 'Homo sapiens',
                    'build': 'B{}'.format(index),
                    'feature_type': 'gene' if index < 5 else 'foo',
                }
            )

            self.data.append(data)

    def _check_filter(self, query_args, expected):
        request = factory.get('/', query_args, format='json')
        force_authenticate(request, self.admin)
        response = self.data_viewset(request)
        expected = [item.pk for item in expected]

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), len(expected))
        for item in response.data:
            self.assertIn(item['id'], expected)

    def test_filter_source(self):
        self._check_filter({'source': 'NCBI1'}, [self.data[1]])
        self._check_filter({'source': 'NCBI5'}, [self.data[5]])
        self._check_filter({'source': 'UCSC'}, [])

    def test_filter_species(self):
        self._check_filter({'species': 'Mus musculus'}, self.data[:5])
        self._check_filter({'species': 'Mus'}, self.data[:5])
        self._check_filter({'species': 'musculus'}, self.data[:5])
        self._check_filter({'species': 'Homo sapiens'}, self.data[5:])
        self._check_filter({'species': 'Homo'}, self.data[5:])
        self._check_filter({'species': 'sapiens'}, self.data[5:])

    def test_filter_build(self):
        self._check_filter({'build': 'B1'}, [self.data[1]])
        self._check_filter({'build': 'B5'}, [self.data[5]])
        self._check_filter({'build': 'XXX'}, [])

    def test_filter_feature_type(self):
        self._check_filter({'feature_type': 'gene'}, self.data[:5])
        self._check_filter({'feature_type': 'foo'}, self.data[5:])

    def test_filter_text(self):
        self._check_filter({'text': 'data'}, self.data)
        self._check_filter({'text': 'mus'}, self.data[:5])
        self._check_filter({'text': 'sapiens'}, self.data[5:])
        self._check_filter({'text': 'contributor'}, self.data)
        self._check_filter({'text': 'blablabla'}, [])
