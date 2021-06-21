from rest_framework import status
from rest_framework.test import APIRequestFactory, force_authenticate

from resolwe.flow.models import Data, DescriptorSchema, Entity, Process
from resolwe.flow.views import DataViewSet, EntityViewSet
from resolwe.test import ProcessTestCase, TestCase

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


# Process test case is required to register descriptor schemas which are used in the tests. Note
# that processes cannot be triggered as tests runs in a transaction.
class TestEntityViewSetFilters(ProcessTestCase, BaseViewSetFiltersTest):
    def setUp(self):
        super().setUp()

        self.viewset = EntityViewSet.as_view(
            actions={
                "get": "list",
            }
        )

        clinical_schema = DescriptorSchema.objects.get(slug="general-clinical")
        sample_schema = DescriptorSchema.objects.get(slug="sample")

        self.entities = [
            Entity.objects.create(
                name="Test entity 1",
                contributor=self.contributor,
                descriptor_schema=clinical_schema,
                descriptor={
                    "general": {"species": "Homo sapiens"},
                    "disease_information": {
                        "disease_type": "Colorectal cancer",
                        "disease_status": "Progresive",
                    },
                    "subject_information": {
                        "batch": 1,
                        "group": "Pre",
                        "subject_id": "P-006",
                        "sample_label": "CRC",
                    },
                    "immuno_oncology_treatment_type": {
                        "io_drug": "D-00A",
                        "io_treatment": "single",
                    },
                    "response_and_survival_analysis": {
                        "pfs_event": "no",
                        "confirmed_bor": "pd",
                    },
                },
            ),
            Entity.objects.create(
                name="Test entity 2",
                contributor=self.contributor,
                descriptor_schema=clinical_schema,
                descriptor={
                    "general": {"species": "Homo sapiens"},
                    "disease_information": {
                        "disease_type": "Mesothelioma",
                        "disease_status": "Regresive",
                    },
                    "subject_information": {
                        "batch": 2,
                        "group": "Post",
                        "subject_id": "P-019",
                        "sample_label": "Meso",
                    },
                    "immuno_oncology_treatment_type": {
                        "io_drug": "D-12A",
                        "io_treatment": "combo",
                    },
                    "response_and_survival_analysis": {
                        "pfs_event": "yes",
                        "confirmed_bor": "sd",
                    },
                },
            ),
            Entity.objects.create(
                name="Test entity 3",
                contributor=self.contributor,
                descriptor_schema=sample_schema,
                descriptor={
                    "general": {
                        "species": "Mus musculus",
                        "description": "First sample",
                        "biosample_source": "lung",
                        "biosample_treatment": "dmso",
                    }
                },
            ),
            Entity.objects.create(
                name="Test entity 4",
                contributor=self.contributor,
                descriptor_schema=sample_schema,
                descriptor={
                    "general": {
                        "species": "Mus musculus",
                        "description": "Second sample",
                        "biosample_source": "liver",
                        "biosample_treatment": "koh",
                    }
                },
            ),
        ]

    def test_filter_species(self):
        self._check_filter({"species": "Mus musculus"}, self.entities[2:])
        self._check_filter({"species": "Mus"}, self.entities[2:])
        self._check_filter({"species": "musculus"}, self.entities[2:])
        self._check_filter({"species": "Homo sapiens"}, self.entities[:2])
        self._check_filter({"species": "Homo"}, self.entities[:2])
        self._check_filter({"species": "sapiens"}, self.entities[:2])

    def test_filter_text(self):
        # By species.
        self._check_filter({"text": "Mus musculus"}, self.entities[2:])
        self._check_filter({"text": "Mus"}, self.entities[2:])
        self._check_filter({"text": "musculus"}, self.entities[2:])
        self._check_filter({"text": "Homo sapiens"}, self.entities[:2])
        self._check_filter({"text": "Homo"}, self.entities[:2])
        self._check_filter({"text": "sapiens"}, self.entities[:2])

        # Check that changes are applied immediately.
        self.entities[0].descriptor["general"]["species"] = "Rattus norvegicus"
        self.entities[0].save()

        self._check_filter({"text": "rattus norvegicus"}, [self.entities[0]])

    def test_filter_descriptor(self):
        # Descriptor / Subject Information / Sample Label
        self._check_filter(
            {"descriptor__subject_information__sample_label__icontains": "CRC"},
            [self.entities[0]],
        )
        self._check_filter(
            {"descriptor__subject_information__sample_label__icontains": "CR"},
            [self.entities[0]],
        )

        # Descriptor / Subject Information / Subject ID
        self._check_filter(
            {"descriptor__subject_information__subject_id__icontains": "P-006"},
            [self.entities[0]],
        )
        self._check_filter(
            {"descriptor__subject_information__subject_id__icontains": "P-00"},
            [self.entities[0]],
        )

        # Descriptor / Subject Information / Batch
        self._check_filter(
            {"descriptor__subject_information__batch__exact": "1"},
            [self.entities[0]],
        )
        self._check_filter(
            {"descriptor__subject_information__batch__exact": "xyz"},
            None,
            expected_status_code=status.HTTP_400_BAD_REQUEST,
        )

        # Descriptor / Subject Information / Group
        self._check_filter(
            {"descriptor__subject_information__group__iexact": "pre"},
            [self.entities[0]],
        )

        # Descriptor / Disease Information / Deisease Type
        self._check_filter(
            {
                "descriptor__disease_information__disease_type__icontains": "Colorectal cancer"
            },
            [self.entities[0]],
        )
        self._check_filter(
            {"descriptor__disease_information__disease_type__icontains": "Cancer"},
            [self.entities[0]],
        )

        # Descriptor / Disease Information / Disease Status
        self._check_filter(
            {"descriptor__disease_information__disease_status__iexact": "progresive"},
            [self.entities[0]],
        )

        # Descriptor / Immuno Oncology Treatment Type / IO Drug
        self._check_filter(
            {"descriptor__immuno_oncology_treatment_type__io_drug__iexact": "D-00A"},
            [self.entities[0]],
        )

        # Descriptor / Immuno Oncology Treatment Type / IO Treatment
        self._check_filter(
            {
                "descriptor__immuno_oncology_treatment_type__io_treatment__iexact": "single"
            },
            [self.entities[0]],
        )

        # Descriptor / Response and Survival Analysis / Confirmed BOR
        self._check_filter(
            {"descriptor__response_and_survival_analysis__confirmed_bor__iexact": "pd"},
            [self.entities[0]],
        )

        # Descriptor / Response and Survival Analysis / PFS Event
        self._check_filter(
            {"descriptor__response_and_survival_analysis__pfs_event__iexact": "NO"},
            [self.entities[0]],
        )

        # Descriptor / General / Description
        self._check_filter(
            {"descriptor__general__description__icontains": "First sample"},
            [self.entities[2]],
        )
        self._check_filter(
            {"descriptor__general__description__icontains": "first"},
            [self.entities[2]],
        )
        self._check_filter(
            {"descriptor__general__description__icontains": "sample"},
            self.entities[2:],
        )

        # Descriptor / General / Biosample source
        self._check_filter(
            {"descriptor__general__biosample_source__icontains": "lung"},
            [self.entities[2]],
        )

        # Descriptor / General / Biosample Treatment
        self._check_filter(
            {"descriptor__general__biosample_treatment__icontains": "dmso"},
            [self.entities[2]],
        )
