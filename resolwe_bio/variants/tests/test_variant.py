from typing import TYPE_CHECKING
from unittest.mock import Mock

from django.contrib.auth import get_user_model
from django.test import TestCase
from rest_framework import status
from rest_framework.test import APIRequestFactory, force_authenticate

from resolwe.flow.executors.socket_utils import Message, MessageType
from resolwe.flow.models import Data
from resolwe.flow.models import Entity as Sample
from resolwe.flow.models import Process
from resolwe.permissions.models import Permission

from resolwe_bio.variants.listener_plugin import (
    VariantAnnotationData,
    VariantCommands,
    VariantData,
)
from resolwe_bio.variants.models import (
    Variant,
    VariantAnnotation,
    VariantAnnotationTranscript,
    VariantCall,
    VariantExperiment,
)
from resolwe_bio.variants.views import (
    VariantAnnotationSerializer,
    VariantAnnotationViewSet,
    VariantCallSerializer,
    VariantCallViewSet,
    VariantExperimentSerializer,
    VariantExperimentViewSet,
    VariantSerializer,
    VariantViewSet,
)

if TYPE_CHECKING:
    from django.contrib.auth.models import User
else:
    User = get_user_model()


class PrepareDataMixin:
    """Prepare the data for all variant tests."""

    variants: list[Variant]
    annotations: list[VariantAnnotation]
    experiments: list[VariantExperiment]
    calls: list[VariantCall]
    contributor: User

    @classmethod
    def setUpTestData(cls):
        """Set up the test data."""
        cls.contributor = User.objects.get_or_create(
            username="contributor", email="contributor@genialis.com"
        )[0]
        sample = Sample.objects.create(contributor=cls.contributor)
        cls.variants = Variant.objects.bulk_create(
            [
                Variant(
                    species="Homo Sapiens",
                    genome_assembly="assembly 1",
                    chromosome="CHR1",
                    position=1234,
                    reference="ref1",
                    alternative="alt1",
                ),
                Variant(
                    species="Mus Musculus",
                    genome_assembly="assembly 2",
                    chromosome="CHR2",
                    position=67891234,
                    reference="ref2",
                    alternative="alt2",
                ),
            ]
        )
        cls.annotations = VariantAnnotation.objects.bulk_create(
            [
                VariantAnnotation(
                    variant=cls.variants[0],
                    type="SNP",
                    clinical_diagnosis="clinical diagnosis 1",
                    clinical_significance="clinical significance 1",
                    dbsnp_id="dbsnp_id 1",
                    clinvar_id="clinical_var_id 1",
                )
            ]
        )
        # Create the transcript data.
        VariantAnnotationTranscript.objects.bulk_create(
            [
                VariantAnnotationTranscript(
                    variant_annotation=cls.annotations[0],
                    annotation="annotation 1",
                    annotation_impact="impact 1",
                    gene="gene 1",
                    protein_impact="protein impact 1",
                    transcript_ids=["f1", "f2"],
                )
            ]
        )

        cls.experiments = VariantExperiment.objects.bulk_create(
            [
                VariantExperiment(
                    variant_data_source="source 1",
                    contributor=cls.contributor,
                ),
                VariantExperiment(
                    variant_data_source="source 2",
                    contributor=cls.contributor,
                ),
            ]
        )
        cls.calls = VariantCall.objects.bulk_create(
            [
                VariantCall(
                    sample=sample,
                    variant=cls.variants[0],
                    quality=0.7,
                    depth_norm_quality=0.7,
                    unfiltered_allele_depth=1,
                    depth=15,
                    filter="filter 1",
                    genotype="1",
                    genotype_quality=1,
                    experiment=cls.experiments[0],
                ),
                VariantCall(
                    sample=sample,
                    variant=cls.variants[1],
                    quality=0.2,
                    depth_norm_quality=0.2,
                    unfiltered_allele_depth=2,
                    depth=5,
                    filter="filter 2",
                    genotype="2",
                    genotype_quality=2,
                    experiment=cls.experiments[1],
                ),
            ]
        )


class ListenerPluginTest(TestCase):

    def setUp(self):
        """Prepare the test data."""
        contributor = User.objects.get_or_create(
            username="contributor", email="contributor@genialis.com"
        )[0]
        self.data = Data.objects.create(
            contributor=contributor,
            entity=Sample.objects.create(contributor=contributor),
            process=Process.objects.create(contributor=contributor),
            status=Data.STATUS_PROCESSING,
        )
        Variant.objects.create(
            species="Homo Sapiens",
            genome_assembly="ENSEMBL",
            chromosome="chr1",
            position=1,
            reference="ref1",
            alternative="alt1",
        )
        return super().setUp()

    def test_add_variants(self):
        """Test listener method."""

        # The first variant already exists in the database and should be re-used.
        # The second one is new and should be created.
        data_source = "process"

        variants_data: list[VariantData] = [
            {
                # Identity the variant.
                "species": "Homo Sapiens",
                "genome_assembly": "ENSEMBL",
                "chromosome": "chr1",
                "position": 1,
                "reference": "ref1",
                "alternative": "alt1",
                # Variant annotation data.
                "quality": 1,
                "depth": 1,
                "genotype": "1",
                "filter": "1",
            },
            {
                # Variant data.
                "species": "Homo Sapiens",
                "genome_assembly": "ENSEMBL",
                "chromosome": "chr2",
                "position": 2,
                "reference": "ref2",
                "alternative": "alt2",
                # Variant annotation data.
                "quality": 2,
                "depth_norm_quality": 0.2,
                "unfiltered_allele_depth": 2,
                "depth": 2,
                "genotype": "2",
                "genotype_quality": 2,
                "filter": "2",
            },
        ]
        message = Message(
            MessageType.COMMAND, "variants_test", [data_source, variants_data]
        )
        manager_mock = Mock(data=Mock(return_value=self.data))
        VariantCommands().handle_add_variants(
            data_id=self.data.pk,
            message=message,
            manager=manager_mock,
        )
        expected = [
            {
                "species": "Homo Sapiens",
                "genome_assembly": "ENSEMBL",
                "chromosome": "chr1",
                "position": 1,
                "reference": "ref1",
                "alternative": "alt1",
            },
            {
                "species": "Homo Sapiens",
                "genome_assembly": "ENSEMBL",
                "chromosome": "chr2",
                "position": 2,
                "reference": "ref2",
                "alternative": "alt2",
            },
        ]

        self.assertCountEqual(
            Variant.objects.all().values(
                "species",
                "genome_assembly",
                "chromosome",
                "position",
                "reference",
                "alternative",
            ),
            expected,
        )
        # Exactly one experiment should be created.
        experiment = VariantExperiment.objects.get()
        self.assertEqual(experiment.contributor, self.data.contributor)
        self.assertEqual(experiment.variant_data_source, "process")

        variant1 = Variant.objects.get(chromosome="chr1")
        variant2 = Variant.objects.get(chromosome="chr2")
        expected_calls = [
            {
                "sample_id": self.data.entity.pk,
                "variant_id": variant1.pk,
                "experiment_id": experiment.pk,
                "quality": 1.0,
                "depth_norm_quality": None,
                "unfiltered_allele_depth": None,
                "depth": 1,
                "filter": "1",
                "genotype": "1",
                "genotype_quality": None,
                "data_id": self.data.pk,
            },
            {
                "sample_id": self.data.entity.pk,
                "variant_id": variant2.pk,
                "experiment_id": experiment.pk,
                "quality": 2.0,
                "depth_norm_quality": 0.2,
                "unfiltered_allele_depth": 2,
                "depth": 2,
                "filter": "2",
                "genotype": "2",
                "genotype_quality": 2,
                "data_id": self.data.pk,
            },
        ]
        self.assertCountEqual(
            VariantCall.objects.all().values(
                "sample_id",
                "variant_id",
                "experiment_id",
                "quality",
                "depth_norm_quality",
                "unfiltered_allele_depth",
                "depth",
                "filter",
                "genotype",
                "genotype_quality",
                "data_id",
            ),
            expected_calls,
        )

    def test_add_variant_annotations(self):

        annotations: list[VariantAnnotationData] = [
            {
                # Identity the variant.
                "species": "Homo Sapiens",
                "genome_assembly": "ENSEMBL",
                "chromosome": "chr1",
                "position": 1,
                "reference": "ref1",
                "alternative": "alt1",
                # Variant annotation data.
                "type": "SNP",
                "clinical_diagnosis": "diagnosis 1",
                "clinical_significance": "significance 1",
                "dbsnp_id": "dbsnp 1",
                "clinvar_id": "clinvar 1",
                "transcripts": [
                    {
                        "annotation": "annotation 1",
                        "annotation_impact": "annotation impact 1",
                        "gene": "gene 1",
                        "protein_impact": "protein impact 1",
                        "transcript_ids": ["t 1", "tt 1"],
                        "canonical": False,
                    }
                ],
            },
            {
                # Identity the variant.
                "species": "Homo Sapiens",
                "genome_assembly": "ENSEMBL",
                "chromosome": "chr2",
                "position": 2,
                "reference": "ref2",
                "alternative": "alt2",
                # Variant annotation data.
                "type": "SNP",
                "clinical_diagnosis": "diagnosis 2",
                "clinical_significance": "significance 2",
                "dbsnp_id": "dbsnp 2",
                "clinvar_id": "clinvar 2",
                "transcripts": [
                    {
                        "annotation": "annotation 2",
                        "annotation_impact": "annotation impact 2",
                        "gene": "gene 2",
                        "protein_impact": "protein impact 2",
                        "transcript_ids": ["t 2", "tt 2"],
                        "canonical": True,
                    }
                ],
            },
        ]
        message = Message(MessageType.COMMAND, "variants_test", annotations)
        manager_mock = Mock(data=Mock(return_value=self.data))
        VariantCommands().handle_add_variants_annotations(
            data_id=self.data.pk,
            message=message,
            manager=manager_mock,
        )
        expected = [
            {
                "species": "Homo Sapiens",
                "genome_assembly": "ENSEMBL",
                "chromosome": "chr1",
                "position": 1,
                "reference": "ref1",
                "alternative": "alt1",
            },
            {
                "species": "Homo Sapiens",
                "genome_assembly": "ENSEMBL",
                "chromosome": "chr2",
                "position": 2,
                "reference": "ref2",
                "alternative": "alt2",
            },
        ]

        self.assertCountEqual(
            Variant.objects.all().values(
                "species",
                "genome_assembly",
                "chromosome",
                "position",
                "reference",
                "alternative",
            ),
            expected,
        )

        # Now checks for annotations.
        variant1, variant2 = Variant.objects.all()

        expected_annotation_1 = {
            "id": variant1.annotation.id,
            "variant_id": variant1.id,
            "type": "SNP",
            "clinical_diagnosis": "diagnosis 1",
            "clinical_significance": "significance 1",
            "dbsnp_id": "dbsnp 1",
            "clinvar_id": "clinvar 1",
            "transcripts": [
                {
                    "id": variant1.annotation.transcripts.get().id,
                    "variant_annotation_id": variant1.annotation.pk,
                    "annotation": "annotation 1",
                    "annotation_impact": "annotation impact 1",
                    "gene": "gene 1",
                    "protein_impact": "protein impact 1",
                    "transcript_ids": ["t 1", "tt 1"],
                    "canonical": False,
                }
            ],
        }
        expected_annotation_2 = {
            "id": variant2.annotation.id,
            "variant_id": variant2.id,
            "type": "SNP",
            "clinical_diagnosis": "diagnosis 2",
            "clinical_significance": "significance 2",
            "dbsnp_id": "dbsnp 2",
            "clinvar_id": "clinvar 2",
            "transcripts": [
                {
                    "id": variant2.annotation.transcripts.get().id,
                    "variant_annotation_id": variant2.annotation.pk,
                    "annotation": "annotation 2",
                    "annotation_impact": "annotation impact 2",
                    "gene": "gene 2",
                    "protein_impact": "protein impact 2",
                    "transcript_ids": ["t 2", "tt 2"],
                    "canonical": True,
                }
            ],
        }

        self.assertDictEqual(
            VariantAnnotationSerializer(variant1.annotation).data, expected_annotation_1
        )
        self.assertDictEqual(
            VariantAnnotationSerializer(variant2.annotation).data, expected_annotation_2
        )


class VariantTest(PrepareDataMixin, TestCase):
    def setUp(self) -> None:
        self.view = VariantViewSet.as_view({"get": "list"})
        return super().setUp()

    def test_filter(self):
        """Test the Variant filter."""
        request = APIRequestFactory().get("/variant")
        # Basic get, no filter.
        expected = VariantSerializer(self.variants, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by id.
        request = APIRequestFactory().get("/variant", {"id": self.variants[0].id})
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by species.
        request = APIRequestFactory().get("/variant", {"species": "Homo Sapiens"})
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by genome_assembly.
        request = APIRequestFactory().get(
            "/variant", {"genome_assembly__icontains": "Embly 1"}
        )
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by chromosome.
        request = APIRequestFactory().get("/variant", {"chromosome__iexact": "chr1"})
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by position.
        request = APIRequestFactory().get("/variant", {"position__lt": "12345"})
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by reference.
        request = APIRequestFactory().get("/variant", {"reference": "ref1"})
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by alternative.
        request = APIRequestFactory().get("/variant", {"alternative": "alt1"})
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by annotation type.
        request = APIRequestFactory().get(
            "/variant", {"annotation__type__isnull": "true"}
        )
        expected = VariantSerializer(self.variants[1:], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variant", {"annotation__type__contains": "SN"}
        )
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by annotation.
        request = APIRequestFactory().get(
            "/variant", {"annotation__transcripts__annotation": "annotation 1"}
        )
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variant", {"annotation__transcripts__annotation__isnull": True}
        )
        expected = VariantSerializer(self.variants[1:], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by annotation impact.
        request = APIRequestFactory().get(
            "/variant",
            {"annotation__transcripts__annotation_impact__icontains": "MpAcT 1"},
        )
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variant", {"annotation__transcripts__annotation_impact__isnull": True}
        )
        expected = VariantSerializer(self.variants[1:], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by gene.
        request = APIRequestFactory().get(
            "/variant", {"annotation__transcripts__gene": "gene 1"}
        )
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variant", {"annotation__transcripts__gene__isnull": True}
        )
        expected = VariantSerializer(self.variants[1:], many=True).data
        response = self.view(request)

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by protein impact.
        request = APIRequestFactory().get(
            "/variant",
            {"annotation__transcripts__protein_impact__icontains": "protein"},
        )
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variant", {"annotation__transcripts__protein_impact__isnull": True}
        )
        expected = VariantSerializer(self.variants[1:], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by clinical diagnosis.
        request = APIRequestFactory().get(
            "/variant",
            {"annotation__clinical_diagnosis__icontains": "clinical diagnosis"},
        )
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variant", {"annotation__clinical_diagnosis__isnull": True}
        )
        expected = VariantSerializer(self.variants[1:], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by clinical significance.
        request = APIRequestFactory().get(
            "/variant",
            {"annotation__clinical_significance__contains": "significance 1"},
        )
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variant", {"annotation__clinical_significance__isnull": True}
        )
        expected = VariantSerializer(self.variants[1:], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by quality.
        request = APIRequestFactory().get(
            "/variant", {"variant_calls__quality__isnull": "True"}
        )
        expected = VariantSerializer([], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variant", {"variant_calls__quality__gt": 0.5}
        )
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by depth.
        request = APIRequestFactory().get(
            "/variant", {"variant_calls__depth__isnull": "True"}
        )
        expected = VariantSerializer([], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get("/variant", {"variant_calls__depth__gt": 10})
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

    def test_ordering(self):
        """Test the Variant ordering."""
        # Order by species.
        request = APIRequestFactory().get("/variant", {"ordering": "species"})
        expected = VariantSerializer(self.variants, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get("/variant", {"ordering": "-species"})
        expected = VariantSerializer(self.variants[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)

        # Order by position.
        request = APIRequestFactory().get("/variant", {"ordering": "position"})
        expected = VariantSerializer(self.variants, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get("/variant", {"ordering": "-position"})
        expected = VariantSerializer(self.variants[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)

        # Order by genome_assembly.
        request = APIRequestFactory().get("/variant", {"ordering": "genome_assembly"})
        expected = VariantSerializer(self.variants, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get("/variant", {"ordering": "-genome_assembly"})
        expected = VariantSerializer(self.variants[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)

        # Order by chromosome.
        request = APIRequestFactory().get("/variant", {"ordering": "chromosome"})
        expected = VariantSerializer(self.variants, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get("/variant", {"ordering": "-chromosome"})
        expected = VariantSerializer(self.variants[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)


class VariantAnnotationTest(PrepareDataMixin, TestCase):
    def setUp(self) -> None:
        self.view = VariantAnnotationViewSet.as_view({"get": "list"})
        self.annotations.append(
            VariantAnnotation.objects.create(
                variant=self.variants[1],
                type="INDEL",
                clinical_diagnosis="clinical diagnosis 2",
                clinical_significance="clinical significance 2",
                dbsnp_id="dbsnp_id 2",
                clinvar_id="clinical_var_id 2",
            )
        )
        VariantAnnotationTranscript.objects.create(
            variant_annotation=self.annotations[-1],
            annotation="annotation 2",
            annotation_impact="impact 2",
            gene="gene 2",
            protein_impact="protein impact 2",
            transcript_ids=["f1", "f2"],
        )
        return super().setUp()

    def test_filter(self):
        # No filter.
        request = APIRequestFactory().get("/variantannotation")
        expected = VariantAnnotationSerializer(self.annotations, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by id.
        request = APIRequestFactory().get(
            "/variantannotation", {"id": self.annotations[0].id}
        )
        expected = VariantAnnotationSerializer(self.annotations[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by type.
        request = APIRequestFactory().get("/variantannotation", {"type": "SNP"})
        expected = VariantAnnotationSerializer(self.annotations[:1], many=True).data
        response = self.view(request)

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by annotation.
        request = APIRequestFactory().get(
            "/variantannotation", {"transcripts__annotation": "annotation 1"}
        )
        expected = VariantAnnotationSerializer(self.annotations[:1], many=True).data
        response = self.view(request)

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by annotation impact.
        request = APIRequestFactory().get(
            "/variantannotation", {"transcripts__annotation_impact": "impact 1"}
        )
        expected = VariantAnnotationSerializer(self.annotations[:1], many=True).data
        response = self.view(request)

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by gene.
        request = APIRequestFactory().get(
            "/variantannotation", {"transcripts__gene": "gene 1"}
        )
        expected = VariantAnnotationSerializer(self.annotations[:1], many=True).data
        response = self.view(request)

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by protein impact.
        request = APIRequestFactory().get(
            "/variantannotation", {"transcripts__protein_impact": "protein impact 1"}
        )
        expected = VariantAnnotationSerializer(self.annotations[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by clinical diagnosis.
        request = APIRequestFactory().get(
            "/variantannotation", {"clinical_diagnosis": "clinical diagnosis 1"}
        )
        expected = VariantAnnotationSerializer(self.annotations[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by clinical significance.
        request = APIRequestFactory().get(
            "/variantannotation", {"clinical_significance": "clinical significance 1"}
        )
        expected = VariantAnnotationSerializer(self.annotations[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by dbsnp_id.
        request = APIRequestFactory().get(
            "/variantannotation", {"dbsnp_id": "dbsnp_id 1"}
        )
        expected = VariantAnnotationSerializer(self.annotations[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by clinical_var_id.
        request = APIRequestFactory().get(
            "/variantannotation", {"clinvar_id": "clinical_var_id 1"}
        )
        expected = VariantAnnotationSerializer(self.annotations[:1], many=True).data
        response = self.view(request)

        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

    def test_ordering(self):
        # Order by gene.
        request = APIRequestFactory().get(
            "/variantannotation", {"ordering": "transcripts__gene"}
        )
        expected = VariantAnnotationSerializer(self.annotations, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variantannotation", {"ordering": "-transcripts__gene"}
        )
        expected = VariantAnnotationSerializer(self.annotations[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)

        # Order by protein impact.
        request = APIRequestFactory().get(
            "/variantannotation", {"ordering": "transcripts__protein_impact"}
        )
        expected = VariantAnnotationSerializer(self.annotations, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variantannotation", {"ordering": "-transcripts__protein_impact"}
        )
        expected = VariantAnnotationSerializer(self.annotations[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)

        # Sort by annotation.
        request = APIRequestFactory().get(
            "/variantannotation", {"ordering": "transcripts__annotation"}
        )
        expected = VariantAnnotationSerializer(self.annotations, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variantannotation", {"ordering": "-transcripts__annotation"}
        )
        expected = VariantAnnotationSerializer(self.annotations[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)

        # Sort by clinical significance.
        request = APIRequestFactory().get(
            "/variantannotation", {"ordering": "clinical_significance"}
        )
        expected = VariantAnnotationSerializer(self.annotations, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variantannotation", {"ordering": "-clinical_significance"}
        )
        expected = VariantAnnotationSerializer(self.annotations[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)


class VariantCallTest(PrepareDataMixin, TestCase):
    def setUp(self) -> None:
        self.view = VariantCallViewSet.as_view({"get": "list"})
        return super().setUp()

    def test_filter(self):
        # No filter no permission for public.
        request = APIRequestFactory().get("/variantcall")
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, [])

        # No filter no permissions for contributor.
        request = APIRequestFactory().get("/variantcall")
        force_authenticate(request, self.contributor)
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, [])

        # No filter read permissions for contributor.
        for call in self.calls:
            call.sample.set_permission(Permission.VIEW, self.contributor)
        request = APIRequestFactory().get("/variantcall")
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by id.
        request = APIRequestFactory().get("/variantcall", {"id": self.calls[0].id})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by quality.
        request = APIRequestFactory().get("/variantcall", {"quality__gt": 0.5})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by depth.
        request = APIRequestFactory().get("/variantcall", {"depth__gt": 10})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by variant id.
        request = APIRequestFactory().get(
            "/variantcall", {"variant__id": self.variants[0].id}
        )
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by variant species.
        request = APIRequestFactory().get(
            "/variantcall", {"variant__species": self.variants[0].species}
        )
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by variant genome_assembly.
        request = APIRequestFactory().get(
            "/variantcall",
            {"variant__genome_assembly": self.variants[0].genome_assembly},
        )
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by variant chromosome.
        request = APIRequestFactory().get(
            "/variantcall", {"variant__chromosome": self.variants[0].chromosome}
        )
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by variant position.
        request = APIRequestFactory().get(
            "/variantcall",
            {
                "variant__position__gte": self.variants[0].position,
                "variant__position__lt": self.variants[1].position,
            },
        )
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by variant reference.
        request = APIRequestFactory().get(
            "/variantcall", {"variant__reference": self.variants[0].reference}
        )
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by variant alternative.
        request = APIRequestFactory().get(
            "/variantcall", {"variant__alternative": self.variants[0].alternative}
        )
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by experiment id.
        request = APIRequestFactory().get(
            "/variantcall", {"experiment__id": self.experiments[0].id}
        )
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

    def test_ordering(self):
        # Set permissions for the contributor.
        for call in self.calls:
            call.sample.set_permission(Permission.VIEW, self.contributor)

        # Order by id.
        request = APIRequestFactory().get("/variantcall", {"ordering": "id"})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get("/variantcall", {"ordering": "-id"})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)

        # Order by quality.
        request = APIRequestFactory().get("/variantcall", {"ordering": "-quality"})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get("/variantcall", {"ordering": "quality"})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)

        # Order by depth.
        request = APIRequestFactory().get("/variantcall", {"ordering": "-depth"})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get("/variantcall", {"ordering": "depth"})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)


class VariantExperimentTest(PrepareDataMixin, TestCase):
    def setUp(self) -> None:
        self.view = VariantExperimentViewSet.as_view({"get": "list"})
        return super().setUp()

    def test_filter(self):
        # No filter.
        request = APIRequestFactory().get("/variantexperiment")
        expected = VariantExperimentSerializer(self.experiments, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by id.
        request = APIRequestFactory().get(
            "/variantexperiment", {"id": self.experiments[0].id}
        )
        expected = VariantExperimentSerializer(self.experiments[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by variant data source.
        request = APIRequestFactory().get(
            "/variantexperiment", {"variant_data_source": "source 1"}
        )
        expected = VariantExperimentSerializer(self.experiments[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by contributor.
        contributor2 = User.objects.create(
            username="contributor2", email="contributor2@genialis.com"
        )
        self.experiments[1].contributor = contributor2
        self.experiments[1].save()
        request = APIRequestFactory().get(
            "/variantexperiment", {"contributor__id": self.contributor.id}
        )
        expected = VariantExperimentSerializer(self.experiments[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variantexperiment", {"contributor__username": self.contributor.username}
        )
        expected = VariantExperimentSerializer(self.experiments[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variantexperiment", {"contributor__email": self.contributor.email}
        )
        expected = VariantExperimentSerializer(self.experiments[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by date.
        request = APIRequestFactory().get(
            "/variantexperiment", {"timestamp__lt": self.experiments[1].timestamp}
        )
        expected = VariantExperimentSerializer(self.experiments[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

    def test_ordering(self):
        # Order by id.
        request = APIRequestFactory().get("/variantexperiment", {"ordering": "id"})
        expected = VariantExperimentSerializer(self.experiments, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get("/variantexperiment", {"ordering": "-id"})
        expected = VariantExperimentSerializer(self.experiments[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)

        # Order by date.
        request = APIRequestFactory().get(
            "/variantexperiment", {"ordering": "timestamp"}
        )
        expected = VariantExperimentSerializer(self.experiments, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variantexperiment", {"ordering": "-timestamp"}
        )
        expected = VariantExperimentSerializer(self.experiments[::-1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)

        # Order by contributor.
        contributor2 = User.objects.create(
            username="contributor2", email="contributor2@genialis.com"
        )
        self.experiments[1].contributor = contributor2
        self.experiments[1].save()

        request = APIRequestFactory().get(
            "/variantexperiment", {"ordering": "-contributor__email"}
        )
        expected = VariantExperimentSerializer(self.experiments, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected)
        request = APIRequestFactory().get(
            "/variantexperiment", {"ordering": "contributor__email"}
        )
        response = self.view(
            APIRequestFactory().get(
                "/variantexperiment", {"ordering": "contributor__email"}
            )
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, expected[::-1])
