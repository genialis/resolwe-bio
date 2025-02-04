import os
from typing import TYPE_CHECKING
from unittest.mock import Mock

from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APIClient, APIRequestFactory, force_authenticate

from resolwe.flow.executors.socket_utils import Message, MessageType
from resolwe.flow.models import Data
from resolwe.flow.models import Entity as Sample
from resolwe.flow.models import Process
from resolwe.permissions.models import Permission
from resolwe.test import ProcessTestCase, tag_process

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
        cls.sample = Sample.objects.create(
            contributor=cls.contributor, name="Test sample"
        )
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
                    transcript_id="f1",
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
                    sample=cls.sample,
                    variant=cls.variants[0],
                    quality=0.7,
                    depth_norm_quality=0.7,
                    alternative_allele_depth=1,
                    depth=15,
                    filter="filter 1",
                    genotype="1",
                    genotype_quality=1,
                    experiment=cls.experiments[0],
                ),
                VariantCall(
                    sample=cls.sample,
                    variant=cls.variants[1],
                    quality=0.2,
                    depth_norm_quality=0.2,
                    alternative_allele_depth=2,
                    depth=5,
                    filter="filter 2",
                    genotype="2",
                    genotype_quality=2,
                    experiment=cls.experiments[1],
                ),
            ]
        )


class VariantProcessTest(ProcessTestCase):
    """Test Variant model can be accesed inside Python process."""

    def setUp(self):
        """Register the processes and prepare necessary objects."""
        super().setUp()
        self._register_schemas(processes_paths=[os.path.dirname(__file__)])
        self.contributor = User.objects.get_or_create(
            username="contributor_variants", email="contributor_variants@genialis.com"
        )[0]

        self.sample1 = Sample.objects.create(contributor=self.contributor)
        self.sample2 = Sample.objects.create(contributor=self.contributor)

        variants = Variant.objects.bulk_create(
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

        VariantCall.objects.bulk_create(
            [
                VariantCall(
                    sample=self.sample1,
                    variant=variants[0],
                    quality=0.7,
                    depth_norm_quality=0.7,
                    alternative_allele_depth=1,
                    depth=15,
                    filter="filter 1",
                    genotype="1",
                    genotype_quality=1,
                ),
                VariantCall(
                    sample=self.sample2,
                    variant=variants[1],
                    quality=0.2,
                    depth_norm_quality=0.2,
                    alternative_allele_depth=2,
                    depth=5,
                    filter="filter 2",
                    genotype="2",
                    genotype_quality=2,
                ),
            ]
        )

    @tag_process("test-variants-simple-get")
    def test_get_variants(self):
        """Test retrieving the Variant model in Python process."""
        data = self.run_process("test-variants-simple-get")
        self.assertCountEqual(
            data.output["variant_species"], ["Mus Musculus", "Homo Sapiens"]
        )

    @tag_process("test-variants-simple-filter")
    def test_filter_variants(self):
        """Test filtering Variant model in Python process."""
        data = self.run_process(
            "test-variants-simple-filter", {"species_filter": "Homo Sapiens"}
        )
        self.assertCountEqual(data.output["variant_species"], ["Homo Sapiens"])

    @tag_process("test-variant-position-get")
    def test_get_variant_position(self):
        """Test retrieving Variant position in Python process."""
        variant = Variant.objects.last()
        data = self.run_process("test-variant-position-get", {"variant_id": variant.pk})
        self.assertEqual(data.output["position"], variant.position)

    @tag_process("test-variantcalls-simple-get")
    def test_get_variant_calls(self):
        """Test retrieving the VariantCall models in Python process."""

        Process.objects.get(slug="test-variantcalls-simple-get").set_permission(
            Permission.VIEW, self.contributor
        )
        data = self.run_process(
            "test-variantcalls-simple-get", contributor=self.contributor
        )
        self.assertCountEqual(data.output["filters"], [])

        self.sample1.set_permission(Permission.VIEW, self.contributor)
        data = self.run_process(
            "test-variantcalls-simple-get", contributor=self.contributor
        )
        self.assertCountEqual(data.output["filters"], ["filter 1"])

        self.sample2.set_permission(Permission.VIEW, self.contributor)
        data = self.run_process(
            "test-variantcalls-simple-get", contributor=self.contributor
        )
        self.assertCountEqual(data.output["filters"], ["filter 1", "filter 2"])


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
                "alternative_allele_depth": 2,
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
                "alternative_allele_depth": None,
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
                "alternative_allele_depth": 2,
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
                "alternative_allele_depth",
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
                        "transcript_id": "t 1",
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
                        "transcript_id": "t 2",
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
                    "transcript_id": "t 1",
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
                    "transcript_id": "t 2",
                    "canonical": True,
                }
            ],
        }

        print(VariantAnnotationSerializer(variant1.annotation).data)

        self.assertDictEqual(
            VariantAnnotationSerializer(variant1.annotation).data, expected_annotation_1
        )
        self.assertDictEqual(
            VariantAnnotationSerializer(variant2.annotation).data, expected_annotation_2
        )


class VariantTest(PrepareDataMixin, TestCase):
    def setUp(self) -> None:
        self.view = VariantViewSet.as_view(
            {
                "get": "list",
                "post": "create",
                "delete": "destroy",
                "patch": "partial_update",
            }
        )
        return super().setUp()

    def test_nested_serialization(self):
        """Test that the annotation data is serialized correctly.

        It must be serialized inside the variant data and include all the transcripts.
        """
        request = APIRequestFactory().get("/variant")
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), 2)
        self.assertTrue("id" in response.data[0]["annotation"]["transcripts"][0])

    def test_create_api(self):
        """Test the Variant creation using the API endpoint."""
        post_data = {
            "species": "Homo sapiens",
            "genome_assembly": "assembly 3",
            "chromosome": "chr3",
            "position": 3,
            "reference": "ref3",
            "alternative": "alt3",
        }
        variant_count = Variant.objects.count()

        # Unauthenticated request must be denied.
        request = APIRequestFactory().post("/variant", post_data, format="json")
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)
        self.assertEqual(
            response.data, {"detail": "Authentication credentials were not provided."}
        )

        # Request from non-staff user must be denied.
        force_authenticate(request, self.contributor)
        response = self.view(request)
        self.assertContains(
            response,
            "You do not have permission to perform this action.",
            status_code=status.HTTP_403_FORBIDDEN,
        )

        # Request from staff users must be granted.
        self.contributor.is_staff = True
        self.contributor.save(update_fields=["is_staff"])
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(Variant.objects.count(), variant_count + 1)
        variant = Variant.objects.last()
        self.assertEqual(variant.species, "Homo sapiens")
        self.assertEqual(variant.genome_assembly, "assembly 3")
        self.assertEqual(variant.chromosome, "chr3")
        self.assertEqual(variant.position, 3)
        self.assertEqual(variant.reference, "ref3")
        self.assertEqual(variant.alternative, "alt3")

        # Test nested objects creation (annotation & transcripts).
        post_data["position"] = 4
        post_data["annotation"] = {
            "type": "SNP",
            "slug": "test_slug",
            "clinical_diagnosis": "clinical diagnosis 1",
            "clinical_significance": "clinical significance 1",
            "dbsnp_id": "dbsnp_id 1",
            "clinvar_id": "clinical_var_id 1",
            "transcripts": [
                {
                    "annotation": "created 1",
                    "annotation_impact": "impact created 1",
                    "gene": "gene created 1",
                    "protein_impact": "protein impact created 1",
                    "transcript_id": "f1",
                    "canonical": True,
                }
            ],
        }

        request = APIRequestFactory().post("/variant", post_data, format="json")
        force_authenticate(request, self.contributor)
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(Variant.objects.count(), variant_count + 2)
        variant = Variant.objects.last()
        annotation = VariantAnnotation.objects.last()
        transcript = VariantAnnotationTranscript.objects.last()
        self.assertEqual(variant.species, "Homo sapiens")
        self.assertEqual(variant.annotation, annotation)
        self.assertTrue(annotation.transcripts.filter(pk=transcript.pk).exists())

        self.contributor.is_staff = False
        self.contributor.save(update_fields=["is_staff"])

    def test_delete_api(self):
        """Test staff user can delete variants."""
        variant_count = Variant.objects.count()
        variant_to_delete = Variant.objects.first()
        # Unauthenticated request must be denied.
        request = APIRequestFactory().delete("/variant")
        response = self.view(request, pk=variant_to_delete.pk)
        self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)
        self.assertEqual(
            response.data, {"detail": "Authentication credentials were not provided."}
        )

        # Request from non-staff user must be denied.
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=variant_to_delete.pk)
        self.assertContains(
            response,
            "You do not have permission to perform this action.",
            status_code=status.HTTP_403_FORBIDDEN,
        )

        # Request from staff users must be granted.
        self.contributor.is_staff = True
        self.contributor.save(update_fields=["is_staff"])
        response = self.view(request, pk=variant_to_delete.pk)
        self.assertEqual(response.status_code, status.HTTP_204_NO_CONTENT)
        self.assertEqual(Variant.objects.count(), variant_count - 1)
        with self.assertRaises(Variant.DoesNotExist):
            variant_to_delete.refresh_from_db()
        self.contributor.is_staff = False
        self.contributor.save(update_fields=["is_staff"])

    def test_patch_api(self):
        """Test updating variant objects over api."""
        # Unauthenticated request must be denied.
        variant = Variant.objects.first()
        request = APIRequestFactory().patch("/variant")
        response = self.view(request, pk=variant.pk)
        self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)
        self.assertEqual(
            response.data, {"detail": "Authentication credentials were not provided."}
        )

        # Request from non-staff user must be denied.
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=variant.pk)
        self.assertContains(
            response,
            "You do not have permission to perform this action.",
            status_code=status.HTTP_403_FORBIDDEN,
        )

        # Request from staff users must be granted.
        self.contributor.is_staff = True
        self.contributor.save(update_fields=["is_staff"])

        patch_data: dict = {"species": "new species"}
        request = APIRequestFactory().patch("/variant", patch_data, format="json")
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=variant.pk)
        variant.refresh_from_db()
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(variant.species, "new species")

        # Update existing annotation.
        patch_data = {
            "annotation": {
                "type": "INDEL",
                "clinical_diagnosis": "new diagnosis",
                "clinical_significance": "new significance",
                "dbsnp_id": "new dbsnp",
                "clinvar_id": "new clinvar",
                "transcripts": [
                    {
                        "annotation": "new annotation",
                        "annotation_impact": "new impact",
                        "gene": "new gene",
                        "protein_impact": "new protein impact",
                        "transcript_id": "new transcript",
                        "canonical": False,
                    }
                ],
            }
        }
        variant_count = Variant.objects.count()
        request = APIRequestFactory().patch("/variant", patch_data, format="json")
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=variant.pk)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(Variant.objects.count(), variant_count)
        variant.refresh_from_db()
        self.assertEqual(variant.annotation.type, "INDEL")
        self.assertEqual(variant.annotation.clinical_diagnosis, "new diagnosis")
        self.assertEqual(variant.annotation.clinical_significance, "new significance")
        self.assertEqual(variant.annotation.dbsnp_id, "new dbsnp")
        self.assertEqual(variant.annotation.clinvar_id, "new clinvar")
        self.assertEqual(variant.annotation.transcripts.count(), 1)
        transcript = variant.annotation.transcripts.first()
        self.assertEqual(transcript.annotation, "new annotation")
        self.assertEqual(transcript.annotation_impact, "new impact")
        self.assertEqual(transcript.gene, "new gene")
        self.assertEqual(transcript.protein_impact, "new protein impact")
        self.assertEqual(transcript.transcript_id, "new transcript")
        self.assertEqual(transcript.canonical, False)

        # Delete the transcripts.
        patch_data = {"annotation": {"transcripts": []}}
        variant_count = Variant.objects.count()
        request = APIRequestFactory().patch("/variant", patch_data, format="json")
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=variant.pk)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(Variant.objects.count(), variant_count)
        variant.refresh_from_db()
        self.assertEqual(variant.annotation.type, "INDEL")
        self.assertEqual(variant.annotation.clinical_diagnosis, "new diagnosis")
        self.assertEqual(variant.annotation.clinical_significance, "new significance")
        self.assertEqual(variant.annotation.dbsnp_id, "new dbsnp")
        self.assertEqual(variant.annotation.clinvar_id, "new clinvar")
        self.assertEqual(variant.annotation.transcripts.count(), 0)

        # Delete the annotation.
        patch_data = {"annotation": None}
        variant_count = Variant.objects.count()
        request = APIRequestFactory().patch("/variant", patch_data, format="json")
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=variant.pk)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(Variant.objects.count(), variant_count)
        variant.refresh_from_db()
        with self.assertRaises(VariantAnnotation.DoesNotExist):
            variant.annotation

        # Create new annotation.
        patch_data = {
            "annotation": {
                "type": "INDEL",
                "clinical_diagnosis": "new diagnosis",
                "clinical_significance": "new significance",
                "dbsnp_id": "new dbsnp",
                "clinvar_id": "new clinvar",
                "transcripts": [
                    {
                        "annotation": "new annotation",
                        "annotation_impact": "new impact",
                        "gene": "new gene",
                        "protein_impact": "new protein impact",
                        "transcript_id": "new transcript",
                        "canonical": False,
                    }
                ],
            }
        }
        variant_count = Variant.objects.count()
        request = APIRequestFactory().patch("/variant", patch_data, format="json")
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=variant.pk)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(Variant.objects.count(), variant_count)
        variant.refresh_from_db()
        self.assertEqual(variant.annotation.type, "INDEL")
        self.assertEqual(variant.annotation.clinical_diagnosis, "new diagnosis")
        self.assertEqual(variant.annotation.clinical_significance, "new significance")
        self.assertEqual(variant.annotation.dbsnp_id, "new dbsnp")
        self.assertEqual(variant.annotation.clinvar_id, "new clinvar")
        self.assertEqual(variant.annotation.transcripts.count(), 1)
        transcript = variant.annotation.transcripts.first()
        self.assertEqual(transcript.annotation, "new annotation")
        self.assertEqual(transcript.annotation_impact, "new impact")
        self.assertEqual(transcript.gene, "new gene")
        self.assertEqual(transcript.protein_impact, "new protein impact")
        self.assertEqual(transcript.transcript_id, "new transcript")
        self.assertEqual(transcript.canonical, False)

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
        self.sample.set_permission(Permission.EDIT, self.contributor)
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
        force_authenticate(request, self.contributor)
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by depth.
        request = APIRequestFactory().get(
            "/variant", {"variant_calls__depth__isnull": "True"}
        )
        force_authenticate(request, self.contributor)
        expected = VariantSerializer([], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)
        request = APIRequestFactory().get("/variant", {"variant_calls__depth__gt": 10})
        force_authenticate(request, self.contributor)
        expected = VariantSerializer(self.variants[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by sample slug.
        request = APIRequestFactory().get(
            "/variant", {"variant_calls__sample__slug": "non_existing"}
        )
        expected = []
        force_authenticate(request, self.contributor)
        response = self.view(request)
        self.assertContains(response, expected, status_code=status.HTTP_200_OK)

        # Filter by sample slug.
        request = APIRequestFactory().get(
            "/variant", {"variant_calls__sample__slug": self.sample.slug}
        )
        force_authenticate(request, self.contributor)
        response = self.view(request)
        expected = VariantSerializer(self.variants, many=True).data
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by sample id.
        request = APIRequestFactory().get("/variant", {"variant_calls__sample": -1})
        force_authenticate(request, self.contributor)
        response = self.view(request)
        self.assertContains(
            response,
            "Select a valid choice. That choice is not one of the available choices.",
            status_code=status.HTTP_400_BAD_REQUEST,
        )

        # Filter by sample.
        request = APIRequestFactory().get(
            "/variant", {"variant_calls__sample__in": [self.sample.pk]}
        )
        force_authenticate(request, self.contributor)
        response = self.view(request)
        expected = VariantSerializer(self.variants, many=True).data
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        sample = Sample.objects.create(contributor=self.contributor)
        sample.set_permission(Permission.EDIT, self.contributor)
        self.calls[0].sample = sample
        self.calls[0].save()
        request = APIRequestFactory().get(
            "/variant", {"variant_calls__sample": self.sample.pk}
        )
        force_authenticate(request, self.contributor)
        response = self.view(request)
        expected = VariantSerializer(self.variants[1:], many=True).data
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        request = APIRequestFactory().get(
            "/variant", {"variant_calls__sample": sample.pk}
        )
        force_authenticate(request, self.contributor)
        response = self.view(request)
        expected = VariantSerializer(self.variants[:1], many=True).data
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by variant call
        request = APIRequestFactory().get(
            "/variant", {"variant_calls": self.calls[0].pk}
        )
        force_authenticate(request, self.contributor)
        response = self.view(request)
        expected = VariantSerializer(self.variants[:1], many=True).data
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by experiment
        request = APIRequestFactory().get(
            "/variant", {"variant_calls__experiment": self.experiments[0].pk}
        )
        force_authenticate(request, self.contributor)
        response = self.view(request)
        expected = VariantSerializer(self.variants[:1], many=True).data
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        self.sample.set_permission(Permission.NONE, self.contributor)

    def test_filter_sample_by_variant(self):
        client = APIClient()
        path = reverse("resolwebio-api:entity-list")

        # Anonymous user, no perms.
        response = response = client.get(
            path, {"variant_calls__variant": self.variants[0].pk}
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, [])

        # Contributor, no perms.
        client.force_authenticate(user=self.contributor)
        response = response = client.get(
            path, {"variant_calls__variant": self.variants[0].pk}
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.data, [])

        # Permission granted.
        self.sample.set_permission(Permission.VIEW, self.contributor)
        response = response = client.get(
            path, {"variant_calls__variant": self.variants[0].pk}
        )
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.data), 1)
        self.assertEqual(response.data[0]["id"], self.sample.pk)

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
            transcript_id="f1",
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
        self.view = VariantCallViewSet.as_view(
            {
                "get": "list",
                "post": "create",
                "delete": "destroy",
                "patch": "partial_update",
            }
        )
        return super().setUp()

    def test_create_api(self):
        """Test the Variant creation using the API endpoint."""
        post_data = {
            "sample": self.sample.pk,
            "variant": self.variants[0].pk,
            "quality": 0.1,
            "depth_norm_quality": 0.2,
            "alternative_allele_depth": 1,
            "depth": 42,
            "filter": "filter created",
            "genotype": "2",
            "genotype_quality": 2,
            "experiment": self.experiments[0].pk,
        }
        variant_calls_count = VariantCall.objects.count()

        # Unauthenticated request must be denied.
        request = APIRequestFactory().post("/variantcalls", post_data, format="json")
        response = self.view(request)
        self.assertContains(
            response,
            "You do not have permission to create variant calls for this sample.",
            status_code=status.HTTP_403_FORBIDDEN,
        )

        # Request from user without permission must be denied.
        self.sample.set_permission(Permission.VIEW, self.contributor)
        request = APIRequestFactory().post("/variantcalls", post_data, format="json")
        force_authenticate(request, self.contributor)
        response = self.view(request)
        self.assertContains(
            response,
            "You do not have permission to create variant calls for this sample.",
            status_code=status.HTTP_403_FORBIDDEN,
        )

        # Request from user with permission must be granted.
        self.sample.set_permission(Permission.EDIT, self.contributor)
        request = APIRequestFactory().post("/variantcalls", post_data, format="json")
        force_authenticate(request, self.contributor)
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(VariantCall.objects.count(), variant_calls_count + 1)
        variant_call = VariantCall.objects.last()
        self.assertEqual(variant_call.sample, self.sample)
        self.assertEqual(variant_call.variant, self.variants[0])
        self.assertEqual(variant_call.quality, 0.1)
        self.assertEqual(variant_call.depth_norm_quality, 0.2)
        self.assertEqual(variant_call.alternative_allele_depth, 1)
        self.assertEqual(variant_call.depth, 42)
        self.assertEqual(variant_call.filter, "filter created")
        self.assertEqual(variant_call.genotype, "2")
        self.assertEqual(variant_call.genotype_quality, 2)
        self.assertEqual(variant_call.experiment, self.experiments[0])

    def test_delete_api(self):
        """Test staff user can delete variant calls."""
        variant_call_count = VariantCall.objects.count()
        variant_call_to_delete = VariantCall.objects.first()
        # Unauthenticated request must be denied.
        request = APIRequestFactory().delete("/variantcalls")
        response = self.view(request, pk=variant_call_to_delete.pk)
        self.assertContains(
            response,
            "No VariantCall matches the given query.",
            status_code=status.HTTP_404_NOT_FOUND,
        )

        # Request from user without permission must be denied.
        variant_call_to_delete.sample.set_permission(Permission.VIEW, self.contributor)
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=variant_call_to_delete.pk)
        self.assertContains(
            response,
            "You do not have permission to perform this action.",
            status_code=status.HTTP_403_FORBIDDEN,
        )
        variant_call_to_delete.refresh_from_db()

        # Request from user with permission must be granted.
        variant_call_to_delete.sample.set_permission(Permission.EDIT, self.contributor)
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=variant_call_to_delete.pk)
        self.assertEqual(response.status_code, status.HTTP_204_NO_CONTENT)
        self.assertEqual(VariantCall.objects.count(), variant_call_count - 1)
        with self.assertRaises(VariantCall.DoesNotExist):
            variant_call_to_delete.refresh_from_db()

    def test_patch_api(self):
        """Test updating variant calls objects over api."""
        # Unauthenticated request must be denied.
        variant_call = VariantCall.objects.first()
        request = APIRequestFactory().patch("/variantcalls")
        response = self.view(request, pk=variant_call.pk)
        self.assertContains(
            response,
            "No VariantCall matches the given query.",
            status_code=status.HTTP_404_NOT_FOUND,
        )

        # Request from user without permission must be denied.
        variant_call.sample.set_permission(Permission.VIEW, self.contributor)
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=variant_call.pk)
        self.assertContains(
            response,
            "You do not have permission to perform this action.",
            status_code=status.HTTP_403_FORBIDDEN,
        )

        # Request from user with permission must be granted.
        variant_call.sample.set_permission(Permission.EDIT, self.contributor)
        patch_data: dict = {
            "quality": 42,
            "depth_norm_quality": 42,
            "alternative_allele_depth": 420,
            "depth": 420,
            "filter": "updated filter",
            "genotype": "UUU",
            # Values below must be ignored.
            "sample": 42,
            "variant": 42,
            "data": 42,
            "experiment": 42,
        }
        request = APIRequestFactory().patch("/variantcalls", patch_data, format="json")
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=variant_call.pk)
        variant_call.refresh_from_db()
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(variant_call.quality, 42)
        self.assertEqual(variant_call.depth_norm_quality, 42)
        self.assertEqual(variant_call.alternative_allele_depth, 420)
        self.assertEqual(variant_call.depth, 420)
        self.assertEqual(variant_call.filter, "updated filter")
        self.assertEqual(variant_call.genotype, "UUU")

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

        # Filter by sample.
        request = APIRequestFactory().get("/variantcall", {"sample": self.sample.id})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Filter by data.
        proc = Process.objects.create(
            type="data:test:process",
            slug="test-process",
            version="1.0.0",
            contributor=self.contributor,
        )

        data = Data.objects.create(contributor=self.contributor, process=proc)
        self.calls[0].data = data
        self.calls[0].save(update_fields=["data"])
        data.set_permission(Permission.NONE, self.contributor)
        request = APIRequestFactory().get("/variantcall", {"data": data.id})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        expected = []
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        data.set_permission(Permission.VIEW, self.contributor)
        request = APIRequestFactory().get("/variantcall", {"data": data.id})
        force_authenticate(request, self.contributor)
        expected = VariantCallSerializer(self.calls[:1], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        data.delete()
        proc.delete()
        self.calls[0].data = None
        self.calls[0].save()

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
            "/variantcall", {"variant": self.variants[0].id}
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
            "/variantcall", {"experiment": self.experiments[0].id}
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
        self.view = VariantExperimentViewSet.as_view(
            {
                "get": "list",
                "post": "create",
                "delete": "destroy",
                "patch": "partial_update",
            }
        )
        return super().setUp()

    def test_create_api(self):
        """Test the VariantExperiment creation using the API endpoint."""
        post_data = {"variant_data_source": "source 1"}

        experiment_count = VariantExperiment.objects.count()

        # Unauthenticated request must be denied.
        request = APIRequestFactory().post(
            "/variantexperiment", post_data, format="json"
        )
        response = self.view(request)
        self.assertContains(
            response,
            "Authentication credentials were not provided.",
            status_code=status.HTTP_403_FORBIDDEN,
        )

        # Request from non-staff user must be denied.
        force_authenticate(request, self.contributor)
        response = self.view(request)
        self.assertContains(
            response,
            "You do not have permission to perform this action.",
            status_code=status.HTTP_403_FORBIDDEN,
        )

        # Request from staff users must be granted.
        self.contributor.is_staff = True
        self.contributor.save(update_fields=["is_staff"])
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_201_CREATED)
        self.assertEqual(VariantExperiment.objects.count(), experiment_count + 1)
        experiment = VariantExperiment.objects.last()
        self.assertEqual(experiment.variant_data_source, "source 1")
        self.assertEqual(experiment.contributor, self.contributor)

    def test_delete_api(self):
        """Test staff user can delete variant experiments."""
        experiment_count = VariantExperiment.objects.count()
        experiment_to_delete = VariantExperiment.objects.first()
        # Unauthenticated request must be denied.
        request = APIRequestFactory().delete("/variantexperiment")
        response = self.view(request, pk=experiment_to_delete.pk)
        self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)
        self.assertEqual(
            response.data, {"detail": "Authentication credentials were not provided."}
        )

        # Request from non-staff user must be denied.
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=experiment_to_delete.pk)
        self.assertContains(
            response,
            "You do not have permission to perform this action.",
            status_code=status.HTTP_403_FORBIDDEN,
        )

        # Request from staff users must be granted.
        self.contributor.is_staff = True
        self.contributor.save(update_fields=["is_staff"])
        response = self.view(request, pk=experiment_to_delete.pk)
        self.assertEqual(response.status_code, status.HTTP_204_NO_CONTENT)
        self.assertEqual(VariantExperiment.objects.count(), experiment_count - 1)
        with self.assertRaises(VariantExperiment.DoesNotExist):
            experiment_to_delete.refresh_from_db()
        self.contributor.is_staff = False
        self.contributor.save(update_fields=["is_staff"])

    def test_patch_api(self):
        """Test updating variant experiments over api."""
        # Unauthenticated request must be denied.
        experiment = VariantExperiment.objects.first()
        request = APIRequestFactory().patch("/variantexperiment")
        response = self.view(request, pk=experiment.pk)
        self.assertEqual(response.status_code, status.HTTP_403_FORBIDDEN)
        self.assertEqual(
            response.data, {"detail": "Authentication credentials were not provided."}
        )

        # Request from non-staff user must be denied.
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=experiment.pk)
        self.assertContains(
            response,
            "You do not have permission to perform this action.",
            status_code=status.HTTP_403_FORBIDDEN,
        )

        # Request from staff users must be granted.
        self.contributor.is_staff = True
        self.contributor.save(update_fields=["is_staff"])

        patch_data: dict = {"variant_data_source": "new source", "contributor": 42}
        request = APIRequestFactory().patch(
            "/variantexperiment", patch_data, format="json"
        )
        force_authenticate(request, self.contributor)
        response = self.view(request, pk=experiment.pk)
        experiment.refresh_from_db()
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(experiment.variant_data_source, "new source")

    def test_distinct(self):
        """Test only distynct experiments are returned."""
        # Create another VariantCall that points to the first sample and experiment.
        self.sample.set_permission(Permission.VIEW, self.contributor)
        VariantCall.objects.create(
            sample=self.sample,
            variant=self.variants[0],
            quality=0.7,
            depth_norm_quality=0.7,
            alternative_allele_depth=1,
            depth=15,
            filter="filter 1",
            genotype="1",
            genotype_quality=1,
            experiment=self.experiments[0],
        )
        # Filter by sample with permissions. Make sure no duplicates are returned.
        request = APIRequestFactory().get(
            "/variantexperiment", {"variant_calls__sample": self.sample.id}
        )
        force_authenticate(request, self.contributor)
        response = self.view(request)
        expected = VariantExperimentSerializer(self.experiments, many=True).data
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

    def test_filter(self):
        # Set the permission on sample.
        self.sample.set_permission(Permission.EDIT, self.contributor)

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

        # Filter by sample, no permissions on sample.
        request = APIRequestFactory().get(
            "/variantexperiment", {"variant_calls__sample": self.sample.id}
        )
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, [])

        # Filter by sample with permissions.
        request = APIRequestFactory().get(
            "/variantexperiment", {"variant_calls__sample": self.sample.id}
        )
        force_authenticate(request, self.contributor)
        expected = VariantExperimentSerializer(self.experiments, many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # No permissions to Sample object.
        sample = Sample.objects.create(contributor=self.contributor, name="New sample")
        self.calls[1].sample = sample
        self.calls[1].save()
        request = APIRequestFactory().get(
            "/variantexperiment", {"variant_calls__sample": sample.id}
        )
        force_authenticate(request, self.contributor)
        expected = []
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, [])

        # Add permissions to Sample object.
        sample.set_permission(Permission.EDIT, self.contributor)
        request = APIRequestFactory().get(
            "/variantexperiment", {"variant_calls__sample": sample.id}
        )
        force_authenticate(request, self.contributor)
        expected = VariantExperimentSerializer(self.experiments[1:], many=True).data
        response = self.view(request)
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertCountEqual(response.data, expected)

        # Restore the test data.
        self.calls[1].sample = self.sample
        self.calls[1].save()
        sample.delete()
        self.sample.set_permission(Permission.NONE, self.contributor)

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
