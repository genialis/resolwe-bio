""".. Ignore pydocstyle D400.

====================
Variants Serializers
====================

"""

from rest_framework import serializers

from resolwe.flow.serializers.fields import DictRelatedField
from resolwe.rest.serializers import SelectiveFieldMixin

from resolwe_bio.variants.models import (
    Variant,
    VariantAnnotation,
    VariantAnnotationTranscript,
    VariantCall,
    VariantExperiment,
)


class VariantTranscriptSerializer(SelectiveFieldMixin, serializers.ModelSerializer):
    """Serializer for VariantAnnotationTranscript objects."""

    class Meta:
        """Serializer configuration."""

        model = VariantAnnotationTranscript
        fields = [
            "id",
            "variant_annotation_id",
            "annotation",
            "annotation_impact",
            "gene",
            "protein_impact",
            "transcript_id",
            "canonical",
        ]


class VariantAnnotationSerializer(SelectiveFieldMixin, serializers.ModelSerializer):
    """Serializer for VariantAnnotation objects."""

    transcripts = DictRelatedField(
        queryset=VariantAnnotationTranscript.objects.all(),
        serializer=VariantTranscriptSerializer,
        allow_null=True,
        required=False,
        many=True,
    )

    class Meta:
        """Serializer configuration."""

        model = VariantAnnotation
        fields = [
            "id",
            "variant_id",
            "type",
            "clinical_diagnosis",
            "clinical_significance",
            "dbsnp_id",
            "clinvar_id",
            "transcripts",
        ]


class VariantSerializer(SelectiveFieldMixin, serializers.ModelSerializer):
    """Serializer for Variant objects."""

    annotation = DictRelatedField(
        queryset=VariantAnnotation.objects.all(),
        serializer=VariantAnnotationSerializer,
        allow_null=True,
        required=False,
    )

    class Meta:
        """Serializer configuration."""

        model = Variant

        read_only_fields = ("id",)
        update_protected_fields = (
            "species",
            "genome_assembly",
            "chromosome",
            "position",
            "reference",
            "alternative",
            "annotation",
        )
        fields = read_only_fields + update_protected_fields


class VariantCallSerializer(SelectiveFieldMixin, serializers.ModelSerializer):
    """Serializer for VariantCall objects."""

    class Meta:
        """Serializer configuration."""

        model = VariantCall
        fields = [
            "id",
            "sample_id",
            "variant_id",
            "experiment_id",
            "quality",
            "genotype_quality",
            "depth",
            "depth_norm_quality",
            "alternative_allele_depth",
            "filter",
            "genotype",
            "data_id",
        ]


class VariantExperimentSerializer(SelectiveFieldMixin, serializers.ModelSerializer):
    """Serializer for VariantExperiment objects."""

    class Meta:
        """Serializer configuration."""

        model = VariantExperiment
        fields = ["id", "timestamp", "contributor", "variant_data_source"]
