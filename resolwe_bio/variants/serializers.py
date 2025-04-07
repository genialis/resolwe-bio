""".. Ignore pydocstyle D400.

====================
Variants Serializers
====================

"""

from resolwe.flow.serializers import ResolweBaseSerializer

from resolwe_bio.variants.models import (
    Variant,
    VariantAnnotation,
    VariantAnnotationTranscript,
    VariantCall,
    VariantExperiment,
)


class VariantTranscriptSerializer(ResolweBaseSerializer):
    """Serializer for VariantAnnotationTranscript objects."""

    class Meta:
        """Serializer configuration."""

        model = VariantAnnotationTranscript
        fields = [
            "id",
            "variant_annotation",
            "annotation",
            "annotation_impact",
            "gene",
            "protein_impact",
            "transcript_id",
            "canonical",
        ]


class VariantAnnotationSerializer(ResolweBaseSerializer):
    """Serializer for VariantAnnotation objects."""

    transcripts = VariantTranscriptSerializer(many=True, required=False)

    class Meta:
        """Serializer configuration."""

        model = VariantAnnotation
        read_only_fields = ("id",)
        update_protected_fields = ("variant",)
        fields = (
            read_only_fields
            + update_protected_fields
            + (
                "type",
                "clinical_diagnosis",
                "clinical_significance",
                "dbsnp_id",
                "clinvar_id",
                "transcripts",
            )
        )


class VariantSerializer(ResolweBaseSerializer):
    """Serializer for Variant objects."""

    annotation = VariantAnnotationSerializer(required=False, allow_null=True)

    class Meta:
        """Serializer configuration."""

        model = Variant
        read_only_fields = ("id",)
        fields = read_only_fields + (
            "species",
            "genome_assembly",
            "chromosome",
            "position",
            "reference",
            "alternative",
            "annotation",
        )


class VariantCallSerializer(ResolweBaseSerializer):
    """Serializer for VariantCall objects."""

    class Meta:
        """Serializer configuration."""

        model = VariantCall
        read_only_fields = ("id",)
        update_protected_fields = (
            "sample",
            "variant",
            "experiment",
            "data",
        )
        fields = (
            read_only_fields
            + update_protected_fields
            + (
                "quality",
                "genotype_quality",
                "depth",
                "depth_norm_quality",
                "alternative_allele_depth",
                "filter",
                "genotype",
            )
        )


class VariantExperimentSerializer(ResolweBaseSerializer):
    """Serializer for VariantExperiment objects."""

    class Meta:
        """Serializer configuration."""

        model = VariantExperiment
        read_only_fields = ("id",)
        update_protected_fields = ("timestamp", "contributor")
        fields = read_only_fields + update_protected_fields + ("variant_data_source",)

    def perform_create(self, serializer):
        """Set the contributor to the current user."""
        serializer.save(contributor=self.context["request"].user)
        super().perform_create(serializer)

    def create(self, validated_data):
        """Create a new VariantExperiment instance."""
        return super().create(validated_data)
