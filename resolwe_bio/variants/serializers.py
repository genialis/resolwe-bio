""".. Ignore pydocstyle D400.

====================
Variants Serializers
====================

"""

from django.db import transaction

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
            "variant_annotation_id",
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
        update_protected_fields = ("variant_id",)
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

    @transaction.atomic
    def create(self, validated_data):
        """Create a new VariantAnnotation instance."""
        transcripts_data = validated_data.pop("transcripts", None)
        annotation = VariantAnnotation.objects.create(**validated_data)
        if transcripts_data:
            VariantAnnotationTranscript.objects.bulk_create(
                VariantAnnotationTranscript(
                    **transcript_data, variant_annotation=annotation
                )
                for transcript_data in transcripts_data
            )
        return annotation

    @transaction.atomic
    def update(self, instance, validated_data):
        """Update the variant annotation."""
        transcripts_data = validated_data.pop("transcripts", None)
        annotation = super().update(instance, validated_data)
        if transcripts_data is not None:
            # Always create new transcripts since updating logic is tedious.
            annotation.transcripts.all().delete()
            VariantAnnotationTranscript.objects.bulk_create(
                VariantAnnotationTranscript(
                    **transcript_data, variant_annotation=annotation
                )
                for transcript_data in transcripts_data
            )
        return annotation


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

    @transaction.atomic
    def create(self, validated_data):
        """Create a new Variant instance."""
        annotation_data = validated_data.pop("annotation", None)
        variant = Variant.objects.create(**validated_data)
        if annotation_data:
            VariantAnnotationSerializer().create(
                {"variant_id": variant.id, **annotation_data}
            )
        return variant

    @transaction.atomic
    def update(self, instance, validated_data):
        """Update the Variant instance."""
        delete_annotation = (
            "annotation" in validated_data and validated_data["annotation"] is None
        )
        annotation_data = validated_data.pop("annotation", None)
        variant = super().update(instance, validated_data)
        if annotation_data:
            # Update existing annotation.
            if instance.has_annotation:
                VariantAnnotationSerializer().update(
                    instance.annotation, annotation_data
                )
            # Create new annotation.
            else:
                VariantAnnotationSerializer().create(
                    {"variant_id": variant.id, **annotation_data}
                )
        elif delete_annotation and instance.annotation:
            instance.annotation.delete()
            instance.annotation = None
        return variant


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
        print("perform create", serializer.validated_data)
        serializer.save(contributor=self.context["request"].user)
        super().perform_create(serializer)

    def create(self, validated_data):
        """Create a new VariantExperiment instance."""
        print("Creating, validated_data", validated_data)
        return super().create(validated_data)
