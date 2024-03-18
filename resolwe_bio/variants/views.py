""".. Ignore pydocstyle D400.

=============================
Expose Variants models on API
=============================

"""

import logging

import django_filters as filters
from rest_framework import mixins, serializers, viewsets

from resolwe.flow.filters import (
    DATE_LOOKUPS,
    NUMBER_LOOKUPS,
    TEXT_LOOKUPS,
    CheckQueryParamsMixin,
    OrderingFilter,
)
from resolwe.rest.serializers import SelectiveFieldMixin

from .models import Variant, VariantAnnotation, VariantCall, VariantExperiment

logger = logging.getLogger(__name__)


class VariantSerializer(SelectiveFieldMixin, serializers.ModelSerializer):
    """Serializer for Variant objects."""

    class Meta:
        """Serializer configuration."""

        model = Variant
        fields = [
            "id",
            "species",
            "genome_assembly",
            "chromosome",
            "position",
            "reference",
            "alternative",
            "annotation",
        ]


class VariantFilter(CheckQueryParamsMixin, filters.FilterSet):
    """Filter the Variant objects endpoint."""

    class Meta:
        """Filter configuration."""

        model = Variant
        fields = {
            "id": NUMBER_LOOKUPS,
            "species": TEXT_LOOKUPS,
            "genome_assembly": TEXT_LOOKUPS,
            "chromosome": TEXT_LOOKUPS,
            "position": NUMBER_LOOKUPS,
            "reference": TEXT_LOOKUPS,
            "alternative": TEXT_LOOKUPS,
            "annotation__type": TEXT_LOOKUPS,
            "annotation__annotation": TEXT_LOOKUPS,
            "annotation__annotation_impact": TEXT_LOOKUPS,
            "annotation__gene": TEXT_LOOKUPS,
            "annotation__protein_impact": TEXT_LOOKUPS,
            "annotation__clinical_diagnosis": TEXT_LOOKUPS,
            "annotation__clinical_significance": TEXT_LOOKUPS,
            "dbsnp_id": TEXT_LOOKUPS,
            "clinical_var_id": TEXT_LOOKUPS,
            "variant_calls__quality": NUMBER_LOOKUPS,
            "variant_calls__depth": NUMBER_LOOKUPS,
        }


class VariantViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """Variant endpoint."""

    queryset = Variant.objects.all()
    serializer_class = VariantSerializer
    filter_backends = [OrderingFilter]

    filterset_class = VariantFilter
    ordering_fields = ("species", "genome_assembly", "position", "chromosome")


class VariantAnnotationSerializer(SelectiveFieldMixin, serializers.ModelSerializer):
    """Serializer for VariantAnnotation objects."""

    class Meta:
        """Serializer configuration."""

        model = VariantAnnotation
        fields = [
            "id",
            "variant_id",
            "type",
            "annotation",
            "annotation_impact",
            "gene",
            "protein_impact",
            "feature_id",
            "clinical_diagnosis",
            "clinical_significance",
            "dbsnp_id",
            "clinical_var_id",
            "data_id",
        ]


class VariantAnnotationFilter(CheckQueryParamsMixin, filters.FilterSet):
    """Filter the VariantAnnotation objects endpoint."""

    class Meta:
        """Filter configuration."""

        model = VariantAnnotation
        fields = {
            "id": NUMBER_LOOKUPS,
            "type": TEXT_LOOKUPS,
            "annotation": TEXT_LOOKUPS,
            "annotation_impact": TEXT_LOOKUPS,
            "gene": TEXT_LOOKUPS,
            "protein_impact": TEXT_LOOKUPS,
            "clinical_diagnosis": TEXT_LOOKUPS,
            "clinical_significance": TEXT_LOOKUPS,
            "dbsnp_id": TEXT_LOOKUPS,
            "clinical_var_id": TEXT_LOOKUPS,
        }


class VariantAnnotationViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """VariantAnnotation endpoint."""

    queryset = VariantAnnotation.objects.all()
    serializer_class = VariantAnnotationSerializer
    filter_backends = [OrderingFilter]

    filterset_class = VariantAnnotationFilter
    ordering_fields = ("gene", "protein_impact", "annotation", "clinical_significance")


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
            "depth",
            "filter",
            "genotype",
            "data_id",
        ]


class VariantCallFilter(CheckQueryParamsMixin, filters.FilterSet):
    """Filter the VariantCall objects endpoint."""

    class Meta:
        """Filter configuration."""

        model = VariantCall
        fields = {
            "id": NUMBER_LOOKUPS,
            "sample__slug": TEXT_LOOKUPS,
            "sample__id": NUMBER_LOOKUPS,
            "data__slug": TEXT_LOOKUPS,
            "data__id": NUMBER_LOOKUPS,
            "variant__id": NUMBER_LOOKUPS,
            "variant__species": TEXT_LOOKUPS,
            "variant__genome_assembly": TEXT_LOOKUPS,
            "variant__chromosome": TEXT_LOOKUPS,
            "variant__position": NUMBER_LOOKUPS,
            "variant__reference": TEXT_LOOKUPS,
            "variant__alternative": TEXT_LOOKUPS,
            "variant__annotation": TEXT_LOOKUPS,
            "experiment__id": NUMBER_LOOKUPS,
        }


class VariantCallViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """VariantCall endpoint."""

    queryset = VariantCall.objects.all()
    serializer_class = VariantCallSerializer
    filter_backends = [OrderingFilter]

    filterset_class = VariantCallFilter
    ordering_fields = ("id", "quality", "depth")


class VariantExperimentSerializer(SelectiveFieldMixin, serializers.ModelSerializer):
    """Serializer for VariantExperiment objects."""

    class Meta:
        """Serializer configuration."""

        model = VariantExperiment
        fields = ["id", "date", "contributor", "variant_data_source"]


class VariantExperimentFilter(CheckQueryParamsMixin, filters.FilterSet):
    """Filter the VariantExperiment objects endpoint."""

    class Meta:
        """Filter configuration."""

        model = VariantExperiment
        fields = {
            "id": NUMBER_LOOKUPS,
            "contributor__username": TEXT_LOOKUPS,
            "date": DATE_LOOKUPS,
            "variant_data_source": TEXT_LOOKUPS,
        }


class VariantExperimentViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """VariantExperiment endpoint."""

    queryset = VariantExperiment.objects.all()
    serializer_class = VariantExperimentSerializer
    filter_backends = [OrderingFilter]

    filterset_class = VariantExperiment
    ordering_fields = ("id", "date", "contributor")
