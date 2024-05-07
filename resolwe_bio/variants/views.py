""".. Ignore pydocstyle D400.

=============================
Expose Variants models on API
=============================

"""

import logging

import django_filters as filters
from rest_framework import mixins, viewsets

from resolwe.flow.filters import OrderingFilter

from resolwe_bio.variants.filters import (
    VariantAnnotationFilter,
    VariantCallFilter,
    VariantExperimentFilter,
    VariantFilter,
)
from resolwe_bio.variants.serializers import (
    VariantAnnotationSerializer,
    VariantCallSerializer,
    VariantExperimentSerializer,
    VariantSerializer,
)

from .models import Variant, VariantAnnotation, VariantCall, VariantExperiment

logger = logging.getLogger(__name__)


class VariantViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """Variant endpoint."""

    queryset = Variant.objects.all()
    serializer_class = VariantSerializer
    filter_backends = [filters.rest_framework.DjangoFilterBackend, OrderingFilter]

    filterset_class = VariantFilter
    ordering_fields = ("species", "genome_assembly", "position", "chromosome")


class VariantAnnotationViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """VariantAnnotation endpoint."""

    queryset = VariantAnnotation.objects.all()
    serializer_class = VariantAnnotationSerializer
    filter_backends = [filters.rest_framework.DjangoFilterBackend, OrderingFilter]

    filterset_class = VariantAnnotationFilter
    ordering_fields = (
        "transcripts__gene",
        "transcripts__protein_impact",
        "transcripts__annotation",
        "clinical_significance",
    )


class VariantCallViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """VariantCall endpoint.

    The default filter backends are used so permissions are respected.
    """

    queryset = VariantCall.objects.all()
    serializer_class = VariantCallSerializer

    filterset_class = VariantCallFilter
    ordering_fields = ("id", "quality", "depth")


class VariantExperimentViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """VariantExperiment endpoint."""

    queryset = VariantExperiment.objects.all()
    serializer_class = VariantExperimentSerializer
    filter_backends = [filters.rest_framework.DjangoFilterBackend, OrderingFilter]

    filterset_class = VariantExperimentFilter
    ordering_fields = ("id", "timestamp", "contributor__email")
