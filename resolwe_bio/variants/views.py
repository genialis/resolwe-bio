""".. Ignore pydocstyle D400.

=============================
Expose Variants models on API
=============================

"""

import logging

import django_filters as filters
from rest_framework import exceptions, mixins, permissions, serializers, viewsets

from resolwe.flow.filters import OrderingFilter
from resolwe.flow.views.mixins import ResolweCreateModelMixin
from resolwe.flow.views.utils import IsStaffOrReadOnly
from resolwe.permissions.models import Permission

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


class VariantViewSet(
    mixins.ListModelMixin,
    ResolweCreateModelMixin,
    mixins.DestroyModelMixin,
    viewsets.GenericViewSet,
):
    """Variant endpoint."""

    queryset = Variant.objects.all()
    serializer_class = VariantSerializer
    filter_backends = [filters.rest_framework.DjangoFilterBackend, OrderingFilter]
    permission_classes = (IsStaffOrReadOnly,)

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


class VariantCallViewSet(
    mixins.ListModelMixin,
    ResolweCreateModelMixin,
    mixins.DestroyModelMixin,
    viewsets.GenericViewSet,
):
    """VariantCall endpoint.

    The default filter backends are used so permissions are respected.
    """

    queryset = VariantCall.objects.all()
    serializer_class = VariantCallSerializer

    filterset_class = VariantCallFilter
    ordering_fields = ("id", "quality", "depth")

    permission_classes = (permissions.IsAuthenticatedOrReadOnly,)

    def perform_create(self, serializer: serializers.BaseSerializer) -> None:
        """Check if user has VIEW permission on the sample."""
        # Check permission on sample.
        sample = serializer.validated_data["sample"]
        if not sample.has_permission(Permission.VIEW, self.request.user):
            raise exceptions.PermissionDenied()

        # Check permission on data (if given).
        if data := serializer.validated_data.get("data"):
            if not data.has_permission(Permission.VIEW, self.request.user):
                raise exceptions.PermissionDenied()

        return super().perform_create(serializer)

    def perform_destroy(self, instance: VariantCall) -> None:
        """Check if user has EDIT permission on the sample."""
        if not instance.sample.has_permission(Permission.EDIT, self.request.user):
            raise exceptions.PermissionDenied()
        if data := instance.data:
            if not data.has_permission(Permission.EDIT, self.request.user):
                raise exceptions.PermissionDenied()

        return super().perform_destroy(instance)


class VariantExperimentViewSet(
    mixins.ListModelMixin,
    ResolweCreateModelMixin,
    mixins.DestroyModelMixin,
    viewsets.GenericViewSet,
):
    """VariantExperiment endpoint."""

    queryset = VariantExperiment.objects.all()
    serializer_class = VariantExperimentSerializer
    filter_backends = [filters.rest_framework.DjangoFilterBackend, OrderingFilter]
    permission_classes = (IsStaffOrReadOnly,)

    filterset_class = VariantExperimentFilter
    ordering_fields = ("id", "timestamp", "contributor__email")
