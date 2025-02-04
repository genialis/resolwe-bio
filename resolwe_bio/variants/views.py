""".. Ignore pydocstyle D400.

=============================
Expose Variants models on API
=============================

"""

import logging

import django_filters as filters
from rest_framework import exceptions, mixins, viewsets

from resolwe.flow.filters import OrderingFilter
from resolwe.flow.models import Data
from resolwe.flow.models import Entity as Sample
from resolwe.flow.views.mixins import ResolweCreateModelMixin
from resolwe.flow.views.utils import IsStaffOrReadOnly
from resolwe.permissions.loader import get_permissions_class
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
    mixins.UpdateModelMixin,
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
    mixins.UpdateModelMixin,
    viewsets.GenericViewSet,
):
    """VariantCall endpoint.

    The default filter backends are used so permissions are respected.
    """

    queryset = VariantCall.objects.all()
    serializer_class = VariantCallSerializer

    filterset_class = VariantCallFilter
    ordering_fields = ("id", "quality", "depth")

    permission_classes = (get_permissions_class(),)

    def _check_create_permissions(self, sample: Sample, data: Data):
        """Check if user can create the variant call.

        :raises PermissionDenied: if user does not have permission to create the variant call.

        :raises ValidationError: if the data object does not belong to the sample object.
        """

        if data and (not data.entity_id == sample.id):
            raise exceptions.ValidationError(
                "The data object does not belong to the sample object."
            )

        if not sample.has_permission(Permission.EDIT, self.request.user):
            raise exceptions.PermissionDenied(
                "You do not have permission to create variant calls for this sample."
            )

    def perform_create(self, serializer):
        """Check if the user can create the variant call object."""
        self._check_create_permissions(
            serializer.validated_data["sample"], serializer.validated_data.get("data")
        )
        return super().perform_create(serializer)


class VariantExperimentViewSet(
    mixins.ListModelMixin,
    ResolweCreateModelMixin,
    mixins.DestroyModelMixin,
    mixins.UpdateModelMixin,
    viewsets.GenericViewSet,
):
    """VariantExperiment endpoint."""

    queryset = VariantExperiment.objects.all()
    serializer_class = VariantExperimentSerializer
    filter_backends = [filters.rest_framework.DjangoFilterBackend, OrderingFilter]

    filterset_class = VariantExperimentFilter
    ordering_fields = ("id", "timestamp", "contributor__email")
    permission_classes = (IsStaffOrReadOnly,)
