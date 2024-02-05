""".. Ignore pydocstyle D400.

=====
Views
=====

"""

import logging

from django.db.models import Q
from rest_framework import mixins, viewsets
from rest_framework.decorators import action
from rest_framework.response import Response

from resolwe.flow.filters import OrderingFilter

from .backends import ResolweBioFilterBackend
from .filters import FeatureFilter, MappingFilter
from .models import Feature, Mapping
from .pagination import LimitOffsetPostPagination
from .serializers import FeatureSerializer, MappingSerializer

logger = logging.getLogger(__name__)


class FeatureViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """Endpoint used for feature search.

    Request:
     - query
     - source
     - species
     - type
     - feature_id

    Response:
     - a list of matching features
    """

    queryset = Feature.objects.all()
    serializer_class = FeatureSerializer
    filter_backends = [ResolweBioFilterBackend, OrderingFilter]
    filterset_class = FeatureFilter
    pagination_class = LimitOffsetPostPagination

    ordering_fields = ("name",)

    def list_with_post(self, request):
        """Endpoint handler."""
        return self.list(request)

    @action(methods=["post"], detail=False)
    def paste(self, request):
        """Endpoint used for exactly matching pasted genes."""
        queryset = self.filter_queryset(self.get_queryset())

        pasted = request.data["pasted"]
        queryset = queryset.filter(Q(name__in=pasted) | Q(feature_id__in=pasted))

        logger.info("Pasted genes: {}".format(", ".join(pasted)))

        page = self.paginate_queryset(queryset)
        if page is not None:
            serializer = self.get_serializer(page, many=True)
            return self.get_paginated_response(serializer.data)

        serializer = self.get_serializer(queryset, many=True)
        return Response(serializer.data)


class MappingSearchViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """
    Endpoint used for mapping search.

    Request:
     - source_id
     - source_db
     - source_species
     - target_id
     - target_db
     - target_species
     - relation_type

    Response:
     - a list of matching mappings
    """

    queryset = Mapping.objects.all()
    serializer_class = MappingSerializer
    filter_backends = [ResolweBioFilterBackend, OrderingFilter]
    filterset_class = MappingFilter
    pagination_class = LimitOffsetPostPagination

    ordering_fields = ("source_id",)
    ordering = "source_id"

    def list_with_post(self, request):
        """Endpoint handler."""
        return self.list(request)
