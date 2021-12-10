""".. Ignore pydocstyle D400.

=====
Views
=====

"""
from rest_framework import mixins, viewsets

from .backends import ResolweBioFilterBackend
from .filters import FeatureFilter, MappingFilter
from .models import Feature, Mapping
from .pagination import LimitOffsetPostPagination
from .serializers import FeatureSerializer, MappingSerializer


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
    filter_backends = [ResolweBioFilterBackend]
    filter_class = FeatureFilter
    pagination_class = LimitOffsetPostPagination

    ordering_fields = ("name",)
    ordering = "name"

    def list_with_post(self, request):
        """Endpoint handler."""
        return self.list(request)


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
    filter_backends = [ResolweBioFilterBackend]
    filter_class = MappingFilter
    pagination_class = LimitOffsetPostPagination

    ordering_fields = ("source_id",)
    ordering = "source_id"

    def list_with_post(self, request):
        """Endpoint handler."""
        return self.list(request)
