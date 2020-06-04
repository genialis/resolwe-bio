""".. Ignore pydocstyle D400.

=====
Views
=====

"""
from rest_framework import mixins, permissions, viewsets
from tests.backends import ResolweBioFilterBackend
from tests.pagination import LimitOffsetPostPagination

from .filters import FeatureAutoCompleteFilter, FeatureFilter, MappingFilter
from .models import Feature, Mapping
from .serializers import FeatureSerializer, MappingSerializer


class FeatureSearchViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """Endpoint used for feature search.

    Request:
     - query
     - source

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


class FeatureAutocompleteViewSet(mixins.ListModelMixin, viewsets.GenericViewSet):
    """Endpoint used for feature autocompletion."""

    queryset = Feature.objects.all()
    serializer_class = FeatureSerializer
    filter_backends = [ResolweBioFilterBackend]
    filter_class = FeatureAutoCompleteFilter
    pagination_class = LimitOffsetPostPagination

    def list_with_post(self, request):
        """Endpoint handler."""
        return self.list(request)


class FeatureViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    mixins.CreateModelMixin,
    mixins.UpdateModelMixin,
    mixins.DestroyModelMixin,
    viewsets.GenericViewSet,
):
    """API view for :class:`Feature` objects."""

    queryset = Feature.objects.all()
    serializer_class = FeatureSerializer
    permission_classes = [permissions.IsAdminUser]

    def create(self, request, *args, **kwargs):
        """Instead of failing, update existing features with a custom create."""
        try:
            feature = Feature.objects.get(
                source=request.data["source"], feature_id=request.data["feature_id"]
            )
            self.kwargs[self.lookup_field] = feature.pk
            return super(FeatureViewSet, self).update(request, *args, **kwargs)
        except (Feature.DoesNotExist, KeyError):
            return super(FeatureViewSet, self).create(request, *args, **kwargs)


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


class MappingViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    mixins.CreateModelMixin,
    mixins.UpdateModelMixin,
    mixins.DestroyModelMixin,
    viewsets.GenericViewSet,
):
    """API view for :class:`Mapping` objects."""

    serializer_class = MappingSerializer
    permission_classes = [permissions.IsAdminUser]
    filter_backends = [DjangoFilterBackend]
    filter_class = MappingFilter
    queryset = Mapping.objects.all()

    def create(self, request, *args, **kwargs):
        """Instead of failing, update existing mappings using a custom create."""
        try:
            mapping = Mapping.objects.get(
                source_db=request.data["source_db"],
                source_id=request.data["source_id"],
                source_species=request.data["source_species"],
                target_db=request.data["target_db"],
                target_id=request.data["target_id"],
                target_species=request.data["target_species"],
            )
            self.kwargs[self.lookup_field] = mapping.pk
            return super(MappingViewSet, self).update(request, *args, **kwargs)
        except (Mapping.DoesNotExist, KeyError):
            return super(MappingViewSet, self).create(request, *args, **kwargs)
