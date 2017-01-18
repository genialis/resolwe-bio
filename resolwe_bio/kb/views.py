""".. Ignore pydocstyle D400.

=====
Views
=====

"""
from django.utils.decorators import classonlymethod

from elasticsearch_dsl.query import Q

from rest_framework import viewsets, mixins, permissions, filters

from haystack.query import SQ

from drf_haystack.viewsets import HaystackViewSet

from resolwe.elastic.viewsets import ElasticSearchBaseViewSet

from .models import Feature, Mapping
from .serializers import (FeatureSerializer, FeatureAutocompleteSerializer,
                          MappingSerializer, MappingSearchSerializer)
from .filters import MappingFilter
from .elastic_indexes import FeatureSearchDocument


class FeatureSearchViewSet(ElasticSearchBaseViewSet):
    """
    Endpoint used for feature search.

    Request:
     - query
     - source

    Response:
     - a list of matching features
    """

    document_class = FeatureSearchDocument
    serializer_class = FeatureSerializer

    filtering_fields = ('name', 'source')
    ordering_fields = ('name',)
    ordering = 'name'

    def custom_filter(self, search):
        """Support general query using the 'query' attribute."""
        query = self.get_query_param('query', None)
        if query:
            if not isinstance(query, list):
                query = [query]

            search = search.filter(
                'bool',
                should=[
                    Q('terms', feature_id=query),
                    Q('terms', name=query),
                    Q('terms', aliases=query),
                ]
            )

        return search

    def filter_permissions(self, search):
        """Feature objects have no permissions."""
        return search


class FeatureAutocompleteViewSet(HaystackViewSet):
    """Endpoint used for feature autocompletion."""

    index_models = [Feature]
    serializer_class = FeatureAutocompleteSerializer

    def filter_queryset(self, queryset):
        """Construct a correct filter query."""
        query = self.request.query_params.get('query', None)

        if not query:
            return queryset.none()

        queryset = queryset.filter(SQ(name_auto=query) | SQ(aliases_auto=query))

        source = self.request.query_params.get('source', None)
        if source:
            queryset = queryset.filter(source=source)

        return queryset


class FeatureViewSet(mixins.ListModelMixin,
                     mixins.RetrieveModelMixin,
                     mixins.CreateModelMixin,
                     mixins.UpdateModelMixin,
                     mixins.DestroyModelMixin,
                     viewsets.GenericViewSet):
    """API view for :class:`Feature` objects."""

    serializer_class = FeatureSerializer
    permission_classes = [permissions.IsAdminUser]
    filter_backends = [filters.DjangoFilterBackend]
    queryset = Feature.objects.all()

    def create(self, request):
        """A custom create, which updates existing features instead of failing."""
        try:
            feature = Feature.objects.get(
                source=request.data['source'],
                feature_id=request.data['feature_id']
            )
            self.kwargs[self.lookup_field] = feature.pk
            return super(FeatureViewSet, self).update(request)  # pylint: disable=no-member
        except (Feature.DoesNotExist, KeyError):  # pylint: disable=no-member
            return super(FeatureViewSet, self).create(request)  # pylint: disable=no-member


class MappingSearchViewSet(HaystackViewSet):
    """
    Endpoint used for mapping search.

    Request:
     - source_id
     - source_db
     - target_id
     - target_db
     - relation_type

    Response:
     - a list of matching mappings
    """

    index_models = [Mapping]
    serializer_class = MappingSearchSerializer

    @classonlymethod
    def as_view(cls, actions=None, **initkwargs):
        """Support POST for searching against a list of genes."""
        if actions.get('get', None) == 'list':
            actions['post'] = 'list_with_post'

        return super(cls, MappingSearchViewSet).as_view(actions, **initkwargs)

    def filter_queryset(self, queryset):
        """Support filtering by a list of genes."""
        queryset = super(MappingSearchViewSet, self).filter_queryset(queryset)
        if 'source_id' in self.request.data:
            queryset = queryset.filter(source_id__in=self.request.data['source_id'])
        if 'source_db' in self.request.data:
            queryset = queryset.filter(source_db=self.request.data['source_db'])
        if 'target_id' in self.request.data:
            queryset = queryset.filter(target_id__in=self.request.data['target_id'])
        if 'target_db' in self.request.data:
            queryset = queryset.filter(target_db=self.request.data['target_db'])
        if 'relation_type' in self.request.data:
            queryset = queryset.filter(relation_type=self.request.data['relation_type'])

        return queryset

    def list_with_post(self, request):
        """Support search via a POST request in addition to GET."""
        return self.list(request)


class MappingViewSet(mixins.ListModelMixin,
                     mixins.RetrieveModelMixin,
                     mixins.CreateModelMixin,
                     mixins.UpdateModelMixin,
                     mixins.DestroyModelMixin,
                     viewsets.GenericViewSet):
    """API view for :class:`Mapping` objects."""

    serializer_class = MappingSerializer
    permission_classes = [permissions.IsAdminUser]
    filter_backends = [filters.DjangoFilterBackend]
    filter_class = MappingFilter
    queryset = Mapping.objects.all()

    def create(self, request):
        """A custom create, which updates existing mappings instead of failing."""
        try:
            mapping = Mapping.objects.get(
                source_db=request.data['source_db'],
                source_id=request.data['source_id'],
                target_db=request.data['target_db'],
                target_id=request.data['target_id']
            )
            self.kwargs[self.lookup_field] = mapping.pk
            return super(MappingViewSet, self).update(request)  # pylint: disable=no-member
        except (Mapping.DoesNotExist, KeyError):  # pylint: disable=no-member
            return super(MappingViewSet, self).create(request)  # pylint: disable=no-member
