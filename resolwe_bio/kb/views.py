""".. Ignore pydocstyle D400.

=====
Views
=====

"""
from elasticsearch_dsl.query import Q

from rest_framework import viewsets, mixins, permissions
from django_filters.rest_framework.backends import DjangoFilterBackend

from resolwe.elastic.viewsets import ElasticSearchBaseViewSet

from .models import Feature, Mapping
from .serializers import FeatureSerializer, MappingSerializer
from .filters import MappingFilter

from .elastic_indexes import FeatureSearchDocument, MappingSearchDocument


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

    filtering_fields = ('name', 'source', 'species', 'feature_id', 'type')
    ordering_fields = ('name',)
    ordering = 'name'

    def get_always_allowed_arguments(self):
        """Return query arguments which are always allowed."""
        return super().get_always_allowed_arguments() + [
            'query',
        ]

    def custom_filter_feature_id(self, value, search):
        """Support exact feature_id queries."""
        if not isinstance(value, list):
            value = [value]

        return search.filter('terms', feature_id=value)

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


class FeatureAutocompleteViewSet(ElasticSearchBaseViewSet):
    """Endpoint used for feature autocompletion."""

    document_class = FeatureSearchDocument
    serializer_class = FeatureSerializer

    filtering_fields = ('source', 'species', 'type')

    def get_always_allowed_arguments(self):
        """Return query arguments which are always allowed."""
        return super().get_always_allowed_arguments() + [
            'query',
        ]

    def custom_filter(self, search):
        """Support autocomplete query using the 'query' attribute."""
        query = self.get_query_param('query', '')
        return search.query(
            'bool',
            should=[
                # Exact matches of name and feature_id.
                Q('match', **{'feature_id.lower': {'query': query, 'boost': 10}}),
                Q('match', **{'name.lower': {'query': query, 'boost': 10}}),
                # Partial matches of name and feature_id.
                Q('match', **{'feature_id.ngrams': {'query': query, 'boost': 5}}),
                Q('match', **{'name.ngrams': {'query': query, 'boost': 5}}),
                # Aliases.
                Q('match', **{'aliases.ngrams': {'query': query}}),
            ]
        )

    def filter_permissions(self, search):
        """Feature objects have no permissions."""
        return search


class FeatureViewSet(mixins.ListModelMixin,
                     mixins.RetrieveModelMixin,
                     mixins.CreateModelMixin,
                     mixins.UpdateModelMixin,
                     mixins.DestroyModelMixin,
                     viewsets.GenericViewSet):
    """API view for :class:`Feature` objects."""

    serializer_class = FeatureSerializer
    permission_classes = [permissions.IsAdminUser]
    filter_backends = [DjangoFilterBackend]
    queryset = Feature.objects.all()

    def create(self, request, *args, **kwargs):
        """Instead of failing, update existing features with a custom create."""
        try:
            feature = Feature.objects.get(
                source=request.data['source'],
                feature_id=request.data['feature_id']
            )
            self.kwargs[self.lookup_field] = feature.pk
            return super(FeatureViewSet, self).update(request, *args, **kwargs)  # pylint: disable=no-member
        except (Feature.DoesNotExist, KeyError):
            return super(FeatureViewSet, self).create(request, *args, **kwargs)  # pylint: disable=no-member


class MappingSearchViewSet(ElasticSearchBaseViewSet):
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

    document_class = MappingSearchDocument
    serializer_class = MappingSerializer

    filtering_fields = ('source_db', 'source_species', 'target_db', 'target_species', 'relation_type')
    ordering_fields = ('source_id',)
    ordering = 'source_id'

    def get_always_allowed_arguments(self):
        """Return query arguments which are always allowed."""
        return super().get_always_allowed_arguments() + [
            'source_id',
            'target_id',
        ]

    def custom_filter(self, search):
        """Support correct searching by ``source_id`` and ``target_id``."""
        for field in ('source_id', 'target_id'):
            query = self.get_query_param(field, None)
            if not query:
                continue
            if not isinstance(query, list):
                query = [query]

            search = search.filter(
                'bool',
                should=[
                    Q('terms', **{field: query}),
                ]
            )

        return search

    def filter_permissions(self, search):
        """Filter permissions since Mapping objects have no permissions."""
        return search


class MappingViewSet(mixins.ListModelMixin,
                     mixins.RetrieveModelMixin,
                     mixins.CreateModelMixin,
                     mixins.UpdateModelMixin,
                     mixins.DestroyModelMixin,
                     viewsets.GenericViewSet):
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
                source_db=request.data['source_db'],
                source_id=request.data['source_id'],
                source_species=request.data['source_species'],
                target_db=request.data['target_db'],
                target_id=request.data['target_id'],
                target_species=request.data['target_species'],
            )
            self.kwargs[self.lookup_field] = mapping.pk
            return super(MappingViewSet, self).update(request, *args, **kwargs)  # pylint: disable=no-member
        except (Mapping.DoesNotExist, KeyError):
            return super(MappingViewSet, self).create(request, *args, **kwargs)  # pylint: disable=no-member
