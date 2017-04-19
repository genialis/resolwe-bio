""".. Ignore pydocstyle D400.

=====
Views
=====

"""
from elasticsearch_dsl.query import Q

from rest_framework import viewsets, mixins, permissions, filters

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

    filtering_fields = ('name', 'source', 'species')
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


class FeatureAutocompleteViewSet(ElasticSearchBaseViewSet):
    """Endpoint used for feature autocompletion."""

    document_class = FeatureSearchDocument
    serializer_class = FeatureSerializer

    filtering_fields = ('source', 'species')

    def custom_filter(self, search):
        """Support autocomplete query using the 'query' attribute."""
        return search.query(
            'bool',
            should=[
                Q('match', autocomplete=self.get_query_param('query', '')),
                # Boost exact name and feature_id matches. Use the 'lower' subfield in order
                # to compare against lowercased terms.
                Q('match', **{'feature_id.lower': {'query': self.get_query_param('query', ''), 'boost': 2}}),
                Q('match', **{'name.lower': {'query': self.get_query_param('query', ''), 'boost': 2}}),
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
    filter_backends = [filters.DjangoFilterBackend]
    queryset = Feature.objects.all()

    def create(self, request):
        """Instead of failing, update existing features with a custom create."""
        try:
            feature = Feature.objects.get(
                source=request.data['source'],
                feature_id=request.data['feature_id']
            )
            self.kwargs[self.lookup_field] = feature.pk
            return super(FeatureViewSet, self).update(request)  # pylint: disable=no-member
        except (Feature.DoesNotExist, KeyError):  # pylint: disable=no-member
            return super(FeatureViewSet, self).create(request)  # pylint: disable=no-member


class MappingSearchViewSet(ElasticSearchBaseViewSet):
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

    document_class = MappingSearchDocument
    serializer_class = MappingSerializer

    filtering_fields = ('source_db', 'target_db', 'relation_type')
    ordering_fields = ('source_id',)
    ordering = 'source_id'

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
    filter_backends = [filters.DjangoFilterBackend]
    filter_class = MappingFilter
    queryset = Mapping.objects.all()

    def create(self, request):
        """Instead of failing, update existing mappings using a custom create."""
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
