""".. Ignore pydocstyle D400.

=====
Views
=====

"""
from django.utils.decorators import classonlymethod

from rest_framework import viewsets, mixins, permissions, filters
from rest_framework.response import Response

from drf_haystack.filters import HaystackAutocompleteFilter
from drf_haystack.serializers import HaystackSerializerMixin
from drf_haystack.viewsets import HaystackViewSet

from .models import Feature, Mapping
from .serializers import FeatureSerializer, MappingSerializer
from .filters import MappingFilter


class FeatureSearchSerializer(HaystackSerializerMixin, FeatureSerializer):
    """Feature search serializer."""

    class Meta(FeatureSerializer.Meta):
        """Meta configuration for the feature search serializer."""

        search_fields = ('genes',)
        field_aliases = {
            'query': 'genes',
        }


class FeatureSearchViewSet(HaystackViewSet):
    """
    Endpoint used for feature search.

    Request:
     - search query

    Response:
     - a list of matching features
    """

    index_models = [Feature]
    serializer_class = FeatureSearchSerializer

    @classonlymethod
    def as_view(cls, actions=None, **initkwargs):
        """Support POST for searching against a list of genes."""
        if actions.get('get', None) == 'list':
            actions['post'] = 'list_with_post'

        return super(cls, FeatureSearchViewSet).as_view(actions, **initkwargs)

    def filter_queryset(self, queryset):
        """Support filtering by a list of genes."""
        queryset = super(FeatureSearchViewSet, self).filter_queryset(queryset)
        if 'query' in self.request.data:
            queryset = queryset.filter(genes__in=self.request.data['query'])

        return queryset

    def list_with_post(self, request):
        """Support search via a POST request in addition to GET."""
        return self.list(request)


class FeatureAutocompleteSerializer(HaystackSerializerMixin, FeatureSerializer):
    """Feature autocomplete serializer."""

    class Meta(FeatureSerializer.Meta):
        """Meta configuration for the feature autocomplete serializer."""

        search_fields = ('genes_auto',)
        field_aliases = {
            'query': 'genes_auto',
        }


class FeatureAutocompleteViewSet(HaystackViewSet):
    """Endpoint used for feature autocompletion."""

    index_models = [Feature]
    serializer_class = FeatureAutocompleteSerializer
    filter_backend = [HaystackAutocompleteFilter]


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


class MappingViewSet(mixins.ListModelMixin,
                     mixins.RetrieveModelMixin,
                     viewsets.GenericViewSet):
    """API view for :class:`Mapping` objects."""

    serializer_class = MappingSerializer
    filter_backends = [filters.DjangoFilterBackend]
    filter_class = MappingFilter
    queryset = Mapping.objects.all()

    @classonlymethod
    def as_view(cls, actions=None, **initkwargs):
        """Support POST for searching against a list of feature ids."""
        if actions.get('get', None) == 'list':
            actions['post'] = 'list_with_post'

        return super(cls, MappingViewSet).as_view(actions, **initkwargs)

    def list_with_post(self, request):
        """Support search via a POST request in addition to GET."""
        queryset = self.queryset
        filter_params = request.data
        for key, value in filter_params.items():
            if isinstance(value, list):
                filter_params[key] = ','.join(value)

        if filter_params:
            queryfilter = self.filter_class(filter_params, queryset=queryset)
            queryset = queryfilter.qs
        serializer = self.get_serializer(queryset, many=True)
        return Response(serializer.data)
