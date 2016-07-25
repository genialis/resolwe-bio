from __future__ import absolute_import, division, print_function, unicode_literals

from django.db.models import Max

from rest_framework import exceptions, status
from rest_framework.decorators import detail_route
from rest_framework.response import Response

from resolwe.flow.models import Collection
from resolwe.flow.views import CollectionViewSet
from .filters import SampleFilter
from .models import Sample
from .serializers import PresampleSerializer, SampleSerializer


class PresampleViewSet(CollectionViewSet):
    filter_class = SampleFilter
    serializer_class = PresampleSerializer

    queryset = Sample.objects.annotate(
        latest_date=Max('data__created')
    ).prefetch_related(
        'descriptor_schema',
        'contributor'
    ).filter(
        presample=True
    ).order_by('-latest_date')


class SampleViewSet(CollectionViewSet):
    filter_class = SampleFilter
    serializer_class = SampleSerializer

    queryset = Sample.objects.annotate(
        latest_date=Max('data__modified')
    ).prefetch_related(
        'descriptor_schema',
        'contributor'
    ).filter(
        presample=False
    ).order_by('-latest_date')

    def _check_collection_permissions(self, collection_id, user):
        """Check that collection exists and user has `add` permission."""
        collection_query = Collection.objects.filter(pk=collection_id)
        if not collection_query.exists():
            raise exceptions.ValidationError('Collection id does not exist')

        collection = collection_query.first()
        if not user.has_perm('add_collection', obj=collection):
            if user.is_authenticated():
                raise exceptions.PermissionDenied()
            else:
                raise exceptions.NotFound()

    @detail_route(methods=[u'post'])
    def add_to_collection(self, request, pk=None):
        sample = self.get_object()

        if 'ids' not in request.data:
            return Response({"error": "`ids` parameter is required"}, status=status.HTTP_400_BAD_REQUEST)

        for collection_id in request.data['ids']:
            self._check_collection_permissions(collection_id, request.user)

        for collection_id in request.data['ids']:
            sample.collections.add(collection_id)

            collection = Collection.objects.get(pk=collection_id)
            for data in sample.data.all():
                collection.data.add(data)

        return Response()

    @detail_route(methods=[u'post'])
    def remove_from_collection(self, request, pk=None):
        sample = self.get_object()

        if 'ids' not in request.data:
            return Response({"error": "`ids` parameter is required"}, status=status.HTTP_400_BAD_REQUEST)

        for collection_id in request.data['ids']:
            self._check_collection_permissions(collection_id, request.user)

        for collection_id in request.data['ids']:
            sample.collections.remove(collection_id)

            collection = Collection.objects.get(pk=collection_id)
            for data in sample.data.all():
                collection.data.remove(data)

        return Response()
