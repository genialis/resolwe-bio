from __future__ import absolute_import, division, print_function, unicode_literals

from django.db.models import Max, Q

from rest_framework import exceptions, status
from rest_framework.decorators import detail_route
from rest_framework.response import Response

from resolwe.flow.models import Collection
from resolwe.flow.views import CollectionViewSet
from .filters import SampleFilter
from .models import Sample
from .serializers import SampleSerializer


class SampleViewSet(CollectionViewSet):
    filter_class = SampleFilter
    serializer_class = SampleSerializer

    def get_queryset(self):
        queryset = Sample.objects.annotate(
            latest_date=Max('data__created')
        ).prefetch_related('descriptor_schema')

        annotated = self.request.query_params.get('annotated', None)

        # Return annotated samples only by default
        if annotated == "0":
            queryset = queryset.filter(~Q(descriptor__has_key='geo') | ~Q(descriptor__geo__has_key='annotator'))
        else:
            queryset = queryset.filter(Q(descriptor__has_key='geo') & Q(descriptor__geo__has_key='annotator'))

        return queryset.order_by('-latest_date')

    @detail_route(methods=[u'post'])
    def add_to_collection(self, request, pk=None):
        sample = self.get_object()

        if 'ids' not in request.data:
            return Response({"error": "`ids`parameter is required"}, status=status.HTTP_400_BAD_REQUEST)

        for collection_id in request.data['ids']:
            collection_query = Collection.objects.filter(pk=collection_id)
            if not collection_query.exists():
                raise exceptions.ValidationError('Collection id does not exist')
            collection = collection_query.first()
            if not request.user.has_perm('add_collection', obj=collection):
                if request.user.is_authenticated():
                    raise exceptions.PermissionDenied()
                else:
                    raise exceptions.NotFound()

        for collection_id in request.data['ids']:
            sample.collections.add(collection_id)

        return Response()

    @detail_route(methods=[u'post'])
    def remove_from_collection(self, request, pk=None):
        sample = self.get_object()

        if 'ids' not in request.data:
            return Response({"error": "`ids`parameter is required"}, status=status.HTTP_400_BAD_REQUEST)

        for collection_id in request.data['ids']:
            collection_query = Collection.objects.filter(pk=collection_id)
            if not collection_query.exists():
                raise exceptions.ValidationError('Collection id does not exist')
            collection = collection_query.first()
            if not request.user.has_perm('add_collection', obj=collection):
                if request.user.is_authenticated():
                    raise exceptions.PermissionDenied()
                else:
                    raise exceptions.NotFound()

        for collection_id in request.data['ids']:
            sample.collections.remove(collection_id)

        return Response()
