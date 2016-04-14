from django.db.models import Max, Q

from rest_framework.decorators import list_route
from rest_framework.response import Response

from resolwe.flow.models import Collection
from resolwe.flow.views import CollectionViewSet


class SampleViewSet(CollectionViewSet):
    queryset = Collection.objects.filter(descriptor_schema__slug='sample').prefetch_related('descriptor_schema')

    @list_route(methods=[u'get'])
    def annotated(self, request):
        """Return list of annotated `Samples`."""
        queryset = self.get_queryset().annotate(
            latest_date=Max('data__created')
        ).filter(
            Q(descriptor__has_key='geo') & Q(descriptor__geo__has_key='annotator')
        ).order_by('-latest_date')
        serializer = self.serializer_class(queryset, many=True, context={'request': request})

        return Response(serializer.data)

    @list_route(methods=[u'get'])
    def unannotated(self, request):
        """Return list of unannotated `Samples`."""
        queryset = self.get_queryset().annotate(
            latest_date=Max('data__created')
        ).filter(
            ~Q(descriptor__has_key='geo') | ~Q(descriptor__geo__has_key='annotator')
        ).order_by('-latest_date')
        serializer = self.serializer_class(queryset, many=True, context={'request': request})

        return Response(serializer.data)
