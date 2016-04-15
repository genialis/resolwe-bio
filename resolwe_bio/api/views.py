from __future__ import absolute_import, division, print_function, unicode_literals

from django.db.models import Max, Q

from resolwe.flow.models import Collection
from resolwe.flow.views import CollectionViewSet


class SampleViewSet(CollectionViewSet):

    def get_queryset(self):
        queryset = Collection.objects.annotate(
            latest_date=Max('data__created')
        ).filter(
            descriptor_schema__slug='sample'
        ).prefetch_related('descriptor_schema')

        annotated = self.request.query_params.get('annotated', None)

        if annotated == "1":
            queryset = queryset.filter(Q(descriptor__has_key='geo') & Q(descriptor__geo__has_key='annotator'))
        elif annotated == "0":
            queryset = queryset.filter(~Q(descriptor__has_key='geo') | ~Q(descriptor__geo__has_key='annotator'))

        return queryset.order_by('-latest_date')
