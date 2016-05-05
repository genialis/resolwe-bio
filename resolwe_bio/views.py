from __future__ import absolute_import, division, print_function, unicode_literals

from django.db.models import Max, Q

from resolwe.flow.views import CollectionViewSet, CollectionFilter
from .models import Sample


class SampleFilter(CollectionFilter):
    class Meta(CollectionFilter.Meta):
        model = Sample


class SampleViewSet(CollectionViewSet):
    filter_class = SampleFilter

    def get_queryset(self):
        queryset = Sample.objects.annotate(
            latest_date=Max('data__created')
        ).prefetch_related('descriptor_schema')

        annotated = self.request.query_params.get('annotated', None)

        if annotated == "1":
            queryset = queryset.filter(Q(descriptor__has_key='geo') & Q(descriptor__geo__has_key='annotator'))
        elif annotated == "0":
            queryset = queryset.filter(~Q(descriptor__has_key='geo') | ~Q(descriptor__geo__has_key='annotator'))

        return queryset.order_by('-latest_date')
