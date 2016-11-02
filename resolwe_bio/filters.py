""".. Ignore pydocstyle D400.

===================
Resolwe Bio Filters
===================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import rest_framework_filters as filters

from resolwe.flow.filters import CollectionFilter, DataFilter
from resolwe.flow.models import Collection
from .models import Sample


class SampleFilter(CollectionFilter):
    """Filter the sample endpoint."""

    collection = filters.ModelChoiceFilter(queryset=Collection.objects.all())

    class Meta(CollectionFilter.Meta):
        """Filter configuration."""

        model = Sample
        fields = CollectionFilter.Meta.fields.update({'collections': ['exact', ]})


class BioCollectionFilter(CollectionFilter):
    """Filter the collection endpoint.

    Enable filtering collections by the sample.

    .. IMPORTANT::

        :class:`CollectionViewSet` must be patched before using it in
        urls to enable this feature:

            .. code:: python

                CollectionViewSet.filter_class = BioCollectionFilter

    """

    sample = filters.ModelChoiceFilter(queryset=Sample.objects.all())

    class Meta(CollectionFilter.Meta):
        """Filter configuration."""

        fields = CollectionFilter.Meta.fields.update({'sample': ['exact', ]})


class BioDataFilter(DataFilter):
    """Filter the data endpoint.

    Enable filtering data by the sample.

    .. IMPORTANT::

        :class:`DataViewSet` must be patched before using it in urls to
        enable this feature:

            .. code:: python

                DataViewSet.filter_class = BioDataFilter

    """

    sample = filters.ModelChoiceFilter(queryset=Sample.objects.all())

    class Meta(DataFilter.Meta):
        """Filter configuration."""

        fields = DataFilter.Meta.fields.update({'sample': ['exact', ]})
