""".. Ignore pydocstyle D400.

===================
Resolwe Bio Filters
===================

"""
import django_filters as filters

from resolwe.flow.filters import CollectionFilter, DataFilter
from resolwe_bio.models import Sample


class BioCollectionFilter(CollectionFilter):
    """Filter the collection endpoint.

    Enable filtering collections by the entity.

    .. IMPORTANT::

        :class:`CollectionViewSet` must be patched before using it in
        urls to enable this feature:

            .. code:: python

                CollectionViewSet.filter_class = BioCollectionFilter

    """

    sample = filters.ModelChoiceFilter(field_name='entity', queryset=Sample.objects.all())


class BioDataFilter(DataFilter):
    """Filter the data endpoint.

    Enable filtering data by the sample.

    .. IMPORTANT::

        :class:`DataViewSet` must be patched before using it in urls to
        enable this feature:

            .. code:: python

                DataViewSet.filter_class = BioDataFilter

    """

    sample = filters.ModelChoiceFilter(field_name='entity', queryset=Sample.objects.all())
