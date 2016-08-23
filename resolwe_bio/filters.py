""".. Ignore pydocstyle D400.

===================
Resolwe Bio Filters
===================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import rest_framework_filters as filters

from resolwe.flow.filters import CollectionFilter
from resolwe.flow.models import Collection
from .models import Sample


class SampleFilter(CollectionFilter):
    """Filter the sample endpoint."""

    collection = filters.ModelChoiceFilter(queryset=Collection.objects.all())

    class Meta(CollectionFilter.Meta):
        """Filter configuration."""

        model = Sample
        fields = CollectionFilter.Meta.fields.update({'collections': ['exact', ]})
