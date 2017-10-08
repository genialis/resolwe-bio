""".. Ignore pydocstyle D400.

=======
Filters
=======

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import rest_framework_filters as filters

from .models import Mapping


class MappingFilter(filters.FilterSet):
    """Filter the mapping endpoint."""

    source_db = filters.AllLookupsFilter()
    source_id = filters.AllLookupsFilter()
    target_db = filters.AllLookupsFilter()
    target_id = filters.AllLookupsFilter()

    class Meta:
        """Filter configuration."""

        model = Mapping
        fields = [
            'source_db', 'source_id', 'target_db', 'target_id'
        ]
