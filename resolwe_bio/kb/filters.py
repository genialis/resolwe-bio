""".. Ignore pydocstyle D400.

=======
Filters
=======

"""
import django_filters as filters

from resolwe.flow.filters import TEXT_LOOKUPS

from .models import Mapping


class MappingFilter(filters.FilterSet):
    """Filter the mapping endpoint."""

    class Meta:
        """Filter configuration."""

        model = Mapping
        fields = {
            'source_db': TEXT_LOOKUPS[:],
            'source_id': TEXT_LOOKUPS[:],
            'source_species': TEXT_LOOKUPS[:],
            'target_db': TEXT_LOOKUPS[:],
            'target_id': TEXT_LOOKUPS[:],
            'target_species': TEXT_LOOKUPS[:],
        }
