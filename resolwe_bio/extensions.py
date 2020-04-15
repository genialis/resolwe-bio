"""Resolwe-bio extensions of Resolwe core."""
from django_filters import rest_framework as filters

from resolwe.composer import composer


class ExtendedDataFilter:
    """Data filter extensions."""

    def filter_output(self, queryset, name, value):
        """Filter queryset by genome build."""
        return queryset.filter(**{"output__{}".format(name): value})

    def filter_output_icontains(self, queryset, name, value):
        """Filter queryset by genome build."""
        return queryset.filter(**{"output__{}__icontains".format(name): value})

    # These filters use custom indexes defined in migrations.
    build = filters.CharFilter(method="filter_output")
    feature_type = filters.CharFilter(method="filter_output")
    source = filters.CharFilter(method="filter_output")
    species = filters.CharFilter(method="filter_output_icontains")


class ExtendedEntityFilter:
    """Data filter extensions."""

    def filter_species(self, queryset, name, value):
        """Filter queryset by genome build."""
        return queryset.filter(descriptor__general__species__icontains=value)

    species = filters.CharFilter(method="filter_species")


composer.add_extension("resolwe.flow.filters.DataFilter", ExtendedDataFilter)
composer.add_extension("resolwe.flow.filters.EntityFilter", ExtendedEntityFilter)
