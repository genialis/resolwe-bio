""".. Ignore pydocstyle D400.

=======
Filters
=======

"""
from functools import reduce

import django_filters as filters
from django.contrib.postgres.search import SearchQuery, SearchRank
from django.db import models
from django.db.models import F

from resolwe.flow.filters import TEXT_LOOKUPS, CheckQueryParamsMixin

from .models import Feature, Mapping


class FullTextFilter(filters.BaseCSVFilter):
    """Filter for full-text search."""

    def filter(self, queryset, values):
        """Filter field by text-search vector.

        Values, given in csv format are joined with the OR operator.

        Filtering by more than 135 values raises recursion error, so we have
        to do it in chunks.
        """
        if not values:
            return queryset

        values = [SearchQuery(v, config="simple") for v in values]

        step = 100
        result_qs = queryset.none()
        for i in range(0, len(values), step):
            query = reduce(SearchQuery.__or__, values[i : i + step])

            result_qs = result_qs.union(
                queryset.filter(**{self.field_name: query})
                # This assumes that field is already a TextSearch vector and thus
                # doesn't need to be transformed. To achieve that F function is
                # required.
                .annotate(rank=SearchRank(F(self.field_name), query))
            )

        return result_qs.order_by("-rank")


class AutoCompleteFilter(filters.Filter):
    """Filter for auto-complete search."""

    def filter(self, queryset, value):
        """Filter field by text-search vector."""
        if not value:
            return queryset.none()

        query = SearchQuery(value, config="simple")

        return (
            queryset.filter(**{self.field_name: query})
            # This assumes that field is already a TextSearch vector and thus
            # doesn't need to be transformed. To achieve that F function is
            # required.
            .annotate(rank=SearchRank(F(self.field_name), query)).order_by("-rank")
        )


class MultichoiceCharFilter(filters.BaseCSVFilter):
    """Filter by comma-separated strings."""

    def filter(self, queryset, value):
        """Perform the filtering."""
        if not value:
            return queryset

        return queryset.filter(**{"{}__in".format(self.field_name): value})


class BaseFeatureFilter(filters.FilterSet):
    """Base filter for feature endpoint."""

    class Meta:
        """Filter configuration."""

        model = Feature
        fields = {
            "source": ["exact"],
            "species": ["exact"],
            "type": ["exact"],
        }
        filter_overrides = {
            models.CharField: {
                "filter_class": MultichoiceCharFilter,
            },
        }


class FeatureFilter(CheckQueryParamsMixin, BaseFeatureFilter):
    """Filter the feature endpoint."""

    query = FullTextFilter(field_name="search")

    class Meta(BaseFeatureFilter.Meta):
        """Filter configuration."""

        fields = {
            **BaseFeatureFilter.Meta.fields,
            "feature_id": ["exact", "in"],
        }


class FeatureAutoCompleteFilter(BaseFeatureFilter):
    """Filter the auto-complete feature endpoint."""

    query = AutoCompleteFilter(field_name="search")


class MappingFilter(CheckQueryParamsMixin, filters.FilterSet):
    """Filter the mapping endpoint."""

    class Meta:
        """Filter configuration."""

        model = Mapping
        fields = {
            "source_db": TEXT_LOOKUPS[:],
            "source_id": TEXT_LOOKUPS[:],
            "source_species": TEXT_LOOKUPS[:],
            "target_db": TEXT_LOOKUPS[:],
            "target_id": TEXT_LOOKUPS[:],
            "target_species": TEXT_LOOKUPS[:],
        }
