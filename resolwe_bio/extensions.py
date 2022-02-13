"""Resolwe-bio extensions of Resolwe core."""
from django.db.models import Q
from django_filters import rest_framework as filters

from resolwe.composer import composer


class ExtendedCollectionFilter:
    """Collection filter extension."""

    descriptor__general__species = filters.CharFilter(
        field_name="entity__descriptor__general__species",
        lookup_expr="icontains",
        distinct=True,
    )
    tissue_type = filters.CharFilter(method="_filter_by_tissue")
    treatment = filters.CharFilter(method="_filter_by_treatment")
    outcome_is_defined = filters.CharFilter(method="_filter_by_outcome")

    def _filter_by_multiple_fields(self, queryset, fields, value):
        """Filter by any of the given fields.

        Return exactly those objects, who icontain `value` in at least one of
        `fields`.
        """

        query = Q()
        for field in fields:
            query |= Q(**{field + "__icontains": value})

        return queryset.filter(query).distinct()

    def _filter_by_tissue(self, queryset, name, value):
        locations = [
            "entity__descriptor__general__organ",
            "entity__descriptor__general__biosample_source",
            "entity__descriptor__disease_information__organ_part",
            "entity__descriptor__disease_information__biopsy_site",
            "entity__descriptor__pathological_information__organ_part",
            "entity__descriptor__pathological_information__biopsy_site",
        ]
        return self._filter_by_multiple_fields(queryset, locations, value)

    def _filter_by_treatment(self, queryset, name, value):
        locations = [
            "entity__descriptor__general__biosample_treatment",
            "entity__descriptor__treatment_type__drug",
            "entity__descriptor__immuno_oncology_treatment_type__io_drug",
        ]
        return self._filter_by_multiple_fields(queryset, locations, value)

    def _filter_by_outcome(self, queryset, name, value):
        locations = [
            "entity__descriptor__response_and_survival_analysis__clinical_benefit",
            "entity__descriptor__response_and_survival_analysis__confirmed_bor",
            "entity__descriptor__response_and_survival_analysis__unconfirmed_bor",
            "entity__descriptor__response_and_survival_analysis__pfs",
            "entity__descriptor__response_and_survival_analysis__os",
            "entity__descriptor__response_and_survival_analysis__dfs",
            "entity__descriptor__response_and_survival_analysis__ttp",
        ]
        return self._filter_by_multiple_fields(queryset, locations, value)


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


composer.add_extension(
    "resolwe.flow.filters.CollectionFilter", ExtendedCollectionFilter
)
composer.add_extension("resolwe.flow.filters.DataFilter", ExtendedDataFilter)
composer.add_extension("resolwe.flow.filters.EntityFilter", ExtendedEntityFilter)
