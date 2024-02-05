""".. Ignore pydocstyle D400.

===================
Resolwe Bio Filters
===================

"""

import django_filters as filters
from rest_framework import exceptions

from resolwe.flow.filters import CollectionFilter, DataFilter, EntityFilter

from resolwe_bio.models import Sample


class BioCollectionFilter(CollectionFilter):
    """Filter the collection endpoint.

    Enable filtering collections by the entity.

    .. IMPORTANT::

        :class:`CollectionViewSet` must be patched before using it in
        urls to enable this feature:

            .. code:: python

                CollectionViewSet.filterset_class = BioCollectionFilter

    """

    sample = filters.ModelChoiceFilter(
        field_name="entity", queryset=Sample.objects.all()
    )


class BioEntityFilter(EntityFilter):
    """Filter the entity endpoint.

    Enable filtering collections by the entity.

    .. IMPORTANT::

        :class:`EntityViewSet` must be patched before using it in
        urls to enable this feature:

            .. code:: python

                EntityViewSet.filterset_class = BioEntityFilter

    """

    descriptor__subject_information__sample_label__icontains = filters.CharFilter(
        field_name="descriptor__subject_information__sample_label",
        lookup_expr="icontains",
    )
    descriptor__subject_information__subject_id__icontains = filters.CharFilter(
        field_name="descriptor__subject_information__subject_id",
        lookup_expr="icontains",
    )
    descriptor__subject_information__batch__exact = filters.CharFilter(
        field_name="descriptor__subject_information__batch",
        method="filter_exact_number",
    )
    descriptor__subject_information__group__iexact = filters.CharFilter(
        field_name="descriptor__subject_information__group", lookup_expr="iexact"
    )
    descriptor__disease_information__disease_type__icontains = filters.CharFilter(
        field_name="descriptor__disease_information__disease_type",
        lookup_expr="icontains",
    )
    descriptor__disease_information__disease_status__iexact = filters.CharFilter(
        field_name="descriptor__disease_information__disease_status",
        lookup_expr="iexact",
    )
    descriptor__immuno_oncology_treatment_type__io_drug__iexact = filters.CharFilter(
        field_name="descriptor__immuno_oncology_treatment_type__io_drug",
        lookup_expr="iexact",
    )
    descriptor__immuno_oncology_treatment_type__io_treatment__iexact = (
        filters.CharFilter(
            field_name="descriptor__immuno_oncology_treatment_type__io_treatment",
            lookup_expr="iexact",
        )
    )
    descriptor__response_and_survival_analysis__confirmed_bor__iexact = (
        filters.CharFilter(
            field_name="descriptor__response_and_survival_analysis__confirmed_bor",
            lookup_expr="iexact",
        )
    )
    descriptor__response_and_survival_analysis__pfs_event__iexact = filters.CharFilter(
        field_name="descriptor__response_and_survival_analysis__pfs_event",
        lookup_expr="iexact",
    )
    descriptor__general__description__icontains = filters.CharFilter(
        field_name="descriptor__general__description", lookup_expr="icontains"
    )
    descriptor__general__biosample_source__icontains = filters.CharFilter(
        field_name="descriptor__general__biosample_source", lookup_expr="icontains"
    )
    descriptor__general__biosample_treatment__icontains = filters.CharFilter(
        field_name="descriptor__general__biosample_treatment", lookup_expr="icontains"
    )

    def filter_exact_number(self, queryset, name, value):
        """Transform value into an integer and filter by exact value."""
        try:
            value = int(value)
        except ValueError:
            raise exceptions.ParseError(f"Value of attribute {name} must be a number.")

        return queryset.filter(**{name: value})


class BioDataFilter(DataFilter):
    """Filter the data endpoint.

    Enable filtering data by the sample.

    .. IMPORTANT::

        :class:`DataViewSet` must be patched before using it in urls to
        enable this feature:

            .. code:: python

                DataViewSet.filterset_class = BioDataFilter

    """

    sample = filters.ModelChoiceFilter(
        field_name="entity", queryset=Sample.objects.all()
    )
