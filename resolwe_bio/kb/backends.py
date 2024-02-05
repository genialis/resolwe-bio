"""Django rest framework filter backend."""

from django_filters import rest_framework as filters


class ResolweBioFilterBackend(filters.DjangoFilterBackend):
    """Filter backend for Resolwe-Bio."""

    def get_filterset_kwargs(self, request, queryset, view):
        """Add non-form POST data to the kwargs."""
        kwargs = super().get_filterset_kwargs(request, queryset, view)
        data = kwargs["data"].copy()
        for key, value in request.data.items():
            if isinstance(value, list):
                for val in value:
                    data.appendlist(key, val)
            else:
                data.appendlist(key, value)
        kwargs["data"] = data
        return kwargs
