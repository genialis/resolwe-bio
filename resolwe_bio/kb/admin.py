""".. Ignore pydocstyle D400.

=====
Admin
=====

"""
from django.contrib import admin

from .models import Feature, Mapping


class FeatureAdmin(admin.ModelAdmin):
    """Admin configuration for Feature model."""

    model = Feature
    search_fields = ['name']
    list_display = ('__str__', 'name', 'source', 'feature_id')


class MappingAdmin(admin.ModelAdmin):
    """Admin configuration for Mapping model."""

    model = Mapping


admin.site.register(Feature, FeatureAdmin)
admin.site.register(Mapping, MappingAdmin)
