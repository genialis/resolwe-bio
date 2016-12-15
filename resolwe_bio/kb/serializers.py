""".. Ignore pydocstyle D400.

===========
Serializers
===========

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from rest_framework import serializers
from drf_haystack.serializers import HaystackSerializerMixin

from .models import Feature, Mapping


class FeatureSerializer(serializers.ModelSerializer):
    """Serializer for feature."""

    class Meta:
        """Serializer configuration."""

        model = Feature
        fields = '__all__'


class FeatureSearchSerializer(HaystackSerializerMixin, FeatureSerializer):
    """Feature search serializer."""

    class Meta(FeatureSerializer.Meta):
        """Meta configuration for the feature search serializer."""

        search_fields = ('genes', 'source')
        field_aliases = {
            'query': 'genes',
        }


class FeatureAutocompleteSerializer(HaystackSerializerMixin, FeatureSerializer):
    """Feature autocomplete serializer."""

    class Meta(FeatureSerializer.Meta):
        """Meta configuration for the feature autocomplete serializer."""

        search_fields = ('genes_auto',)
        field_aliases = {
            'query': 'genes_auto',
        }


class MappingSerializer(serializers.ModelSerializer):
    """Serializer for mapping."""

    class Meta:
        """Serializer configuration."""

        model = Mapping
        fields = '__all__'


class MappingSearchSerializer(HaystackSerializerMixin, MappingSerializer):
    """Mapping search serializer."""

    class Meta(MappingSerializer.Meta):
        """Meta configuration for the mapping search serializer."""

        search_fields = ('relation_type', 'source_db', 'source_id', 'target_db', 'target_id')
