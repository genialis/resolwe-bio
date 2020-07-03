""".. Ignore pydocstyle D400.

===========
Serializers
===========

"""
from rest_framework import serializers

from resolwe.rest.serializers import SelectiveFieldMixin

from .models import Feature, Mapping


class FeatureSerializer(SelectiveFieldMixin, serializers.ModelSerializer):
    """Serializer for feature."""

    class Meta:
        """Serializer configuration."""

        model = Feature
        fields = [
            "aliases",
            "description",
            "feature_id",
            "full_name",
            "name",
            "source",
            "species",
            "sub_type",
            "type",
        ]


class MappingSerializer(SelectiveFieldMixin, serializers.ModelSerializer):
    """Serializer for mapping."""

    class Meta:
        """Serializer configuration."""

        model = Mapping
        fields = [
            "relation_type",
            "source_db",
            "source_id",
            "source_species",
            "target_db",
            "target_id",
            "target_species",
        ]
