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
        fields = '__all__'


class MappingSerializer(SelectiveFieldMixin, serializers.ModelSerializer):
    """Serializer for mapping."""

    class Meta:
        """Serializer configuration."""

        model = Mapping
        fields = '__all__'
