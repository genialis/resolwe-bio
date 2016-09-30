""".. Ignore pydocstyle D400.

=======================
Resolwe-Bio Serializers
=======================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from rest_framework import serializers

from resolwe.flow.serializers import CollectionSerializer
from .models import Sample


class SampleSerializer(CollectionSerializer):
    """Serializer for sample."""

    collections = serializers.PrimaryKeyRelatedField(many=True, read_only=True)

    def create(self, validated_data):
        """Create ``Sample``and set ``presample``to ``False``."""
        validated_data['presample'] = False

        return super(SampleSerializer, self).create(validated_data)

    class Meta(CollectionSerializer.Meta):
        """Serializer configuration."""

        model = Sample
        fields = CollectionSerializer.Meta.fields + ('collections',)
        read_only_fields = CollectionSerializer.Meta.read_only_fields + ('presample',)


class PresampleSerializer(CollectionSerializer):
    """Serializer for presample."""

    collections = serializers.PrimaryKeyRelatedField(many=True, read_only=True)

    class Meta(CollectionSerializer.Meta):
        """Serializer configuration."""

        model = Sample
        fields = CollectionSerializer.Meta.fields + ('collections', 'presample')
