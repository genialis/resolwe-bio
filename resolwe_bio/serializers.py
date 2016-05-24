"""
=======================
Resolwe-Bio Serializers
=======================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from rest_framework import serializers

from resolwe.flow.serializers import CollectionSerializer
from .models import Sample


class SampleSerializer(CollectionSerializer):
    collections = serializers.PrimaryKeyRelatedField(many=True, read_only=True)

    class Meta(CollectionSerializer.Meta):
        model = Sample
        fields = CollectionSerializer.Meta.fields + ('collections',)
        read_only_fields = CollectionSerializer.Meta.read_only_fields + ('presample',)


class PresampleSerializer(CollectionSerializer):
    collections = serializers.PrimaryKeyRelatedField(many=True, read_only=True)

    class Meta(CollectionSerializer.Meta):
        model = Sample
        fields = CollectionSerializer.Meta.fields + ('collections', 'presample')
