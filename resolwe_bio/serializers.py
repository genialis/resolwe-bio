"""
=======================
Resolwe-Bio Serializers
=======================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from rest_framework import serializers

from resolwe.flow.serializers import CollectionSerializer


class SampleSerializer(CollectionSerializer):
    collections = serializers.PrimaryKeyRelatedField(many=True, read_only=True)

    class Meta(CollectionSerializer.Meta):
        fields = CollectionSerializer.Meta.fields + ('collections',)
