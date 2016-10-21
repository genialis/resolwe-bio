""".. Ignore pydocstyle D400.

===========
Serializers
===========

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from rest_framework import serializers

from .models import Feature


class FeatureSerializer(serializers.ModelSerializer):
    """Serializer for feature."""

    class Meta:
        """Serializer configuration."""

        model = Feature
        fields = '__all__'
