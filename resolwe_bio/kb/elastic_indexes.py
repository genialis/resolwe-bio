"""Elastic search indexes for knowledge base."""
from __future__ import absolute_import, division, print_function, unicode_literals

import elasticsearch_dsl as dsl

from resolwe.elastic.indices import BaseDocument, BaseIndex

from .models import Feature


class FeatureSearchDocument(BaseDocument):
    """Index for feature search."""

    # pylint: disable=no-member
    source = dsl.String(index='not_analyzed')
    feature_id = dsl.String(index='not_analyzed')
    species = dsl.String()
    type = dsl.String()  # pylint: disable=invalid-name
    sub_type = dsl.String()
    name = dsl.String(index='not_analyzed')
    full_name = dsl.String()
    description = dsl.String()
    aliases = dsl.String(multi=True, index='not_analyzed')

    class Meta:
        """Meta class for feature search document."""

        index = 'feature_search'


class FeatureSearchIndex(BaseIndex):
    """Index for feature objects used in ``FeatureSearchDocument``."""

    queryset = Feature.objects.all()
    object_type = Feature
    document_class = FeatureSearchDocument
