"""Elastic search indexes for knowledge base."""
from __future__ import absolute_import, division, print_function, unicode_literals

import elasticsearch_dsl as dsl

from resolwe.elastic.indices import BaseDocument, BaseIndex

from .models import Feature, Mapping

# Analyzer for feature identifiers and names, used during boosting.
# pylint: disable=invalid-name
identifier_analyzer = dsl.analyzer('identifier_analyzer', tokenizer='keyword', filter=['lowercase'])
# pylint: enable=invalid-name


class FeatureSearchDocument(BaseDocument):
    """Index for feature search."""

    # pylint: disable=no-member
    source = dsl.Keyword()
    feature_id = dsl.Keyword(
        # Additional subfield used for boosting during autocomplete.
        fields={'lower': {'type': 'text', 'analyzer': identifier_analyzer}},
    )
    species = dsl.Keyword()
    type = dsl.Keyword(index=False)  # pylint: disable=invalid-name
    sub_type = dsl.Keyword(index=False)
    name = dsl.Keyword(
        # Additional subfield used for boosting during autocomplete.
        fields={'lower': {'type': 'text', 'analyzer': identifier_analyzer}},
    )
    full_name = dsl.Text(index=False)
    description = dsl.Text(index=False)
    aliases = dsl.Keyword(
        multi=True,
        # Additional subfield used for boosting during autocomplete.
        fields={'lower': {'type': 'text', 'analyzer': identifier_analyzer}},
    )

    # Autocomplete.
    autocomplete = dsl.Text(
        multi=True,
        # During indexing, we lowercase terms and tokenize using edge_ngram.
        analyzer=dsl.analyzer(
            'autocomplete_index',
            tokenizer='keyword',
            filter=[
                'lowercase',
                dsl.token_filter(
                    'autocomplete_filter',
                    type='edgeNGram',
                    min_gram=1,
                    max_gram=15
                )
            ],
        ),
        # During search, we only lowercase terms.
        search_analyzer=dsl.analyzer(
            'autocomplete_search',
            tokenizer='keyword',
            filter=[
                'lowercase'
            ],
        ),
    )

    class Meta:
        """Meta class for feature search document."""

        index = 'feature_search'


class FeatureSearchIndex(BaseIndex):
    """Index for feature objects used in ``FeatureSearchDocument``."""

    queryset = Feature.objects.all()
    object_type = Feature
    document_class = FeatureSearchDocument

    def get_autocomplete_value(self, obj):
        """Return autocomplete value."""
        return [
            obj.feature_id,
            obj.name
        ] + obj.aliases

    def get_permissions(self, obj):
        """Skip since Feature objects have no permissions."""
        return {
            'users': [],
            'groups': [],
            'public': False,
        }


class MappingSearchDocument(BaseDocument):
    """Index for mapping search."""

    # pylint: disable=no-member
    relation_type = dsl.Keyword(index=False)
    source_db = dsl.Keyword()
    source_id = dsl.Keyword()
    source_species = dsl.Keyword()
    target_db = dsl.Keyword()
    target_id = dsl.Keyword()
    target_species = dsl.Keyword()
    relation_type = dsl.Keyword()

    class Meta:
        """Meta class for mapping search document."""

        index = 'mapping_search'


class MappingSearchIndex(BaseIndex):
    """Index for mapping objects used in ``MappingSearchDocument``."""

    queryset = Mapping.objects.all()
    object_type = Mapping
    document_class = MappingSearchDocument

    def get_permissions(self, obj):
        """Skip since Mapping objects have no permissions."""
        return {
            'users': [],
            'groups': [],
            'public': False,
        }
