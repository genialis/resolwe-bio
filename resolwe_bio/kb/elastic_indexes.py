"""Elastic search indexes for knowledge base."""
from __future__ import absolute_import, division, print_function, unicode_literals

import elasticsearch_dsl as dsl

from resolwe.elastic.indices import BaseDocument, BaseIndex

from .models import Feature, Mapping

# pylint: disable=invalid-name
# Analyzer for feature identifiers and names, used during boosting.
identifier_analyzer = dsl.analyzer('identifier_analyzer', tokenizer='keyword', filter=['lowercase'])
# During indexing, we lowercase terms and tokenize using edge_ngram.
autocomplete_analyzer = dsl.analyzer(
    'autocomplete_index',
    tokenizer='keyword',
    filter=[
        'lowercase',
        dsl.token_filter('autocomplete_filter', type='edgeNGram', min_gram=1, max_gram=15)
    ],
)
# During search, we only lowercase terms.
autocomplete_search_analyzer = dsl.analyzer('autocomplete_search', tokenizer='keyword', filter=['lowercase'])
# pylint: enable=invalid-name


class FeatureSearchDocument(BaseDocument):
    """Index for feature search."""

    # pylint: disable=no-member
    source = dsl.Keyword()
    feature_id = dsl.Keyword(
        # Additional subfield used for boosting during autocomplete.
        fields={
            'lower': {'type': 'text', 'analyzer': identifier_analyzer},
            'ngrams': {
                'type': 'text',
                'analyzer': autocomplete_analyzer,
                'search_analyzer': autocomplete_search_analyzer,
            },
        },
    )
    species = dsl.Keyword()
    type = dsl.Keyword()  # pylint: disable=invalid-name
    sub_type = dsl.Keyword(index=False)
    name = dsl.Keyword(
        # Additional subfield used for boosting during autocomplete.
        fields={
            'lower': {'type': 'text', 'analyzer': identifier_analyzer},
            'ngrams': {
                'type': 'text',
                'analyzer': autocomplete_analyzer,
                'search_analyzer': autocomplete_search_analyzer,
            },
        },
    )
    full_name = dsl.Text(index=False)
    description = dsl.Text(index=False)
    aliases = dsl.Keyword(
        multi=True,
        # Additional subfield used for boosting during autocomplete.
        fields={
            'ngrams': {
                'type': 'text',
                'analyzer': autocomplete_analyzer,
                'search_analyzer': autocomplete_search_analyzer,
            },
        },

    )

    class Meta:
        """Meta class for feature search document."""

        index = 'feature_search'


class FeatureSearchIndex(BaseIndex):
    """Index for feature objects used in ``FeatureSearchDocument``."""

    queryset = Feature.objects.all()
    object_type = Feature
    document_class = FeatureSearchDocument

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
