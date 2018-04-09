"""Resolwe-bio extensions of Resolwe core."""
import elasticsearch_dsl as dsl
from elasticsearch_dsl.query import Q

from resolwe.elastic.composer import composer


class ExtendedDataDocument(object):
    """Data ES document extensions."""

    source = dsl.Keyword()
    species = dsl.Text()
    build = dsl.Keyword()
    feature_type = dsl.Keyword()


class ExtendedDataIndex(object):
    """Data ES index extensions."""

    mapping = {
        'source': 'output.source',
        'species': 'output.species',
        'build': 'output.build',
        'feature_type': 'output.feature_type',
    }


class ExtendedDataViewSet(object):
    """Data viewset extensions."""

    filtering_fields = ('source', 'species', 'build', 'feature_type')

    @staticmethod
    def text_filter(value):
        """Extend full-text data filter."""
        return [
            Q('match', species={'query': value, 'operator': 'and', 'boost': 2.0}),
            Q('match', source={'query': value, 'operator': 'and', 'boost': 2.0}),
            Q('match', build={'query': value, 'operator': 'and', 'boost': 2.0}),
            Q('match', feature_type={'query': value, 'operator': 'and', 'boost': 1.0}),
        ]


composer.add_extension('resolwe.flow.elastic_indexes.data.DataDocument', ExtendedDataDocument)
composer.add_extension('resolwe.flow.elastic_indexes.data.DataIndex', ExtendedDataIndex)
composer.add_extension('resolwe.flow.views.data.DataViewSet', ExtendedDataViewSet)
