""".. Ignore pydocstyle D400.

================================
Knowledge Base App Configuration
================================

"""
from django.apps import AppConfig


class KnowledgeBaseConfig(AppConfig):
    """App configuration."""

    name = 'resolwe_bio.kb'
    label = 'resolwe_bio_kb'
    verbose_name = 'Resolwe Bioinformatics Knowledge Base'

    def ready(self):
        """Perform application initialization."""
        from haystack import connections

        for connection in connections.all():
            if connection.__class__.__name__ == 'ElasticsearchSearchEngine':
                # Modify elastic search schema. The default mapping for edge_ngram is incorrect
                # as it also uses the edgengram_analyzer during querying.
                from haystack.backends.elasticsearch_backend import FIELD_MAPPINGS
                FIELD_MAPPINGS['edge_ngram'] = {
                    'type': 'string',
                    'analyzer': 'edgengram_analyzer',
                    'search_analyzer': 'standard',
                }
                break
