""".. Ignore pydocstyle D400.

==============
Search Indexes
==============

"""
from haystack import indexes

from .models import Feature


class FeatureIndex(indexes.SearchIndex, indexes.Indexable):
    """Feature search index definition."""

    text = indexes.MultiValueField(document=True)
    genes_auto = indexes.EdgeNgramField()

    def get_model(self):
        """Model to index."""
        return Feature

    def prepare_text(self, obj):
        """Prepare the value for the 'text' field during indexing."""
        return [obj.name, obj.feature_id] + obj.aliases

    def prepare_genes_auto(self, obj):
        """Prepare the value for the 'genes_auto' field during indexing."""
        return ' '.join(self.prepare_text(obj))
