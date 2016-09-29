""".. Ignore pydocstyle D400.

==============
Search Indexes
==============

"""
from haystack import indexes

from .models import Feature


class FeatureIndex(indexes.SearchIndex, indexes.Indexable):
    """Feature search index definition."""

    # This is a workaround for the Haystack limitation that all document=True fields
    # on all indices in the whole application must be of the same type. Therefore, we
    # use a CharField and have a separate MultiValueField.
    text = indexes.CharField(document=True)
    genes = indexes.MultiValueField()
    genes_auto = indexes.EdgeNgramField()

    def get_model(self):
        """Model to index."""
        return Feature

    def prepare_text(self, obj):
        """Prepare the value for the 'text' field during indexing."""
        return '\n'.join(self.prepare_genes(obj))

    def prepare_genes(self, obj):
        """Prepare the value for the 'genes' field during indexing."""
        return [obj.name, obj.feature_id] + obj.aliases

    def prepare_genes_auto(self, obj):
        """Prepare the value for the 'genes_auto' field during indexing."""
        return ' '.join(self.prepare_genes(obj))
