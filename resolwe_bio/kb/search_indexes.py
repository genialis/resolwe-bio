""".. Ignore pydocstyle D400.

==============
Search Indexes
==============

"""
from haystack import indexes

from .models import Feature, Mapping


class FeatureIndex(indexes.SearchIndex, indexes.Indexable):
    """Feature search index definition."""

    # This is a workaround for the Haystack limitation that all document=True fields
    # on all indices in the whole application must be of the same type. Therefore, we
    # use a CharField and have a separate MultiValueField.
    text = indexes.CharField(document=True)
    genes = indexes.MultiValueField()
    source = indexes.CharField(model_attr='source')

    name_auto = indexes.EdgeNgramField(boost=10.0)
    aliases_auto = indexes.EdgeNgramField()

    def get_model(self):
        """Model to index."""
        return Feature

    def prepare_text(self, obj):
        """Prepare the value for the 'text' field during indexing."""
        return '\n'.join(self.prepare_genes(obj))

    def prepare_genes(self, obj):
        """Prepare the value for the 'genes' field during indexing."""
        return [obj.name, obj.feature_id] + obj.aliases

    def prepare_name_auto(self, obj):
        """Prepare the value for the 'name_auto' field during indexing."""
        return ' '.join([obj.name, obj.feature_id])

    def prepare_aliases_auto(self, obj):
        """Prepare the value for the 'aliases_auto' field during indexing."""
        return ' '.join(obj.aliases)


class MappingIndex(indexes.SearchIndex, indexes.Indexable):
    """Mapping search index definition."""

    text = indexes.CharField(document=True)
    # TODO: All these fields should not use the 'snowball' analyzer (Haystack limitation!).
    relation_type = indexes.CharField(model_attr='relation_type')
    source_db = indexes.CharField(model_attr='source_db')
    source_id = indexes.CharField(model_attr='source_id')
    target_db = indexes.CharField(model_attr='target_db')
    target_id = indexes.CharField(model_attr='target_id')

    def get_model(self):
        """Model to index."""
        return Mapping

    def prepare_text(self, obj):
        """Prepare the value for the 'text' field during indexing."""
        return '\n'.join([obj.source_db, obj.source_id, obj.target_db, obj.target_id])
