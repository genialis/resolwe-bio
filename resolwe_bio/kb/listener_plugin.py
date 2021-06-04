"""Register knowledgebase tables in listener."""
from django.contrib.auth import get_user_model
from django.db.models import QuerySet

from resolwe.flow.managers.listener.python_process_plugin import ExposeObjectPlugin

UserClass = get_user_model()


class ExposeFeature(ExposeObjectPlugin):
    """Expose the Feature knowledge base model in listener."""

    full_model_name = "resolwe_bio_kb.Feature"

    def filter_objects(self, user: UserClass, queryset: QuerySet, data) -> QuerySet:
        """Filter the objects for the given user."""
        return queryset


class ExposeMapping(ExposeObjectPlugin):
    """Expose the Mapping knowledge base model in listener."""

    full_model_name = "resolwe_bio_kb.Mapping"

    def filter_objects(self, user: UserClass, queryset: QuerySet, data) -> QuerySet:
        """Filter the objects for the given user."""
        return queryset
