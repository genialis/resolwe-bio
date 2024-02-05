"""Register knowledgebase tables in listener."""

from typing import Tuple

from django.contrib.auth import get_user_model
from django.db.models import QuerySet

from resolwe.flow.executors.socket_utils import Message, Response
from resolwe.flow.managers.listener.listener import Processor
from resolwe.flow.managers.listener.permission_plugin import ExposeObjectPlugin
from resolwe.flow.managers.listener.plugin import (
    ListenerPlugin,
    listener_plugin_manager,
)
from resolwe.flow.managers.listener.redis_cache import CachedObjectPlugin, cache_manager

from resolwe_bio.kb.models import Feature

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


class FeatureCache(CachedObjectPlugin):
    """Cache the Feature knowledge base model in Redis."""

    model = Feature
    cached_fields = [
        "id",
        "source",
        "feature_id",
        "species",
        "type",
        "sub_type",
        "name",
        "full_name",
        "description",
        "aliases",
    ]
    identifier_fields = ("feature_id", "source", "species")
    expiration_time = None  # Never expire.


class KnowledgeBasePlugin(ListenerPlugin):
    """Handler methods for KnowledgeBase methods."""

    name = "KnowledgeBase plugin"
    plugin_manager = listener_plugin_manager

    def handle_filter_features(
        self,
        data_id: int,
        message: Message[Tuple[str, str, list[str], list[str]]],
        manager: "Processor",
    ) -> Response[dict]:
        """Get the feature data based on the given identifiers."""
        source, species, feature_ids, requested_fields = message.message_data
        identifiers_list = list()
        for feature_id in feature_ids:
            # Identifiers must be ordered same as FeatureCache identifier_fields.
            data = {"feature_id": feature_id, "source": source, "species": species}
            identifiers = (data.get(name) for name in FeatureCache.identifier_fields)
            identifiers_list.append(tuple(identifiers))

        # The get itself will wait for up to one minute when lock is set to any entry.
        cached_data = cache_manager.mget(Feature, identifiers_list)
        to_return = [
            [entry[field] for field in requested_fields]
            for entry in cached_data
            if entry is not None
        ]

        missing_identifiers = []
        for position, (identifier, value) in enumerate(
            zip(identifiers_list, cached_data)
        ):
            if value is None:
                missing_identifiers.append(identifier)

        # Read the missing data from the database.
        if missing_identifiers:
            try:
                cache_manager.lock(Feature, missing_identifiers)
                feature_id_index = FeatureCache.identifier_fields.index("feature_id")
                missing_feature_ids = [
                    identifier[feature_id_index] for identifier in missing_identifiers
                ]
                missing_features = Feature.objects.filter(
                    source=source, species=species, feature_id__in=missing_feature_ids
                )
                # Cache the missing data.
                cache_manager.mcache(missing_features)
                # Fill the missing values to the to_return list.
                to_return.extend(missing_features.values_list(*requested_fields))
            finally:
                # Unlock the entries.
                cache_manager.unlock(Feature, missing_identifiers)

        return message.respond(to_return)
