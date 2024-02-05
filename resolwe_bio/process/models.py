"""Additional models used by ProcessBio."""

from typing import Any, Dict, List

from resolwe.process.communicator import communicator
from resolwe.process.models import Model


class Feature(Model):
    """Expose Feature model in Python Processes."""

    _app_name = "resolwe_bio_kb"
    _model_name = "Feature"
    _filter_response_fields = [
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

    @classmethod
    def filter(cls, **filters: Dict[str, Any]) -> List["Model"]:
        """Filter features from the database."""
        # Make sure attributes have 'id' in the first place.
        attributes = filters.pop("__fields", None)
        attributes = attributes or cls._filter_response_fields
        attributes = ["id"] + [
            attribute for attribute in attributes if attribute != "id"
        ]

        if set(filters.keys()) == {"source", "species", "feature_id__in"}:
            objects = communicator.filter_features(
                (
                    filters["source"],
                    filters["species"],
                    filters["feature_id__in"],
                    attributes,
                )
            )
        else:
            objects = communicator.filter_objects(
                cls._app_name, cls._model_name, filters, attributes
            )
        models = []
        for entry in objects:
            model = cls(entry[0])
            for field_name, value in zip(attributes[1:], entry[1:]):
                field = model.fields[field_name]
                model._cache[field_name] = field.clean(value)
            models.append(model)
        return models


class Mapping(Model):
    """Expose Feature model in Python Processes."""

    _app_name = "resolwe_bio_kb"
    _model_name = "Mapping"
    _filter_response_fields = [
        "relation_type",
        "source_db",
        "source_id",
        "source_species",
        "target_db",
        "target_id",
        "target_species",
    ]
