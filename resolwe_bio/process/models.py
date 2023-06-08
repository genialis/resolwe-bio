"""Additional models used by ProcessBio."""
from typing import Any

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

    def get_features(
        self,
        source,
        species,
        feature_ids,
        required_fields=(
            "source",
            "feature_id",
            "species",
            "type",
            "sub_type",
            "name",
            "full_name",
            "description",
            "aliases",
        ),
    ) -> list:
        """Get features from the database."""
        return communicator.filter_features(
            (source, species, feature_ids, required_fields)
        )


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
