""".. Ignore pydocstyle D400.

================================
Mutations Base App Configuration
================================

"""
from django.apps import AppConfig


class MutationsConfig(AppConfig):
    """App configuration."""

    name = "resolwe_bio.mutations"
    label = "resolwe_bio_mutations"
    verbose_name = "Resolwe Bioinformatics Mutations Base"
