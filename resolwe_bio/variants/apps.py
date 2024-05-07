""".. Ignore pydocstyle D400.

===============================
Variants Base App Configuration
===============================

"""

from django.apps import AppConfig


class VariantsConfig(AppConfig):
    """App configuration."""

    name = "resolwe_bio.variants"
    label = "resolwe_bio_variants"
    verbose_name = "Resolwe Bioinformatics Variants Base"
