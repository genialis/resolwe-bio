"""Resolwe Bioinformatics configuration."""
from django.apps import AppConfig


class BaseConfig(AppConfig):
    """App configuration."""

    name = 'resolwe_bio'
    verbose_name = 'Resolwe Bioinformatics'

    def ready(self):
        """Perform application initialization."""
        # Register signals handlers
        from . import signals  # pylint: disable=unused-variable
