"""Resolwe Bioinformatics configuration."""
from django.apps import AppConfig

from resolwe.composer import composer
from resolwe_bio.extensions import ExtendedDataFilter


class BaseConfig(AppConfig):
    """App configuration."""

    name = 'resolwe_bio'
    verbose_name = 'Resolwe Bioinformatics'

    composer.add_extension("resolwe.flow.filters.DataFilter", ExtendedDataFilter)
