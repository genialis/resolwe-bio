"""Resolwe Bioinformatics configuration."""
from django.apps import AppConfig


class BaseConfig(AppConfig):
    """App configuration."""

    name = "resolwe_bio"
    verbose_name = "Resolwe Bioinformatics"
