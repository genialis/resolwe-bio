""".. Ignore pydocstyle D400.

================================
Knowledge Base App Configuration
================================

"""
from django.apps import AppConfig


class KnowledgeBaseConfig(AppConfig):
    """App configuration."""

    name = "resolwe_bio.kb"
    label = "resolwe_bio_kb"
    verbose_name = "Resolwe Bioinformatics Knowledge Base"
