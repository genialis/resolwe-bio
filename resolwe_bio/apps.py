from django.apps import AppConfig


class BaseConfig(AppConfig):
    name = 'resolwe_bio'
    verbose_name = 'Resolwe Bioinformatics'

    def ready(self):
        """
        Performs application initialization.
        """

        # Register signals handlers
        from . import signals  # pylint: disable=unused-import
