import os
from unittest import mock

from django.core.management import call_command
from django.test import TestCase

from resolwe_bio.utils.test import TEST_FILES_DIR


class ImportKnowledgeBaseTestCase(TestCase):
    @mock.patch("resolwe_bio.kb.management.commands.insert_features.logger")
    def test_insert_features(self, mock_logger):
        call_command("insert_features", os.path.join(TEST_FILES_DIR, "features.tab"))
        mock_logger.info.assert_called_with(
            "Total features: 3. Inserted 3, updated 0, unchanged 0."
        )

        call_command(
            "insert_features", os.path.join(TEST_FILES_DIR, "features_update.tab.gz")
        )
        mock_logger.info.assert_called_with(
            "Total features: 4. Inserted 1, updated 1, unchanged 2."
        )

    @mock.patch("resolwe_bio.kb.management.commands.insert_mappings.logger")
    def test_insert_mappings(self, mock_logger):
        call_command(
            "insert_mappings", os.path.join(TEST_FILES_DIR, "mappings.tab.zip")
        )
        mock_logger.info.assert_called_with(
            "Total mappings: 5. Inserted 5, unchanged 0."
        )

        call_command(
            "insert_mappings", os.path.join(TEST_FILES_DIR, "mappings_update.tab")
        )
        mock_logger.info.assert_called_with(
            "Total mappings: 6. Inserted 2, unchanged 4."
        )
