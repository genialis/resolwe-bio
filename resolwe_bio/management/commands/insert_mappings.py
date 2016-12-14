""".. Ignore pydocstyle D400.

==============================
Insert Knowledge Base Mappings
==============================

"""
from __future__ import absolute_import, division, print_function, unicode_literals
import csv
import logging
from tqdm import tqdm

from django.core.management.base import BaseCommand
from django.db import IntegrityError

from resolwe.utils import BraceMessage as __

from resolwe_bio.kb.models import Mapping
from .utils import decompress


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


class Command(BaseCommand):
    """Insert knowledge base mappings."""

    help = "Insert knowledge base mappings"

    def add_arguments(self, parser):
        """Command arguments."""
        parser.add_argument('file_name', type=str, help="Tab-separated file with mappings (supports tab, gz or zip)")

    def handle(self, *args, **options):
        """Command handle."""
        count_inserted, count_failed = 0, 0

        for tab_file_name, line_count, tab_file in decompress(options['file_name']):
            logger.info(__("Importing mappings from \"{}\":", tab_file_name))

            reader = csv.DictReader(tab_file, delimiter=str('\t'))
            bar_format = '{desc}{percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'

            for row in tqdm(reader, total=line_count, bar_format=bar_format):
                mapping = Mapping(relation_type=row['relation_type'],
                                  source_db=row['source_db'],
                                  source_id=row['source_id'],
                                  target_db=row['target_db'],
                                  target_id=row['target_id'])
                try:
                    mapping.save()
                    count_inserted += 1
                except IntegrityError as exc:
                    if 'duplicate key' in exc.message:
                        count_failed += 1
                    else:
                        raise

        logger.info(__("Total mappings: {}. Inserted {}, failed {}.",
                       count_inserted + count_failed, count_inserted, count_failed))
