""".. Ignore pydocstyle D400.

==============================
Insert Knowledge Base Features
==============================

"""
from __future__ import absolute_import, division, print_function, unicode_literals
import csv
import logging
from tqdm import tqdm

from django.core.management.base import BaseCommand

from resolwe.elastic.builder import index_builder
from resolwe.utils import BraceMessage as __

from resolwe_bio.kb.elastic_indexes import FeatureSearchIndex
from resolwe_bio.kb.models import Feature
from .utils import decompress


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


SUBTYPE_MAP = {
    'processed_pseudogene': 'pseudo',
    'unprocessed_pseudogene': 'pseudo',
    'polymorphic_pseudogene': 'pseudo',
    'transcribed_unprocessed_pseudogene': 'pseudo',
    'unitary_pseudogene': 'pseudo',
    'transcribed_processed_pseudogene': 'pseudo',
    'transcribed_unitary_pseudogene': 'pseudo',
    'TR_J_pseudogene': 'pseudo',
    'IG_pseudogene': 'pseudo',
    'IG_D_pseudogene': 'pseudo',
    'IG_C_pseudogene': 'pseudo',
    'TR_V_pseudogene': 'pseudo',
    'IG_V_pseudogene': 'pseudo',
    'pseudogene': 'pseudo',
    'pseudo': 'pseudo',
    'asRNA': 'asRNA',
    'antisense': 'asRNA',
    'protein_coding': 'protein-coding',
    'protein-coding': 'protein-coding',
    'IG_V_gene': 'protein-coding',
    'IG_LV_gene': 'protein-coding',
    'TR_C_gene': 'protein-coding',
    'TR_V_gene': 'protein-coding',
    'TR_J_gene': 'protein-coding',
    'IG_J_gene': 'protein-coding',
    'TR_D_gene': 'protein-coding',
    'IG_C_gene': 'protein-coding',
    'IG_D_gene': 'protein-coding',
    'miRNA': 'ncRNA',
    'lincRNA': 'ncRNA',
    'processed_transcript': 'ncRNA',
    'sense_intronic': 'ncRNA',
    'sense_overlapping': 'ncRNA',
    'bidirectional_promoter_lncRNA': 'ncRNA',
    'ribozyme': 'ncRNA',
    'Mt_tRNA': 'ncRNA',
    'Mt_rRNA': 'ncRNA',
    'misc_RNA': 'ncRNA',
    'macro_lncRNA': 'ncRNA',
    '3prime_overlapping_ncRNA': 'ncRNA',
    'sRNA': 'ncRNA',
    'snRNA': 'snRNA',
    'scaRNA': 'snoRNA',
    'snoRNA': 'snoRNA',
    'rRNA': 'rRNA',
    'ncRNA': 'ncRNA',
    'tRNA': 'tRNA',
    'other': 'other',
    'unknown': 'unknown'
}


class Command(BaseCommand):
    """Insert knowledge base features."""

    help = "Insert knowledge base features"

    def add_arguments(self, parser):
        """Command arguments."""
        parser.add_argument('file_name', type=str, help="Tab-separated file with features (supports tab, gz or zip)")

    def handle(self, *args, **options):
        """Command handle."""
        count_inserted, count_updated, count_unchanged, count_failed = 0, 0, 0, 0

        for tab_file_name, line_count, tab_file in decompress(options['file_name']):
            logger.info(__("Importing features from \"{}\":", tab_file_name))

            reader = csv.DictReader(tab_file, delimiter=str('\t'))
            bar_format = '{desc}{percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'

            for row in tqdm(reader, total=line_count, bar_format=bar_format):
                aliases_text = row['Aliases'].strip()
                aliases = []
                if aliases_text and aliases_text != '-':
                    aliases = aliases_text.split(',')

                sub_type = SUBTYPE_MAP.get(row['Gene type'], 'other')

                values = {
                    'species': row['Species'],
                    'type': row['Type'],
                    'sub_type': sub_type,
                    'name': row['Name'],
                    'full_name': row['Full name'],
                    'description': row['Description'],
                    'aliases': aliases,
                }

                feature, created = Feature.objects.get_or_create(source=row['Source'],
                                                                 feature_id=row['ID'],
                                                                 defaults=values)
                if created:
                    count_inserted += 1
                else:
                    is_update = False
                    for attr, value in values.items():
                        if getattr(feature, attr) != value:
                            setattr(feature, attr, value)
                            is_update = True

                    if is_update:
                        feature.save()
                        count_updated += 1
                    else:
                        count_unchanged += 1

        index_builder.push(index=FeatureSearchIndex)

        count_total = count_inserted + count_updated + count_unchanged + count_failed
        logger.info("Total features: %d. Inserted %d, updated %d, "  # pylint: disable=logging-not-lazy
                    "unchanged %d, failed %d." %
                    (count_total, count_inserted, count_updated, count_unchanged, count_failed))
