""".. Ignore pydocstyle D400.

==============================
Insert Knowledge Base Features
==============================

"""
import csv
import json
import logging

from django.core.exceptions import ValidationError
from django.core.management.base import BaseCommand
from django.db import connection

from resolwe.utils import BraceMessage as __

from resolwe_bio.kb.models import Feature

from .utils import decompress

logger = logging.getLogger(__name__)


SUBTYPE_MAP = {
    "processed_pseudogene": "pseudo",
    "unprocessed_pseudogene": "pseudo",
    "polymorphic_pseudogene": "pseudo",
    "transcribed_unprocessed_pseudogene": "pseudo",
    "unitary_pseudogene": "pseudo",
    "transcribed_processed_pseudogene": "pseudo",
    "transcribed_unitary_pseudogene": "pseudo",
    "TR_J_pseudogene": "pseudo",
    "IG_pseudogene": "pseudo",
    "IG_D_pseudogene": "pseudo",
    "IG_C_pseudogene": "pseudo",
    "TR_V_pseudogene": "pseudo",
    "IG_V_pseudogene": "pseudo",
    "pseudogene": "pseudo",
    "pseudo": "pseudo",
    "asRNA": "asRNA",
    "antisense": "asRNA",
    "protein_coding": "protein-coding",
    "protein-coding": "protein-coding",
    "IG_V_gene": "protein-coding",
    "IG_LV_gene": "protein-coding",
    "TR_C_gene": "protein-coding",
    "TR_V_gene": "protein-coding",
    "TR_J_gene": "protein-coding",
    "IG_J_gene": "protein-coding",
    "TR_D_gene": "protein-coding",
    "IG_C_gene": "protein-coding",
    "IG_D_gene": "protein-coding",
    "miRNA": "ncRNA",
    "lincRNA": "ncRNA",
    "processed_transcript": "ncRNA",
    "sense_intronic": "ncRNA",
    "sense_overlapping": "ncRNA",
    "bidirectional_promoter_lncRNA": "ncRNA",
    "ribozyme": "ncRNA",
    "Mt_tRNA": "ncRNA",
    "Mt_rRNA": "ncRNA",
    "misc_RNA": "ncRNA",
    "macro_lncRNA": "ncRNA",
    "3prime_overlapping_ncRNA": "ncRNA",
    "sRNA": "ncRNA",
    "snRNA": "snRNA",
    "scaRNA": "snoRNA",
    "snoRNA": "snoRNA",
    "rRNA": "rRNA",
    "ncRNA": "ncRNA",
    "tRNA": "tRNA",
    "other": "other",
    "unknown": "unknown",
}


class Command(BaseCommand):
    """Insert knowledge base features."""

    help = "Insert knowledge base features"

    def add_arguments(self, parser):
        """Command arguments."""
        parser.add_argument(
            "file_name",
            type=str,
            help="Tab-separated file with features (supports tab, gz or zip)",
        )

    def handle(self, *args, **options):
        """Command handle."""
        count_total, count_inserted, count_updated = 0, 0, 0
        to_index = []

        type_choices = list(zip(*Feature.TYPE_CHOICES))[0]
        subtype_choices = list(zip(*Feature.SUBTYPE_CHOICES))[0]

        for tab_file_name, tab_file in decompress(options["file_name"]):
            logger.info(__('Importing features from "{}"...', tab_file_name))

            features = []
            unique_features = set()
            for row in csv.DictReader(tab_file, delimiter=str("\t")):
                sub_type = SUBTYPE_MAP.get(row["Gene type"], "other")

                if row["Type"] not in type_choices:
                    raise ValidationError("Unknown type: {}".format(row["Type"]))
                if sub_type not in subtype_choices:
                    raise ValidationError("Unknown subtype: {}".format(sub_type))

                aliases_text = row["Aliases"].strip()
                aliases = []
                if aliases_text and aliases_text != "-":
                    aliases = aliases_text.split(",")

                if (row["Source"], row["ID"]) in unique_features:
                    raise ValidationError(
                        "Duplicated feature (source: '{}', id: '{}') found in '{}'".format(
                            row["Source"], row["ID"], tab_file_name
                        )
                    )

                # NOTE: For performance reasons this is a list instead of a dict.
                #       Make sure that any changes also reflect in the SQL query
                #       below.
                features.append(
                    [
                        row["Source"],
                        row["ID"],
                        row["Species"],
                        row["Type"],
                        sub_type,
                        row["Name"],
                        row["Full name"],
                        row["Description"],
                        aliases,
                    ]
                )
                unique_features.add((row["Source"], row["ID"]))

            with connection.cursor() as cursor:
                cursor.execute(
                    """
                    WITH tmp AS (
                        INSERT INTO {table_name} (
                            source, feature_id, species, type,
                            sub_type, name, full_name, description,
                            aliases
                        )
                        SELECT
                            value->>0, value->>1, value->>2, value->>3,
                            value->>4, value->>5, value->>6, value->>7,
                            ARRAY(SELECT json_array_elements_text(value->8))
                        FROM json_array_elements(%s)
                        ON CONFLICT (species, source, feature_id) DO UPDATE
                        SET
                            type = EXCLUDED.type,
                            sub_type = EXCLUDED.sub_type,
                            name = EXCLUDED.name,
                            full_name = EXCLUDED.full_name,
                            description = EXCLUDED.description,
                            aliases = EXCLUDED.aliases
                        WHERE (
                            {table_name}.type, {table_name}.sub_type, {table_name}.name,
                            {table_name}.full_name, {table_name}.description, {table_name}.aliases
                        ) IS DISTINCT FROM (
                            EXCLUDED.type, EXCLUDED.sub_type, EXCLUDED.name,
                            EXCLUDED.full_name, EXCLUDED.description, EXCLUDED.aliases
                        )
                        RETURNING id, xmax
                    )
                    SELECT
                        COALESCE(array_agg(id), ARRAY[]::INTEGER[]) AS ids,
                        COUNT(CASE WHEN xmax = 0 THEN 1 END) AS count_inserted,
                        COUNT(CASE WHEN xmax != 0 THEN 1 END) AS count_updated
                    FROM tmp;
                    """.format(
                        table_name=Feature._meta.db_table,
                    ),
                    params=[json.dumps(features)],
                )
                result = cursor.fetchone()

            to_index.extend(result[0])

            count_total += len(features)
            count_inserted += result[1]
            count_updated += result[2]

        logger.info(
            "Total features: %d. Inserted %d, updated %d, unchanged %d."
            % (
                count_total,
                count_inserted,
                count_updated,
                count_total - count_inserted - count_updated,
            )
        )
